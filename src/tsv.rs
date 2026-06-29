use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::ffi::CString;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::os::raw::{c_char, c_double};

unsafe extern "C" {
    fn strtod(nptr: *const c_char, endptr: *mut *mut c_char) -> c_double;
}

use anyhow::{bail, Result};

pub const CHOICES_MAX: usize = 7;
pub const AGGR_SEP: char = ',';

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DateFormat {
    Year,
    Ymd,
    None,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Date {
    pub year: i16,
    pub month: i8,
    pub day: i8,
}

#[derive(Debug, Clone, PartialEq)]
pub enum JsonValue {
    Null,
    String(String),
    Int(i64),
    Double {
        n: f64,
        decimals: usize,
        scientific: bool,
    },
    Boolean(bool),
    Array(JsonArray),
    Map(JsonMap),
}

#[derive(Debug, Clone, PartialEq)]
pub struct JsonArray {
    pub values: Vec<JsonValue>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct JsonMap {
    pub values: BTreeMap<String, JsonValue>,
}

#[derive(Debug, Clone, PartialEq)]
enum JsonTokenType {
    Empty,
    Name,
    Delimiter,
    Text,
    Integer,
    Double,
}

#[derive(Debug, Clone)]
struct JsonToken {
    token_type: JsonTokenType,
    name: String,
    n: i64,
    d: f64,
    decimals: usize,
    scientific: bool,
}

struct JsonCharInput {
    text: Vec<char>,
    pos: usize,
}

impl JsonCharInput {
    fn new(text: &str) -> Self {
        JsonCharInput {
            text: text.chars().collect(),
            pos: 0,
        }
    }

    fn get(&mut self) -> Option<char> {
        let c = self.text.get(self.pos).copied();
        if c.is_some() {
            self.pos += 1;
        }
        c
    }

    fn unget(&mut self) {
        assert!(self.pos > 0);
        self.pos -= 1;
    }

    fn error<T>(&self, what: &str) -> Result<T> {
        bail!("Json parse error at character {}: {}", self.pos, what)
    }
}

impl JsonToken {
    fn new() -> Self {
        JsonToken {
            token_type: JsonTokenType::Empty,
            name: String::new(),
            n: 0,
            d: 0.0,
            decimals: 0,
            scientific: false,
        }
    }

    fn read_input(in_arg: &mut JsonCharInput) -> Result<Self> {
        let mut token = JsonToken::new();
        let mut c = loop {
            let Some(c) = in_arg.get() else {
                return Ok(token);
            };
            if !c.is_whitespace() {
                break c;
            }
        };

        if c == '\'' || c == '"' {
            token.token_type = JsonTokenType::Text;
            let quote = c;
            loop {
                let Some(c_text) = in_arg.get() else {
                    return in_arg.error("Text is not finished: end of file");
                };
                if c_text == '\n' {
                    continue;
                }
                if c_text == quote {
                    break;
                }
                token.name.push(c_text);
            }
        } else if c.is_ascii_digit() || c == '-' {
            let mut minus_possible = true;
            loop {
                if c.is_ascii_digit()
                    || c == '.'
                    || c == 'e'
                    || c == 'E'
                    || (c == '-' && minus_possible)
                    || c == 'x'
                    || c.is_ascii_hexdigit()
                {
                    token.name.push(c);
                } else {
                    in_arg.unget();
                    break;
                }
                let Some(next) = in_arg.get() else {
                    break;
                };
                minus_possible = next == 'e' || next == 'E';
                c = next;
            }
            if token.name == "-" {
                token.token_type = JsonTokenType::Delimiter;
            } else {
                token.token_type = JsonTokenType::Text;
                token.to_number_date();
            }
        } else if c.is_ascii_alphabetic() || c == '_' {
            token.token_type = JsonTokenType::Name;
            loop {
                if c.is_ascii_alphabetic() || c.is_ascii_digit() || c == '_' {
                    token.name.push(c);
                } else {
                    in_arg.unget();
                    break;
                }
                let Some(next) = in_arg.get() else {
                    break;
                };
                c = next;
            }
        } else if !c.is_whitespace()
            && c.is_ascii_graphic()
            && !(c.is_ascii_alphanumeric() || c == '_')
        {
            token.token_type = JsonTokenType::Delimiter;
            token.name.push(c);
        } else {
            return in_arg.error(&format!("Unknown token starting with ASCII {}", c as u32));
        }

        Ok(token)
    }

    fn is_delimiter(&self, c: char) -> bool {
        self.token_type == JsonTokenType::Delimiter && self.name.starts_with(c)
    }

    fn to_number_date(&mut self) {
        if self.token_type != JsonTokenType::Name && self.token_type != JsonTokenType::Text {
            return;
        }
        if self.name.contains(' ') {
            return;
        }
        if let Some(hex) = self.name.strip_prefix("0x") {
            if let Ok(n) = i64::from_str_radix(hex, 16) {
                self.n = n;
                self.token_type = JsonTokenType::Integer;
            }
        } else if self.name.contains('.') || self.name.contains('e') || self.name.contains('E') {
            if let Ok(d) = self.name.parse::<f64>() {
                self.d = d;
                self.token_type = JsonTokenType::Double;
                if let Some(point) = self.name.find('.') {
                    self.decimals = self.name[point + 1..]
                        .chars()
                        .take_while(|c| c.is_ascii_digit())
                        .count();
                }
                self.scientific = self.name.contains('e');
            }
        } else if let Ok(n) = self.name.parse::<i64>() {
            self.n = n;
            self.token_type = JsonTokenType::Integer;
        }
    }
}

impl JsonArray {
    pub fn new(_name: Option<&str>) -> Self {
        JsonArray { values: Vec::new() }
    }

    fn parse(in_arg: &mut JsonCharInput) -> Result<Self> {
        let mut array = JsonArray::new(None);
        let mut first = true;
        loop {
            let mut token = JsonToken::read_input(in_arg)?;
            if token.is_delimiter(']') {
                break;
            }
            if !first {
                if !token.is_delimiter(',') {
                    return in_arg.error("','");
                }
                token = JsonToken::read_input(in_arg)?;
            }
            array.values.push(JsonValue::parse(in_arg, token)?);
            first = false;
        }
        Ok(array)
    }

    pub fn add_null(&mut self) {
        self.values.push(JsonValue::Null);
    }

    pub fn add_int(&mut self, value: i64) {
        self.values.push(JsonValue::Int(value));
    }

    pub fn add_double(&mut self, value: f64, decimals: usize) {
        self.values.push(JsonValue::Double {
            n: value,
            decimals,
            scientific: false,
        });
    }

    pub fn add_string(&mut self, value: &str) {
        self.values.push(JsonValue::String(value.to_string()));
    }

    pub fn add_boolean(&mut self, value: bool) {
        self.values.push(JsonValue::Boolean(value));
    }

    pub fn add_array(&mut self, value: JsonArray) {
        self.values.push(JsonValue::Array(value));
    }

    pub fn add_map(&mut self, value: JsonMap) {
        self.values.push(JsonValue::Map(value));
    }

    pub fn at(&self, index: usize) -> Option<&JsonValue> {
        Some(
            self.values
                .get(index)
                .unwrap_or_else(|| panic!("Index out of range: {}", index)),
        )
    }

    pub fn get_size(&self) -> usize {
        self.values.len()
    }

    pub fn to_json_text(&self) -> String {
        let fields = self
            .values
            .iter()
            .map(JsonValue::to_json_text)
            .collect::<Vec<_>>();
        format!("[{}]", fields.join(","))
    }
}

impl JsonMap {
    pub fn new(_name: Option<&str>) -> Self {
        JsonMap {
            values: BTreeMap::new(),
        }
    }

    pub fn from_file(fname: &str) -> Result<Self> {
        let text = fs::read_to_string(fname)?;
        let mut in_arg = JsonCharInput::new(&text);
        let token = JsonToken::read_input(&mut in_arg)?;
        if !token.is_delimiter('{') {
            return in_arg.error(&format!(
                "Json file '{}': text should start with '{{'",
                fname
            ));
        }
        let mut json = JsonMap::new(None);
        json.parse(&mut in_arg)?;
        Ok(json)
    }

    fn parse(&mut self, in_arg: &mut JsonCharInput) -> Result<()> {
        assert!(self.values.is_empty());
        let mut first = true;
        loop {
            let mut token = JsonToken::read_input(in_arg)?;
            if token.is_delimiter('}') {
                break;
            }
            if !first {
                if !token.is_delimiter(',') {
                    return in_arg.error("','");
                }
                token = JsonToken::read_input(in_arg)?;
            }
            if token.token_type != JsonTokenType::Name && token.token_type != JsonTokenType::Text {
                return in_arg.error("name or text");
            }
            let name = token.name;
            let colon = JsonToken::read_input(in_arg)?;
            if !colon.is_delimiter(':') {
                return in_arg.error("':'");
            }
            let token = JsonToken::read_input(in_arg)?;
            let value = JsonValue::parse(in_arg, token)?;
            if self.values.insert(name.clone(), value).is_some() {
                bail!("Duplicate name in JSON map: \"{}\"", name);
            }
            first = false;
        }
        Ok(())
    }

    pub fn add_null(&mut self, name: &str) {
        if self
            .values
            .insert(name.to_string(), JsonValue::Null)
            .is_some()
        {
            panic!("Duplicate name in JSON map: \"{}\"", name);
        }
    }

    pub fn add_int(&mut self, name: &str, value: i64) {
        if self
            .values
            .insert(name.to_string(), JsonValue::Int(value))
            .is_some()
        {
            panic!("Duplicate name in JSON map: \"{}\"", name);
        }
    }

    pub fn add_double(&mut self, name: &str, value: f64, decimals: usize) {
        if self
            .values
            .insert(
                name.to_string(),
                JsonValue::Double {
                    n: value,
                    decimals,
                    scientific: false,
                },
            )
            .is_some()
        {
            panic!("Duplicate name in JSON map: \"{}\"", name);
        }
    }

    pub fn add_string(&mut self, name: &str, value: &str) {
        if self
            .values
            .insert(name.to_string(), JsonValue::String(value.to_string()))
            .is_some()
        {
            panic!("Duplicate name in JSON map: \"{}\"", name);
        }
    }

    pub fn add_boolean(&mut self, name: &str, value: bool) {
        if self
            .values
            .insert(name.to_string(), JsonValue::Boolean(value))
            .is_some()
        {
            panic!("Duplicate name in JSON map: \"{}\"", name);
        }
    }

    pub fn add_array(&mut self, name: &str, value: JsonArray) {
        if self
            .values
            .insert(name.to_string(), JsonValue::Array(value))
            .is_some()
        {
            panic!("Duplicate name in JSON map: \"{}\"", name);
        }
    }

    pub fn add_map(&mut self, name: &str, value: JsonMap) {
        if self
            .values
            .insert(name.to_string(), JsonValue::Map(value))
            .is_some()
        {
            panic!("Duplicate name in JSON map: \"{}\"", name);
        }
    }

    pub fn at(&self, name_arg: &str) -> Option<&JsonValue> {
        self.values.get(name_arg)
    }

    pub fn get_keys(&self) -> Vec<String> {
        self.values.keys().cloned().collect()
    }

    pub fn to_json_text(&self) -> String {
        let fields = self
            .values
            .iter()
            .map(|(name, value)| format!("{}:{}", json_to_str(name), value.to_json_text()))
            .collect::<Vec<_>>();
        format!("{{{}}}", fields.join(","))
    }
}

impl JsonValue {
    fn parse(in_arg: &mut JsonCharInput, first_token: JsonToken) -> Result<Self> {
        match first_token.token_type {
            JsonTokenType::Name => {
                let token_name = first_token.name.to_ascii_lowercase();
                if token_name == "null" || token_name == "nan" {
                    Ok(JsonValue::Null)
                } else if token_name == "true" {
                    Ok(JsonValue::Boolean(true))
                } else if token_name == "false" {
                    Ok(JsonValue::Boolean(false))
                } else {
                    Ok(JsonValue::String(first_token.name))
                }
            }
            JsonTokenType::Text => Ok(JsonValue::String(first_token.name)),
            JsonTokenType::Integer => Ok(JsonValue::Int(first_token.n)),
            JsonTokenType::Double => Ok(JsonValue::Double {
                n: first_token.d,
                decimals: first_token.decimals,
                scientific: first_token.scientific,
            }),
            JsonTokenType::Delimiter => {
                if first_token.is_delimiter('{') {
                    let mut json = JsonMap::new(None);
                    json.parse(in_arg)?;
                    Ok(JsonValue::Map(json))
                } else if first_token.is_delimiter('[') {
                    Ok(JsonValue::Array(JsonArray::parse(in_arg)?))
                } else {
                    in_arg.error("Bad delimiter")
                }
            }
            JsonTokenType::Empty => in_arg.error("Bad token"),
        }
    }

    pub fn as_json_null(&self) -> Option<&JsonValue> {
        if matches!(self, JsonValue::Null) {
            Some(self)
        } else {
            None
        }
    }

    pub fn as_json_string(&self) -> Option<&String> {
        if let JsonValue::String(value) = self {
            Some(value)
        } else {
            None
        }
    }

    pub fn as_json_int(&self) -> Option<i64> {
        if let JsonValue::Int(value) = self {
            Some(*value)
        } else {
            None
        }
    }

    pub fn as_json_double(&self) -> Option<f64> {
        if let JsonValue::Double { n, .. } = self {
            Some(*n)
        } else {
            None
        }
    }

    pub fn as_json_boolean(&self) -> Option<bool> {
        if let JsonValue::Boolean(value) = self {
            Some(*value)
        } else {
            None
        }
    }

    pub fn as_json_array(&self) -> Option<&JsonArray> {
        if let JsonValue::Array(value) = self {
            Some(value)
        } else {
            None
        }
    }

    pub fn as_json_map(&self) -> Option<&JsonMap> {
        if let JsonValue::Map(value) = self {
            Some(value)
        } else {
            None
        }
    }

    pub fn get_string(&self) -> String {
        if matches!(self, JsonValue::Null) {
            panic!("undefined");
        }
        if let JsonValue::String(value) = self {
            return value.clone();
        }
        panic!("Not a JsonString");
    }

    pub fn get_int(&self) -> i64 {
        if matches!(self, JsonValue::Null) {
            panic!("undefined");
        }
        if let JsonValue::Int(value) = self {
            return *value;
        }
        panic!("Not a JsonInt");
    }

    pub fn get_double(&self) -> f64 {
        if matches!(self, JsonValue::Null) {
            return f64::NAN;
        }
        if let JsonValue::Double { n, .. } = self {
            return *n;
        }
        panic!("Not a JsonDouble");
    }

    pub fn get_boolean(&self) -> bool {
        if matches!(self, JsonValue::Null) {
            panic!("undefined");
        }
        if let JsonValue::Boolean(value) = self {
            return *value;
        }
        panic!("Not a JsonBoolean");
    }

    pub fn at(&self, name_arg: &str) -> Option<&JsonValue> {
        if matches!(self, JsonValue::Null) {
            return None;
        }
        if let JsonValue::Map(value) = self {
            return value.at(name_arg);
        }
        panic!("Not a JsonMap");
    }

    pub fn at_index(&self, index: usize) -> Option<&JsonValue> {
        if matches!(self, JsonValue::Null) {
            return None;
        }
        if let JsonValue::Array(value) = self {
            return value.at(index);
        }
        panic!("Not a JsonArray");
    }

    pub fn get_size(&self) -> usize {
        if matches!(self, JsonValue::Null) {
            return 0;
        }
        if let JsonValue::Array(value) = self {
            return value.get_size();
        }
        panic!("Not a JsonArray");
    }

    pub fn to_json_text(&self) -> String {
        match self {
            JsonValue::Null => "null".to_string(),
            JsonValue::String(value) => json_to_str(value),
            JsonValue::Int(value) => value.to_string(),
            JsonValue::Double {
                n,
                decimals,
                scientific,
            } => {
                if n.is_finite() {
                    ONumber::new(*decimals, *scientific).to_string_f64(*n)
                } else {
                    "null".to_string()
                }
            }
            JsonValue::Boolean(value) => {
                if *value {
                    "true".to_string()
                } else {
                    "false".to_string()
                }
            }
            JsonValue::Array(value) => value.to_json_text(),
            JsonValue::Map(value) => value.to_json_text(),
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct JsonContainer {
    pub values: BTreeMap<String, JsonValue>,
}

impl JsonContainer {
    pub fn add_null(&mut self, name: &str) {
        if self
            .values
            .insert(name.to_string(), JsonValue::Null)
            .is_some()
        {
            panic!("Duplicate name in JSON map: \"{}\"", name);
        }
    }

    pub fn add_int(&mut self, name: &str, value: i64) {
        if self
            .values
            .insert(name.to_string(), JsonValue::Int(value))
            .is_some()
        {
            panic!("Duplicate name in JSON map: \"{}\"", name);
        }
    }

    pub fn add_double(&mut self, name: &str, value: f64, decimals: usize) {
        if self
            .values
            .insert(
                name.to_string(),
                JsonValue::Double {
                    n: value,
                    decimals,
                    scientific: false,
                },
            )
            .is_some()
        {
            panic!("Duplicate name in JSON map: \"{}\"", name);
        }
    }

    pub fn add_string(&mut self, name: &str, value: &str) {
        if self
            .values
            .insert(name.to_string(), JsonValue::String(value.to_string()))
            .is_some()
        {
            panic!("Duplicate name in JSON map: \"{}\"", name);
        }
    }

    pub fn add_boolean(&mut self, name: &str, value: bool) {
        if self
            .values
            .insert(name.to_string(), JsonValue::Boolean(value))
            .is_some()
        {
            panic!("Duplicate name in JSON map: \"{}\"", name);
        }
    }

    pub fn add_array(&mut self, name: &str, value: JsonArray) {
        if self
            .values
            .insert(name.to_string(), JsonValue::Array(value))
            .is_some()
        {
            panic!("Duplicate name in JSON map: \"{}\"", name);
        }
    }

    pub fn add_map(&mut self, name: &str, value: JsonMap) {
        if self
            .values
            .insert(name.to_string(), JsonValue::Map(value))
            .is_some()
        {
            panic!("Duplicate name in JSON map: \"{}\"", name);
        }
    }

    pub fn at(&self, name_arg: &str) -> Option<&JsonValue> {
        self.values.get(name_arg)
    }

    pub fn get_keys(&self) -> Vec<String> {
        self.values.keys().cloned().collect()
    }

    pub fn to_json_text(&self) -> String {
        let fields = self
            .values
            .iter()
            .map(|(name, value)| format!("{}:{}", json_to_str(name), value.to_json_text()))
            .collect::<Vec<_>>();
        format!("{{{}}}", fields.join(","))
    }
}

impl Date {
    pub fn new(year: i16, month: i8, day: i8) -> Self {
        Date { year, month, day }
    }

    pub fn is_year(n: i16) -> bool {
        n > 1000 && n < 2500
    }

    pub fn is_month(n: i16) -> bool {
        (0..12).contains(&n)
    }

    pub fn is_day(n: i16) -> bool {
        (0..31).contains(&n)
    }

    pub fn parse(s: &str, fmt: DateFormat) -> Date {
        match fmt {
            DateFormat::Year => {
                let mut words = s.split_whitespace();
                if let (Some(year), None) = (words.next(), words.next()) {
                    if let Ok(year) = year.parse::<i16>() {
                        if Date::is_year(year) {
                            return Date::new(year, 0, 0);
                        }
                    }
                } else if let Ok(year) = s.parse::<i16>() {
                    if Date::is_year(year) {
                        return Date::new(year, 0, 0);
                    }
                }
            }
            DateFormat::Ymd => {
                let mut pos = 0usize;
                let mut values = [0i16; 3];
                let mut delimiters = ['\0'; 2];
                let mut ok = true;
                for i in 0..3 {
                    while pos < s.len() {
                        let c = s[pos..].chars().next().unwrap();
                        if !c.is_whitespace() {
                            break;
                        }
                        pos += c.len_utf8();
                    }
                    let number_start = pos;
                    if pos < s.len() {
                        let c = s[pos..].chars().next().unwrap();
                        if c == '+' || c == '-' {
                            pos += c.len_utf8();
                        }
                    }
                    let digits_start = pos;
                    while pos < s.len() {
                        let c = s[pos..].chars().next().unwrap();
                        if !c.is_ascii_digit() {
                            break;
                        }
                        pos += c.len_utf8();
                    }
                    if digits_start == pos {
                        ok = false;
                        break;
                    }
                    match s[number_start..pos].parse::<i16>() {
                        Ok(value) => values[i] = value,
                        Err(_) => {
                            ok = false;
                            break;
                        }
                    }
                    if i < 2 {
                        while pos < s.len() {
                            let c = s[pos..].chars().next().unwrap();
                            if !c.is_whitespace() {
                                break;
                            }
                            pos += c.len_utf8();
                        }
                        if pos >= s.len() {
                            ok = false;
                            break;
                        }
                        let c = s[pos..].chars().next().unwrap();
                        delimiters[i] = c;
                        pos += c.len_utf8();
                    }
                }
                if ok && s[pos..].trim().is_empty() && delimiters[0] == delimiters[1] {
                    let year = values[0];
                    let month = values[1] - 1;
                    let day = values[2] - 1;
                    if Date::is_year(year) && Date::is_month(month) && Date::is_day(day) {
                        return Date::new(year, month as i8, day as i8);
                    }
                }
            }
            DateFormat::None => panic!("Unknown date format"),
        }
        Date::default()
    }

    pub fn empty(&self) -> bool {
        self.year == 0 && self.month == 0 && self.day == 0
    }

    pub fn save_text(&self, os: &mut dyn Write) -> Result<()> {
        write!(
            os,
            "{:04}-{:02}-{:02}",
            self.year,
            i16::from(self.month) + 1,
            i16::from(self.day) + 1
        )?;
        Ok(())
    }

    pub fn to_json(&self, parent: Option<&mut JsonContainer>, name: Option<&str>) -> JsonMap {
        let mut json = JsonMap::new(name);
        json.add_int("year", i64::from(self.year));
        json.add_int("month", i64::from(self.month));
        json.add_int("day", i64::from(self.day));
        if let (Some(parent), Some(name)) = (parent, name) {
            parent.add_map(name, json.clone());
        }
        json
    }

    pub fn less(&self, other: &Date, equal: bool) -> bool {
        if self.year != other.year {
            return self.year < other.year;
        }
        if self.month != other.month {
            return self.month < other.month;
        }
        if self.day != other.day {
            return self.day < other.day;
        }
        equal
    }

    pub fn year_divisible(&self) -> bool {
        self.month == 0 && self.day == 0
    }

    pub fn quarter_divisible(&self) -> bool {
        self.month % 3 == 0 && self.day == 0
    }

    pub fn month_divisible(&self) -> bool {
        self.day == 0
    }
}

impl Default for Date {
    fn default() -> Self {
        Date {
            year: 0,
            month: 0,
            day: 0,
        }
    }
}

impl PartialOrd for Date {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Date {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.year
            .cmp(&other.year)
            .then(self.month.cmp(&other.month))
            .then(self.day.cmp(&other.day))
    }
}

impl std::ops::Sub for Date {
    type Output = Date;

    fn sub(self, other: Date) -> Date {
        let mut d = Date::new(
            self.year - other.year,
            self.month - other.month,
            self.day - other.day,
        );
        while d.month < 0 {
            d.month += 12;
            d.year -= 1;
        }
        d
    }
}

// --- TsvOut ---

/// Tab-separated output writer. Matches the C++ TsvOut class.
pub struct TsvOut<'a> {
    os: Option<&'a mut dyn Write>,
    on: Option<ONumber>,
    lines: usize,
    fields_max: usize,
    fields: usize,
    pub use_pound: bool,
}

impl<'a> TsvOut<'a> {
    pub fn new(os: Option<&'a mut dyn Write>) -> Self {
        Self::with_format(os, 6, false)
    }

    pub fn with_format(os: Option<&'a mut dyn Write>, precision: usize, scientific: bool) -> Self {
        let on = os.as_ref().map(|_| ONumber::new(precision, scientific));
        TsvOut {
            os,
            on,
            lines: 0,
            fields_max: 0,
            fields: 0,
            use_pound: true,
        }
    }

    pub fn live(&self) -> bool {
        self.os.is_some()
    }

    pub fn empty(&self) -> bool {
        self.lines == 0 && self.fields_max == 0 && self.fields == 0
    }

    pub fn write_field<T: TsvOutField + ?Sized>(&mut self, field: &T) -> Result<()> {
        if let Some(ref mut os) = self.os {
            if self.lines > 0 && self.fields >= self.fields_max {
                bail!("TsvOut: fields_max = {}", self.fields_max);
            }
            if self.fields > 0 {
                write!(os, "\t")?;
            } else if self.lines == 0 && self.use_pound {
                write!(os, "#")?;
            }
            if let Some(on) = &self.on {
                field.write_tsv_out(os, on)?;
            }
            self.fields += 1;
        }
        Ok(())
    }

    pub fn new_ln(&mut self) -> Result<()> {
        if let Some(ref mut os) = self.os {
            writeln!(os)?;
            if self.lines == 0 {
                self.fields_max = self.fields;
            }
            self.lines += 1;
            if self.fields != self.fields_max {
                bail!(
                    "TsvOut: fields_max = {}, but fields = {}",
                    self.fields_max,
                    self.fields
                );
            }
            self.fields = 0;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ONumber {
    precision: usize,
    scientific: bool,
}

impl ONumber {
    pub fn new(precision: usize, scientific: bool) -> Self {
        ONumber {
            precision,
            scientific,
        }
    }

    pub fn write_f64(&self, os: &mut dyn Write, value: f64) -> Result<()> {
        if self.scientific {
            write!(os, "{}", self.format_scientific(value))?;
        } else {
            write!(os, "{value:.*}", self.precision)?;
        }
        Ok(())
    }

    pub fn to_string_f64(&self, value: f64) -> String {
        if self.scientific {
            self.format_scientific(value)
        } else {
            format!("{value:.*}", self.precision)
        }
    }

    fn format_scientific(&self, value: f64) -> String {
        let s = format!("{value:.*e}", self.precision);
        let Some((mantissa, exponent)) = s.split_once('e') else {
            return s;
        };
        let Ok(exponent) = exponent.parse::<i32>() else {
            return s;
        };
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{mantissa}e{sign}{:02}", exponent.abs())
    }
}

pub trait TsvOutField {
    fn write_tsv_out(&self, os: &mut dyn Write, on: &ONumber) -> Result<()>;
}

impl TsvOutField for str {
    fn write_tsv_out(&self, os: &mut dyn Write, _on: &ONumber) -> Result<()> {
        write!(os, "{self}")?;
        Ok(())
    }
}

impl TsvOutField for &str {
    fn write_tsv_out(&self, os: &mut dyn Write, on: &ONumber) -> Result<()> {
        (*self).write_tsv_out(os, on)
    }
}

impl TsvOutField for String {
    fn write_tsv_out(&self, os: &mut dyn Write, on: &ONumber) -> Result<()> {
        self.as_str().write_tsv_out(os, on)
    }
}

impl TsvOutField for Date {
    fn write_tsv_out(&self, os: &mut dyn Write, _on: &ONumber) -> Result<()> {
        self.save_text(os)
    }
}

impl TsvOutField for Header {
    fn write_tsv_out(&self, os: &mut dyn Write, _on: &ONumber) -> Result<()> {
        self.save_text(os)
    }
}

impl TsvOutField for TextTable {
    fn write_tsv_out(&self, os: &mut dyn Write, _on: &ONumber) -> Result<()> {
        self.save_text(os)
    }
}

macro_rules! impl_tsv_out_display {
    ($($t:ty),* $(,)?) => {
        $(
            impl TsvOutField for $t {
                fn write_tsv_out(
                    &self,
                    os: &mut dyn Write,
                    _on: &ONumber,
                ) -> Result<()> {
                    write!(os, "{self}")?;
                    Ok(())
                }
            }
        )*
    };
}

macro_rules! impl_tsv_out_float {
    ($($t:ty),* $(,)?) => {
        $(
            impl TsvOutField for $t {
                fn write_tsv_out(
                    &self,
                    os: &mut dyn Write,
                    on: &ONumber,
                ) -> Result<()> {
                    on.write_f64(os, *self as f64)
                }
            }
        )*
    };
}

impl_tsv_out_display!(bool, char, i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize);
impl_tsv_out_float!(f32, f64);

impl Drop for TsvOut<'_> {
    fn drop(&mut self) {
        if self.fields > 0 && !std::thread::panicking() {
            panic!("TsvOut: unfinished line with {} fields", self.fields);
        }
    }
}

// --- TextTable ---

/// Column header metadata
#[derive(Debug, Clone)]
pub struct Header {
    pub name: String,
    pub len_max: usize,
    pub numeric: bool,
    pub scientific: bool,
    pub decimals: usize,
    pub null: bool,
    pub choices: BTreeSet<String>,
}

impl Header {
    pub fn new(name: &str) -> Self {
        Header {
            name: name.to_string(),
            len_max: 0,
            numeric: true,
            scientific: false,
            decimals: 0,
            null: false,
            choices: BTreeSet::new(),
        }
    }

    pub fn qc(&self) -> Result<()> {
        if self.name.is_empty() {
            bail!("Empty name");
        }
        if self.scientific && !self.numeric {
            bail!("scientific implies numeric");
        }
        if self.decimals > 0 && !self.numeric {
            bail!("decimals implies numeric");
        }
        Ok(())
    }

    pub fn save_text(&self, os: &mut dyn Write) -> Result<()> {
        write!(
            os,
            "{}\t{}\t{}\t{}",
            self.name,
            self.len_max,
            if self.numeric {
                format!(
                    "{}({})",
                    if self.scientific { "float" } else { "int" },
                    self.decimals
                )
            } else {
                "char".to_string()
            },
            if self.null { "null" } else { "not null" }
        )?;
        Ok(())
    }

    pub fn save_sql(&self, os: &mut dyn Write) -> Result<()> {
        let name = self.name.replace(' ', "_").replace('%', "Percent");
        write!(os, "{} ", name)?;
        if self.numeric {
            if self.scientific {
                write!(os, "float")?;
            } else {
                write!(os, "numeric({},{})", self.len_max.max(1), self.decimals)?;
            }
        } else {
            write!(os, "varchar({})", self.len_max.max(1))?;
        }
        if !self.null {
            write!(os, " not null")?;
        }
        writeln!(os)?;
        Ok(())
    }
}

/// Tab-separated value table with header.
/// Matches the C++ TextTable class.
pub struct TextTable {
    pub name: String,
    pub pound: bool,
    pub save_header: bool,
    pub header: Vec<Header>,
    pub rows: Vec<Vec<String>>,
}

impl TextTable {
    pub fn str2header(s: &str, sep: char) -> Vec<Header> {
        if s.is_empty() {
            return Vec::new();
        }
        s.split(sep).map(|name| Header::new(name.trim())).collect()
    }

    /// Parse a TSV file
    pub fn from_file(fname: &str) -> Result<Self> {
        Self::from_file_with_options(fname, "", true, 0)
    }

    pub fn from_file_with_synonyms(fname: &str, column_synonyms_fname: &str) -> Result<Self> {
        Self::from_file_with_options(fname, column_synonyms_fname, true, 0)
    }

    pub fn from_file_with_header(fname: &str, header_p: bool) -> Result<Self> {
        Self::from_file_with_options(fname, "", header_p, 0)
    }

    pub fn from_file_with_options(
        fname: &str,
        column_synonyms_fname: &str,
        header_p: bool,
        _display_period: u32,
    ) -> Result<Self> {
        let file = File::open(fname)?;
        let reader = BufReader::new(file);
        let mut header: Vec<Header> = Vec::new();
        let mut rows: Vec<Vec<String>> = Vec::new();
        let mut pound = false;
        let mut data_line: Option<String> = None;

        let mut lines = reader.lines();
        while let Some(line) = lines.next() {
            let mut line = line?;
            trim_trailing_in_place(&mut line);
            if line.is_empty() {
                continue;
            }

            let this_pound = line.starts_with('#');
            if this_pound {
                pound = true;
                line.remove(0);
            }
            if line.is_empty() {
                continue;
            }

            if header.is_empty() || this_pound {
                header = TextTable::str2header(&line, '\t');
            }
            assert!(!header.is_empty());

            if header_p {
                if !this_pound {
                    if !pound {
                        data_line = match lines.next() {
                            Some(line) => Some(line?),
                            None => None,
                        };
                    } else {
                        data_line = Some(line);
                    }
                    break;
                }
            } else {
                for (i, h) in header.iter_mut().enumerate() {
                    h.name = (i + 1).to_string();
                }
                data_line = Some(line);
                pound = true;
                break;
            }
        }
        if header.is_empty() {
            bail!("Cannot read the table header\nIn table file: {}", fname);
        }

        while let Some(mut line) = data_line {
            trim_trailing_in_place(&mut line);
            if !line.is_empty() {
                let mut row = split_tsv_trimmed(&line);
                while row.len() < header.len() {
                    row.push(String::new());
                }
                rows.push(row);
            }
            data_line = match lines.next() {
                Some(line) => Some(line?),
                None => None,
            };
        }

        let mut table = TextTable {
            name: fname.to_string(),
            pound,
            save_header: true,
            header,
            rows,
        };
        if !column_synonyms_fname.is_empty() {
            table.apply_column_synonyms(column_synonyms_fname)?;
        }
        table.set_header()?;
        Ok(table)
    }

    /// Create an empty table with the given header
    pub fn with_header(header: Vec<Header>) -> Self {
        TextTable {
            name: String::new(),
            pound: false,
            save_header: true,
            header,
            rows: Vec::new(),
        }
    }

    pub fn with_pound_header(pound: bool, header: Vec<Header>) -> Self {
        TextTable {
            name: String::new(),
            pound,
            save_header: true,
            header,
            rows: Vec::new(),
        }
    }

    fn apply_column_synonyms(&mut self, column_synonyms_fname: &str) -> Result<()> {
        let file = File::open(column_synonyms_fname)?;
        let reader = BufReader::new(file);
        let mut main_syn = String::new();
        for line in reader.lines() {
            let line = line?;
            let syn = line.trim();
            if syn.is_empty() {
                main_syn.clear();
            } else if main_syn.is_empty() {
                main_syn = syn.to_string();
            } else if main_syn != syn {
                if let Some(i) = self.col2num_(syn) {
                    if self.has_column(&main_syn) {
                        bail!(
                            "Table \"{}\": Column \"{}\" already exists",
                            self.name,
                            main_syn
                        );
                    } else {
                        self.header[i].name = main_syn.clone();
                    }
                }
            }
        }
        Ok(())
    }

    /// Find column index by name, returns no column if not found.
    pub fn col2num_(&self, column_name: &str) -> Option<usize> {
        self.header.iter().position(|h| h.name == column_name)
    }

    /// Find column index by name, returns error if not found.
    pub fn col2num(&self, column_name: &str) -> Result<usize> {
        self.col2num_(column_name).ok_or_else(|| {
            anyhow::anyhow!(
                "Table has no column \"{}\"\nIn table file: {}",
                column_name,
                self.name
            )
        })
    }

    /// Check if table has a column
    pub fn has_column(&self, column_name: &str) -> bool {
        self.col2num_(column_name).is_some()
    }

    pub fn columns2nums(&self, columns: &[&str]) -> Result<Vec<usize>> {
        columns.iter().map(|name| self.col2num(name)).collect()
    }

    pub fn qc(&self) -> Result<()> {
        let mut names: Vec<&str> = Vec::with_capacity(self.header.len());
        for (i, h) in self.header.iter().enumerate() {
            h.qc()
                .map_err(|e| anyhow::anyhow!("Header column #{}: {}", i + 1, e))?;
            names.push(h.name.as_str());
        }
        names.sort_unstable();
        for pair in names.windows(2) {
            if pair[0] == pair[1] {
                bail!(
                    "Duplicate column name: \"{}\"\nIn table file: {}",
                    pair[0],
                    self.name
                );
            }
        }
        for (row_num, row) in self.rows.iter().enumerate() {
            if row.len() != self.header.len() {
                bail!(
                    "Row {} contains {} columns whereas table has {} columns\nIn table file: {}",
                    row_num + 1,
                    row.len(),
                    self.header.len(),
                    self.name
                );
            }
            for field in row {
                if field.contains('\t') {
                    bail!(
                        "Field \"{}\" of row {} contains a tab character\nIn table file: {}",
                        self.header[row_num].name,
                        row_num + 1,
                        self.name
                    );
                }
                if field.contains('\n') {
                    bail!(
                        "Field \"{}\" of row {} contains an EOL character\nIn table file: {}",
                        self.header[row_num].name,
                        row_num + 1,
                        self.name
                    );
                }
            }
        }
        Ok(())
    }

    pub fn set_header(&mut self) -> Result<()> {
        for (row_num, row) in self.rows.iter().enumerate() {
            if row.len() != self.header.len() {
                bail!(
                    "Row {} contains {} columns whereas header has {} columns",
                    row_num + 1,
                    row.len(),
                    self.header.len()
                );
            }
            for (i, field) in row.iter().enumerate() {
                let field = field.trim();
                let h = &mut self.header[i];
                if crate::common::str_null(field) {
                    h.null = true;
                    continue;
                }
                h.len_max = h.len_max.max(field.len());
                if h.choices.len() <= CHOICES_MAX {
                    h.choices.insert(field.to_string());
                }
                if h.numeric {
                    let mut endptr: *mut c_char = std::ptr::null_mut();
                    let numeric = CString::new(field).is_ok_and(|field_c| unsafe {
                        strtod(field_c.as_ptr(), &mut endptr);
                        endptr == field_c.as_ptr().add(field.len()) as *mut c_char
                    });
                    if !numeric {
                        h.numeric = false;
                        h.scientific = false;
                        h.decimals = 0;
                    }
                }
                if h.numeric {
                    let (scientific, _has_point, decimals) = crate::common::get_scientific(field);
                    if scientific {
                        h.scientific = true;
                    }
                    h.decimals = h.decimals.max(decimals);
                }
            }
        }
        for row in &self.rows {
            for (i, field) in row.iter().enumerate() {
                let field = field.trim();
                if crate::common::str_null(field) {
                    continue;
                }
                let h = &mut self.header[i];
                if !h.numeric {
                    continue;
                }
                let (_scientific, has_point, decimals) = crate::common::get_scientific(field);
                h.len_max = h
                    .len_max
                    .max(field.len() + (h.decimals - decimals) + usize::from(!has_point));
            }
        }
        Ok(())
    }

    pub fn compare(&self, row1: &[String], row2: &[String], column: usize) -> i32 {
        let s1 = row1.get(column).map(String::as_str).unwrap_or("");
        let s2 = row2.get(column).map(String::as_str).unwrap_or("");
        if self.header[column].numeric {
            let a = if crate::common::str_null(s1) {
                0.0
            } else {
                let s1_c = CString::new(s1).expect("numeric TextTable column");
                let mut endptr: *mut c_char = std::ptr::null_mut();
                let value = unsafe { strtod(s1_c.as_ptr(), &mut endptr) };
                if endptr == s1_c.as_ptr() as *mut c_char {
                    panic!("numeric TextTable column");
                }
                value
            };
            let b = if crate::common::str_null(s2) {
                0.0
            } else {
                let s2_c = CString::new(s2).expect("numeric TextTable column");
                let mut endptr: *mut c_char = std::ptr::null_mut();
                let value = unsafe { strtod(s2_c.as_ptr(), &mut endptr) };
                if endptr == s2_c.as_ptr() as *mut c_char {
                    panic!("numeric TextTable column");
                }
                value
            };
            if a < b {
                return -1;
            }
            if a > b {
                return 1;
            }
            return 0;
        }
        if s1 < s2 {
            -1
        } else if s1 > s2 {
            1
        } else {
            0
        }
    }

    /// Sort rows by given columns
    pub fn sort(&mut self, by: &[&str]) -> Result<()> {
        let col_nums: Vec<usize> = by
            .iter()
            .map(|name| self.col2num(name))
            .collect::<Result<Vec<_>>>()?;

        let mut rows = std::mem::take(&mut self.rows);
        rows.sort_by(|a, b| {
            for &col in &col_nums {
                let cmp = self.compare(a, b, col);
                if cmp < 0 {
                    return std::cmp::Ordering::Less;
                }
                if cmp > 0 {
                    return std::cmp::Ordering::Greater;
                }
            }
            for col in 0..self.header.len() {
                let cmp = self.compare(a, b, col);
                if cmp < 0 {
                    return std::cmp::Ordering::Less;
                }
                if cmp > 0 {
                    return std::cmp::Ordering::Greater;
                }
            }
            std::cmp::Ordering::Equal
        });
        self.rows = rows;

        Ok(())
    }

    pub fn deredundify<F>(&mut self, equiv_cols: &[&str], equiv_better: F) -> Result<()>
    where
        F: Fn(&[String], &[String]) -> i32,
    {
        let col_nums: Vec<usize> = equiv_cols
            .iter()
            .map(|name| self.col2num(name))
            .collect::<Result<Vec<_>>>()?;

        let mut rows = std::mem::take(&mut self.rows);
        let lt = |row1: &[String], row2: &[String]| -> bool {
            for &col in &col_nums {
                match self.compare(row1, row2, col) {
                    -1 => return true,
                    1 => return false,
                    _ => {}
                }
            }
            false
        };

        rows.sort_by(|a, b| {
            if lt(a, b) {
                std::cmp::Ordering::Less
            } else if lt(b, a) {
                std::cmp::Ordering::Greater
            } else {
                std::cmp::Ordering::Equal
            }
        });

        let mut to_delete = vec![false; rows.len()];
        let mut i = 0;
        while i < rows.len() {
            let mut j = i + 1;
            while j <= rows.len() {
                if j == rows.len() {
                    i = rows.len();
                    break;
                }
                if lt(&rows[i], &rows[j]) {
                    i = j;
                    break;
                }
                let mut k = i;
                while k < j {
                    if !to_delete[k] {
                        let mut stop = false;
                        let mut l = k + 1;
                        while l <= j {
                            if !to_delete[l] {
                                if equiv_better(&rows[k], &rows[l]) == 1 {
                                    to_delete[l] = true;
                                } else if equiv_better(&rows[l], &rows[k]) == 1 {
                                    to_delete[k] = true;
                                    stop = true;
                                    break;
                                }
                            }
                            l += 1;
                        }
                        if stop {
                            break;
                        }
                    }
                    k += 1;
                }
                j += 1;
            }
        }

        if to_delete.iter().any(|delete| *delete) {
            let mut keep = Vec::with_capacity(rows.len());
            for (row, delete) in rows.into_iter().zip(to_delete) {
                if !delete {
                    keep.push(row);
                }
            }
            rows = keep;
        }
        self.rows = rows;

        Ok(())
    }

    /// Filter to only the specified columns (in order)
    pub fn filter_columns(&mut self, new_column_names: &[&str]) -> Result<()> {
        let col_nums: Vec<usize> = self.columns2nums(new_column_names)?;

        let new_header: Vec<Header> = col_nums.iter().map(|&i| self.header[i].clone()).collect();

        let new_rows: Vec<Vec<String>> = self
            .rows
            .iter()
            .map(|row| {
                col_nums
                    .iter()
                    .map(|&i| row.get(i).cloned().unwrap_or_default())
                    .collect()
            })
            .collect();

        self.header = new_header;
        self.rows = new_rows;

        Ok(())
    }

    pub fn duplicate_column(&mut self, column_name_from: &str, column_name_to: &str) -> Result<()> {
        assert!(!column_name_to.is_empty());
        let from = self.col2num(column_name_from)?;
        if self.has_column(column_name_to) {
            bail!("Table already has column \"{}\"", column_name_to);
        }
        self.header.push(self.header[from].clone());
        self.header.last_mut().unwrap().name = column_name_to.to_string();
        for row in &mut self.rows {
            row.push(row[from].clone());
        }
        self.qc()?;
        Ok(())
    }

    pub fn substitue_column(
        &mut self,
        column_name_from: &mut String,
        column_name_to: &str,
    ) -> Result<()> {
        self.duplicate_column(column_name_from, column_name_to)?;
        *column_name_from = column_name_to.to_string();
        Ok(())
    }

    pub fn is_key(&self, col_num: usize) -> bool {
        assert!(col_num < self.header.len());

        let h = &self.header[col_num];
        if h.null {
            return false;
        }
        if h.numeric && (h.scientific || h.decimals > 0) {
            return false;
        }

        let mut values = HashSet::with_capacity(self.rows.len());
        for row in &self.rows {
            assert!(!row[col_num].is_empty());
            if !values.insert(row[col_num].as_str()) {
                return false;
            }
        }

        true
    }

    pub fn aggr2values(aggr: &str) -> Vec<String> {
        if aggr.is_empty() {
            return Vec::new();
        }
        let mut values: Vec<String> = aggr.split(AGGR_SEP).map(|s| s.trim().to_string()).collect();
        values.sort();
        values.dedup();
        values
    }

    pub fn group(
        &mut self,
        by: &[&str],
        sum: &[&str],
        min_v: &[&str],
        max_v: &[&str],
        aggr: &[&str],
    ) -> Result<()> {
        let by_index = self.columns2nums(by)?;
        let sum_index = self.columns2nums(sum)?;
        let min_index = self.columns2nums(min_v)?;
        let max_index = self.columns2nums(max_v)?;
        let aggr_index = self.columns2nums(aggr)?;

        let mut name2oper: HashMap<&str, &str> = HashMap::new();
        for (columns, oper) in [
            (by, "group by"),
            (sum, "sum"),
            (min_v, "min"),
            (max_v, "max"),
            (aggr, "aggregation"),
        ] {
            for &column in columns {
                if let Some(prev) = name2oper.insert(column, oper) {
                    bail!(
                        "Column \"{}\" is used for operations \"{}\" and \"{}\"",
                        column,
                        prev,
                        oper
                    );
                }
            }
        }
        for &column in sum {
            if !self.header[self.col2num(column)?].numeric {
                bail!("Summation column \"{}\" is not numeric", column);
            }
        }

        self.sort(by)?;

        let mut i = 0usize;
        for j in 1..self.rows.len() {
            assert!(i < j);
            if rows_same(&self.rows[i], &self.rows[j], &by_index) {
                self.merge(i, j, &sum_index, &min_index, &max_index, &aggr_index)?;
            } else {
                i += 1;
                if i < j {
                    self.rows[i] = std::mem::take(&mut self.rows[j]);
                }
            }
        }
        if !self.rows.is_empty() {
            i += 1;
        }
        self.rows.truncate(i);

        let mut new_columns =
            Vec::with_capacity(by.len() + sum.len() + min_v.len() + max_v.len() + aggr.len());
        new_columns.extend_from_slice(by);
        new_columns.extend_from_slice(sum);
        new_columns.extend_from_slice(min_v);
        new_columns.extend_from_slice(max_v);
        new_columns.extend_from_slice(aggr);
        self.filter_columns(&new_columns)?;
        Ok(())
    }

    fn merge(
        &mut self,
        to_row_num: usize,
        from_row_num: usize,
        sum: &[usize],
        min_v: &[usize],
        max_v: &[usize],
        aggr: &[usize],
    ) -> Result<()> {
        assert!(to_row_num < from_row_num);

        for &i in sum {
            let h = &self.header[i];
            assert!(h.numeric);
            let s1 = &self.rows[to_row_num][i];
            let s2 = &self.rows[from_row_num][i];
            let d1 = if crate::common::str_null(s1) {
                0.0
            } else {
                let s1_c = CString::new(s1.as_str())?;
                let mut endptr: *mut c_char = std::ptr::null_mut();
                let value = unsafe { strtod(s1_c.as_ptr(), &mut endptr) };
                if endptr == s1_c.as_ptr() as *mut c_char {
                    bail!("Cannot convert \"{}\" to number", s1);
                }
                value
            };
            let d2 = if crate::common::str_null(s2) {
                0.0
            } else {
                let s2_c = CString::new(s2.as_str())?;
                let mut endptr: *mut c_char = std::ptr::null_mut();
                let value = unsafe { strtod(s2_c.as_ptr(), &mut endptr) };
                if endptr == s2_c.as_ptr() as *mut c_char {
                    bail!("Cannot convert \"{}\" to number", s2);
                }
                value
            };
            self.rows[to_row_num][i] =
                ONumber::new(h.decimals, h.scientific).to_string_f64(d1 + d2);
        }

        for &i in min_v {
            if self.rows[to_row_num][i].is_empty()
                || (!self.rows[from_row_num][i].is_empty()
                    && self.compare(&self.rows[to_row_num], &self.rows[from_row_num], i) == 1)
            {
                self.rows[to_row_num][i] = self.rows[from_row_num][i].clone();
            }
        }

        for &i in max_v {
            if self.rows[to_row_num][i].is_empty()
                || (!self.rows[from_row_num][i].is_empty()
                    && self.compare(&self.rows[to_row_num], &self.rows[from_row_num], i) == -1)
            {
                self.rows[to_row_num][i] = self.rows[from_row_num][i].clone();
            }
        }

        for &i in aggr {
            if self.rows[from_row_num][i].is_empty() {
                continue;
            }
            if self.rows[from_row_num][i].contains(AGGR_SEP) {
                bail!(
                    "Cannot aggregate column {} for row {} because it contains \"{}\"",
                    self.header[i].name,
                    from_row_num + 1,
                    AGGR_SEP
                );
            }
            if self.rows[to_row_num][i].is_empty() {
                self.rows[to_row_num][i] = self.rows[from_row_num][i].clone();
            } else {
                let addition = self.rows[from_row_num][i].clone();
                aggregate(&mut self.rows[to_row_num][i], &addition, AGGR_SEP);
            }
        }

        Ok(())
    }

    pub fn col_nums_row2values(&self, col_nums: &[usize], row_num: usize) -> Vec<String> {
        let row = &self.rows[row_num];
        col_nums.iter().map(|&i| row[i].clone()).collect()
    }

    pub fn find(
        &self,
        col_nums: &[usize],
        target_values: &[String],
        row_num_start: usize,
    ) -> Option<usize> {
        assert_eq!(col_nums.len(), target_values.len());
        for i in row_num_start..self.rows.len() {
            if self.col_nums_row2values(col_nums, i) == target_values {
                return Some(i);
            }
        }
        None
    }

    pub fn col2values(&self, col: usize) -> Vec<String> {
        assert!(col < self.header.len());
        let mut values = BTreeSet::new();
        for row in &self.rows {
            if !row[col].is_empty() {
                values.insert(row[col].clone());
            }
        }
        values.into_iter().collect()
    }

    pub fn find_date(&self, fmt: &mut DateFormat) -> Option<usize> {
        for (date_col, h) in self.header.iter().enumerate() {
            if h.null || h.scientific {
                continue;
            }
            for fmt_ in [DateFormat::Year, DateFormat::Ymd] {
                *fmt = fmt_;
                let mut is_date = true;
                for row in &self.rows {
                    if Date::parse(&row[date_col], *fmt).empty() {
                        is_date = false;
                        break;
                    }
                }
                if is_date {
                    return Some(date_col);
                }
            }
        }
        None
    }

    pub fn save_text(&self, out: &mut dyn Write) -> Result<()> {
        if self.save_header {
            if self.pound {
                write!(out, "#")?;
            }
            let names: Vec<&str> = self.header.iter().map(|h| h.name.as_str()).collect();
            writeln!(out, "{}", names.join("\t"))?;
        }
        for row in &self.rows {
            writeln!(out, "{}", row.join("\t"))?;
        }
        Ok(())
    }

    pub fn print_header(&self, out: &mut dyn Write) -> Result<()> {
        for (i, h) in self.header.iter().enumerate() {
            write!(out, "{}\t", i + 1)?;
            h.save_text(out)?;
            writeln!(out)?;
        }
        Ok(())
    }

    /// Build a key index on given columns
    pub fn build_key(&self, columns: &[&str]) -> Result<TextTableKey> {
        let col_nums: Vec<usize> = self.columns2nums(columns)?;

        let mut data: HashMap<Vec<String>, usize> = HashMap::new();
        for (row_num, row) in self.rows.iter().enumerate() {
            let key: Vec<String> = col_nums
                .iter()
                .map(|&i| row.get(i).cloned().unwrap_or_default())
                .collect();
            if key.iter().any(|value| value.is_empty()) {
                bail!(
                    "Empty value in key, in row {}\nIn table file: {}",
                    row_num + 1,
                    self.name
                );
            }
            if data.contains_key(&key) {
                bail!(
                    "Duplicate key {} for the key on {}\nIn table file: {}",
                    key.join(","),
                    columns.join(","),
                    self.name
                );
            }
            data.insert(key, row_num);
        }

        Ok(TextTableKey { col_nums, data })
    }

    pub fn build_index(&self, columns: &[&str]) -> Result<TextTableIndex> {
        let col_nums: Vec<usize> = self.columns2nums(columns)?;

        let mut data: HashMap<Vec<String>, Vec<usize>> = HashMap::new();
        for row_num in 0..self.rows.len() {
            let key = self.col_nums_row2values(&col_nums, row_num);
            data.entry(key).or_default().push(row_num);
        }

        Ok(TextTableIndex { col_nums, data })
    }
}

/// Key index for fast row lookup
pub struct TextTableKey {
    pub col_nums: Vec<usize>,
    pub data: HashMap<Vec<String>, usize>,
}

impl TextTableKey {
    pub fn find(&self, values: &[String]) -> Option<usize> {
        self.data.get(values).copied()
    }
}

/// Non-unique index for fast row lookup
pub struct TextTableIndex {
    pub col_nums: Vec<usize>,
    pub data: HashMap<Vec<String>, Vec<usize>>,
}

impl TextTableIndex {
    pub fn find(&self, values: &[String]) -> Option<&Vec<usize>> {
        self.data.get(values)
    }
}

fn rows_same(row1: &[String], row2: &[String], col_nums: &[usize]) -> bool {
    col_nums.iter().all(|&i| row1[i] == row2[i])
}

fn trim_trailing_in_place(s: &mut String) {
    let len = s.trim_end().len();
    s.truncate(len);
}

fn split_tsv_trimmed(s: &str) -> Vec<String> {
    s.split('\t')
        .map(|field| field.trim().to_string())
        .collect()
}

fn aggregate(aggregation: &mut String, addition: &str, aggr_sep: char) {
    let mut values: Vec<String> = aggregation
        .split(aggr_sep)
        .map(|s| s.trim().to_string())
        .collect();
    values.push(addition.to_string());
    values.sort();
    values.dedup();
    *aggregation = values.join(&aggr_sep.to_string());
}

fn json_to_str(s: &str) -> String {
    if !s.is_empty()
        && s.chars().all(|c| c.is_ascii_digit())
        && (s.len() == 1 || !s.starts_with('0'))
    {
        return s.to_string();
    }
    let mut escaped = String::with_capacity(s.len() + 2);
    escaped.push('"');
    for c in s.chars() {
        match c {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\t"),
            '\t' => escaped.push_str("\\t"),
            c => escaped.push(c),
        }
    }
    escaped.push('"');
    escaped
}

#[cfg(test)]
mod tests {
    use super::{Date, JsonContainer, JsonValue};

    #[test]
    fn json_container_get_keys_matches_original_container_surface() {
        let mut parent = JsonContainer::default();
        Date::new(2024, 0, 0).to_json(Some(&mut parent), Some("year"));
        Date::new(2024, 6, 27).to_json(Some(&mut parent), Some("date"));

        assert_eq!(
            parent.get_keys(),
            vec!["date".to_string(), "year".to_string()]
        );
        assert!(matches!(parent.at("date"), Some(JsonValue::Map(_))));
        assert!(matches!(parent.at("year"), Some(JsonValue::Map(_))));
    }
}
