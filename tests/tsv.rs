use std::fs;

use amrfinder::tsv::{
    Date, DateFormat, Header, JsonArray, JsonContainer, JsonMap, JsonValue, TextTable, TsvOut,
};

#[test]
fn tsv_out() {
    let mut buf = Vec::new();
    {
        let mut tsv = TsvOut::new(Some(&mut buf));
        tsv.write_field(&"col1").unwrap();
        tsv.write_field(&"col2").unwrap();
        tsv.write_field(&"col3").unwrap();
        tsv.new_ln().unwrap();
        tsv.write_field(&"a").unwrap();
        tsv.write_field(&"b").unwrap();
        tsv.write_field(&"c").unwrap();
        tsv.new_ln().unwrap();
    }
    let output = String::from_utf8(buf).unwrap();
    assert_eq!(output, "#col1\tcol2\tcol3\na\tb\tc\n");
}

#[test]
fn tsv_out_constructor_format_controls_numbers() {
    let mut buf = Vec::new();
    {
        let mut tsv = TsvOut::with_format(Some(&mut buf), 2, true);
        tsv.write_field(&"value").unwrap();
        tsv.write_field(&1.25_f64).unwrap();
        tsv.new_ln().unwrap();
    }
    assert_eq!(String::from_utf8(buf).unwrap(), "#value\t1.25e+00\n");
}

#[test]
fn tsv_out_new_ln_and_drop_match_original_line_completion_rules() {
    let mut buf = Vec::new();
    {
        let mut tsv = TsvOut::new(Some(&mut buf));
        tsv.write_field(&"value").unwrap();
        tsv.new_ln().unwrap();
    }
    assert_eq!(String::from_utf8(buf).unwrap(), "#value\n");

    let unfinished = std::panic::catch_unwind(|| {
        let mut buf = Vec::new();
        let mut tsv = TsvOut::new(Some(&mut buf));
        tsv.write_field(&"unfinished").unwrap();
    });
    assert!(unfinished.is_err());
}

#[test]
fn tsv_out_field_writes_date_with_translated_save_text() {
    let mut buf = Vec::new();
    {
        let mut tsv = TsvOut::new(Some(&mut buf));
        tsv.write_field(&"date").unwrap();
        tsv.new_ln().unwrap();
        tsv.write_field(&Date::new(2024, 6, 27)).unwrap();
        tsv.new_ln().unwrap();
    }
    assert_eq!(String::from_utf8(buf).unwrap(), "#date\n2024-07-28\n");
}

#[test]
fn tsv_out_field_writes_header_with_translated_save_text() {
    let mut header = Header::new("score");
    header.len_max = 8;
    header.numeric = true;
    header.scientific = true;
    header.decimals = 3;
    header.null = true;

    let mut buf = Vec::new();
    {
        let mut tsv = TsvOut::new(Some(&mut buf));
        tsv.write_field(&header).unwrap();
        tsv.new_ln().unwrap();
    }

    assert_eq!(
        String::from_utf8(buf).unwrap(),
        "#score\t8\tfloat(3)\tnull\n"
    );
}

#[test]
fn tsv_out_field_writes_text_table_with_translated_save_text() {
    let mut table = TextTable::with_header(vec![Header::new("name"), Header::new("value")]);
    table.rows = vec![
        vec!["alpha".to_string(), "1".to_string()],
        vec!["beta".to_string(), "2".to_string()],
    ];

    let mut buf = Vec::new();
    {
        let mut tsv = TsvOut::new(Some(&mut buf));
        tsv.write_field(&table).unwrap();
        tsv.new_ln().unwrap();
    }

    assert_eq!(
        String::from_utf8(buf).unwrap(),
        "#name\tvalue\nalpha\t1\nbeta\t2\n\n"
    );
}

#[test]
fn text_table_basic() {
    let table = TextTable::with_header(vec![Header::new("Name"), Header::new("Value")]);
    assert_eq!(table.header.len(), 2);
    assert_eq!(table.col2num("Name").unwrap(), 0);
    assert_eq!(table.col2num("Value").unwrap(), 1);
    assert!(table.col2num("Missing").is_err());
    assert_eq!(table.col2num_("Missing"), None);
}

#[test]
fn text_table_str2header_empty_input_matches_original_string_vector() {
    assert!(TextTable::str2header("", ',').is_empty());
    assert_eq!(
        TextTable::str2header(" gene , value ", ',')
            .iter()
            .map(|h| h.name.as_str())
            .collect::<Vec<_>>(),
        vec!["gene", "value"]
    );
}

#[test]
fn date_parse_matches_stream_style_ymd() {
    assert_eq!(
        Date::parse("2024-1-2", DateFormat::Ymd),
        Date::new(2024, 0, 1)
    );
    assert_eq!(
        Date::parse(" 2024 - 1 - 2 ", DateFormat::Ymd),
        Date::new(2024, 0, 1)
    );
    assert_eq!(
        Date::parse("+2024-+1-+2", DateFormat::Ymd),
        Date::new(2024, 0, 1)
    );
    assert_eq!(
        Date::parse("2024/12/31", DateFormat::Ymd),
        Date::new(2024, 11, 30)
    );
    assert_eq!(
        Date::parse(" 2024 ", DateFormat::Year),
        Date::new(2024, 0, 0)
    );
    assert!(Date::parse("999-1-1", DateFormat::Ymd).empty());
}

#[test]
fn date_to_json_matches_cpp_hook_shape() {
    let mut parent = JsonContainer::default();
    let json = Date::new(2024, 6, 27).to_json(Some(&mut parent), Some("date"));
    assert_eq!(
        json.to_json_text(),
        "{\"day\":27,\"month\":6,\"year\":2024}"
    );
    assert_eq!(parent.values.len(), 1);
    assert!(matches!(parent.values.get("date"), Some(JsonValue::Map(_))));
}

#[test]
fn json_map_rejects_duplicate_names_like_original() {
    let mut json = Date::new(2024, 6, 27).to_json(None, None);
    let duplicate = std::panic::catch_unwind(move || {
        json.add_int("year", 2025);
    });
    assert!(duplicate.is_err());

    let mut parent = JsonContainer::default();
    let date = Date::new(2024, 6, 27).to_json(None, None);
    parent.add_map("date", date.clone());
    let duplicate = std::panic::catch_unwind(move || {
        parent.add_map("date", date);
    });
    assert!(duplicate.is_err());
}

#[test]
fn json_map_at_and_get_keys_match_original_json_map_surface() {
    let mut parent = JsonContainer::default();
    let json = Date::new(2024, 6, 27).to_json(Some(&mut parent), Some("date"));

    assert_eq!(
        json.get_keys(),
        vec!["day".to_string(), "month".to_string(), "year".to_string()]
    );
    assert!(matches!(json.at("year"), Some(JsonValue::Int(2024))));
    assert!(json.at("missing").is_none());

    assert!(matches!(parent.at("date"), Some(JsonValue::Map(_))));
    assert!(parent.at("missing").is_none());
}

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

#[test]
fn json_map_and_container_accept_scalar_and_map_children_like_original_container() {
    let mut date = JsonMap::new(None);
    date.add_int("year", 2024);
    date.add_int("month", 6);
    date.add_int("day", 27);

    let mut record = JsonMap::new(None);
    record.add_int("id", 7);
    record.add_map("date", date);

    let mut parent = JsonContainer::default();
    parent.add_int("count", 1);
    parent.add_map("record", record);

    assert_eq!(
        parent.get_keys(),
        vec!["count".to_string(), "record".to_string()]
    );
    assert_eq!(
        parent.to_json_text(),
        "{\"count\":1,\"record\":{\"date\":{\"day\":27,\"month\":6,\"year\":2024},\"id\":7}}"
    );

    let duplicate = std::panic::catch_unwind(move || {
        let mut parent = parent;
        parent.add_int("count", 2);
    });
    assert!(duplicate.is_err());
}

#[test]
fn json_map_and_container_accept_string_and_boolean_children_like_original_scalars() {
    let mut record = JsonMap::new(None);
    record.add_string("name", "a\"b\\c");
    record.add_boolean("active", true);

    assert!(matches!(
        record.at("name"),
        Some(JsonValue::String(value)) if value == "a\"b\\c"
    ));
    assert!(matches!(
        record.at("active"),
        Some(JsonValue::Boolean(true))
    ));
    assert_eq!(
        record.to_json_text(),
        "{\"active\":true,\"name\":\"a\\\"b\\\\c\"}"
    );

    let mut parent = JsonContainer::default();
    parent.add_string("label", "x\ny");
    parent.add_boolean("ok", false);
    parent.add_map("record", record);

    assert_eq!(
        parent.to_json_text(),
        "{\"label\":\"x\\ny\",\"ok\":false,\"record\":{\"active\":true,\"name\":\"a\\\"b\\\\c\"}}"
    );

    let duplicate = std::panic::catch_unwind(move || {
        let mut parent = parent;
        parent.add_boolean("ok", true);
    });
    assert!(duplicate.is_err());
}

#[test]
fn json_null_double_and_array_match_original_scalar_and_array_surface() {
    let mut array = JsonArray::new(None);
    array.add_int(7);
    array.add_null();
    array.add_double(1.25, 2);
    array.add_double(f64::NAN, 2);
    array.add_string("x");
    array.add_boolean(false);

    assert_eq!(array.get_size(), 6);
    assert_eq!(array.at(0).unwrap().get_int(), 7);
    assert!(array.at(1).unwrap().get_double().is_nan());
    assert_eq!(array.at(2).unwrap().get_double(), 1.25);
    assert_eq!(array.at(4).unwrap().get_string(), "x");
    assert!(!array.at(5).unwrap().get_boolean());
    assert_eq!(array.to_json_text(), "[7,null,1.25,null,\"x\",false]");

    let mut record = JsonMap::new(None);
    record.add_null("missing");
    record.add_double("ratio", 1.25, 2);
    record.add_array("values", array);

    assert_eq!(record.at("ratio").unwrap().as_json_double(), Some(1.25));
    assert!(record.at("missing").unwrap().as_json_null().is_some());
    assert_eq!(record.at("missing").unwrap().get_size(), 0);
    assert!(record.at("missing").unwrap().at("anything").is_none());
    assert_eq!(record.at("values").unwrap().get_size(), 6);
    assert_eq!(
        record.to_json_text(),
        "{\"missing\":null,\"ratio\":1.25,\"values\":[7,null,1.25,null,\"x\",false]}"
    );
}

#[test]
fn json_getters_reject_wrong_types_like_original_as_type_surface() {
    let mut record = JsonMap::new(None);
    record.add_int("id", 1);
    record.add_string("name", "abc");

    let wrong_double = std::panic::catch_unwind(|| {
        record.at("id").unwrap().get_double();
    });
    assert!(wrong_double.is_err());

    let wrong_array = std::panic::catch_unwind(|| {
        record.at("name").unwrap().get_size();
    });
    assert!(wrong_array.is_err());

    let missing_index = std::panic::catch_unwind(|| {
        let mut array = JsonArray::new(None);
        array.add_int(1);
        array.at(1);
    });
    assert!(missing_index.is_err());
}

#[test]
fn json_map_file_constructor_parses_original_token_surface() {
    let dir = tempfile::tempdir().unwrap();
    let json_path = dir.path().join("record.json");
    fs::write(
        &json_path,
        r#"{z:[1,1.250,null,nan,true,false,"1",bare,{x:0x10}],a:"text"}"#,
    )
    .unwrap();

    let json = JsonMap::from_file(json_path.to_str().unwrap()).unwrap();
    assert_eq!(json.get_keys(), vec!["a".to_string(), "z".to_string()]);
    assert_eq!(json.at("a").unwrap().get_string(), "text");

    let values = json.at("z").unwrap().as_json_array().unwrap();
    assert_eq!(values.get_size(), 9);
    assert_eq!(values.at(0).unwrap().get_int(), 1);
    assert_eq!(values.at(1).unwrap().get_double(), 1.25);
    assert!(values.at(2).unwrap().as_json_null().is_some());
    assert!(values.at(3).unwrap().as_json_null().is_some());
    assert!(values.at(4).unwrap().get_boolean());
    assert!(!values.at(5).unwrap().get_boolean());
    assert_eq!(values.at(6).unwrap().get_string(), "1");
    assert_eq!(values.at(7).unwrap().get_string(), "bare");
    assert_eq!(values.at(8).unwrap().at("x").unwrap().get_int(), 16);
    assert_eq!(
        json.to_json_text(),
        "{\"a\":\"text\",\"z\":[1,1.250,null,null,true,false,1,\"bare\",{\"x\":16}]}"
    );
}

#[test]
fn json_map_file_constructor_accepts_identifier_digits_and_underscores_like_original_token() {
    let dir = tempfile::tempdir().unwrap();
    let json_path = dir.path().join("identifier.json");
    fs::write(&json_path, "{gene_1:abc2,_flag:true}").unwrap();

    let json = JsonMap::from_file(json_path.to_str().unwrap()).unwrap();
    assert_eq!(
        json.get_keys(),
        vec!["_flag".to_string(), "gene_1".to_string()]
    );
    assert_eq!(json.at("gene_1").unwrap().get_string(), "abc2");
    assert!(json.at("_flag").unwrap().get_boolean());
    assert_eq!(json.to_json_text(), "{\"_flag\":true,\"gene_1\":\"abc2\"}");
}

#[test]
fn json_map_file_constructor_rejects_bad_root_and_duplicate_names() {
    let dir = tempfile::tempdir().unwrap();
    let bad_root_path = dir.path().join("array.json");
    fs::write(&bad_root_path, "[1]").unwrap();
    let err = JsonMap::from_file(bad_root_path.to_str().unwrap()).unwrap_err();
    assert!(err.to_string().contains("text should start with '{'"));

    let duplicate_path = dir.path().join("duplicate.json");
    fs::write(&duplicate_path, "{a:1,a:2}").unwrap();
    let err = JsonMap::from_file(duplicate_path.to_str().unwrap()).unwrap_err();
    assert!(err
        .to_string()
        .contains("Duplicate name in JSON map: \"a\""));
}

#[test]
fn json_map_file_constructor_keeps_backslash_text_like_original_token() {
    let dir = tempfile::tempdir().unwrap();
    let json_path = dir.path().join("backslash.json");
    fs::write(
        &json_path,
        r#"{"path":"c:\\tmp","items":["one\,two","\]"]}"#,
    )
    .unwrap();

    let json = JsonMap::from_file(json_path.to_str().unwrap()).unwrap();
    assert_eq!(json.at("path").unwrap().get_string(), r#"c:\\tmp"#);

    let values = json.at("items").unwrap().as_json_array().unwrap();
    assert_eq!(values.at(0).unwrap().get_string(), r#"one\,two"#);
    assert_eq!(values.at(1).unwrap().get_string(), r#"\]"#);
    assert_eq!(
        json.to_json_text(),
        "{\"items\":[\"one\\\\,two\",\"\\\\]\"],\"path\":\"c:\\\\\\\\tmp\"}"
    );
}

#[test]
fn json_map_file_constructor_rejects_backslash_escaped_quote_like_original_token() {
    let dir = tempfile::tempdir().unwrap();
    let json_path = dir.path().join("escaped_quote.json");
    fs::write(&json_path, r#"{"a\"b":"x"}"#).unwrap();

    let err = JsonMap::from_file(json_path.to_str().unwrap()).unwrap_err();
    assert!(err.to_string().contains("':'"));
}

#[test]
fn json_map_file_constructor_uses_original_general_delimiter_token_path() {
    let dir = tempfile::tempdir().unwrap();
    let json_path = dir.path().join("semicolon.json");
    fs::write(&json_path, "{a;1}").unwrap();

    let err = JsonMap::from_file(json_path.to_str().unwrap()).unwrap_err();
    assert!(err.to_string().contains("':'"));
    assert!(!err.to_string().contains("Unknown token"));
}

#[test]
fn json_map_file_constructor_matches_original_exponent_minus_tokenization() {
    let dir = tempfile::tempdir().unwrap();
    let accepted_path = dir.path().join("positive_exponent.json");
    fs::write(&accepted_path, "{value:1e2}").unwrap();

    let json = JsonMap::from_file(accepted_path.to_str().unwrap()).unwrap();
    assert_eq!(json.at("value").unwrap().get_double(), 100.0);
    assert_eq!(json.to_json_text(), "{\"value\":1e+02}");

    let rejected_path = dir.path().join("negative_exponent.json");
    fs::write(&rejected_path, "{value:1e-2}").unwrap();

    let err = JsonMap::from_file(rejected_path.to_str().unwrap()).unwrap_err();
    assert!(err.to_string().contains("','"));
}

#[test]
fn json_map_file_constructor_preserves_lowercase_scientific_double_format_flag() {
    let dir = tempfile::tempdir().unwrap();
    let json_path = dir.path().join("scientific.json");
    fs::write(&json_path, "{lower:1.25e2,upper:1.25E2}").unwrap();

    let json = JsonMap::from_file(json_path.to_str().unwrap()).unwrap();

    assert_eq!(json.to_json_text(), "{\"lower\":1.25e+02,\"upper\":125.00}");
}

#[test]
fn json_text_output_matches_original_to_str_for_natural_strings() {
    let mut record = JsonMap::new(None);
    record.add_string("1", "123");
    record.add_string("code", "0123");
    record.add_string("line", "a\rb");

    assert_eq!(
        record.to_json_text(),
        "{1:123,\"code\":\"0123\",\"line\":\"a\\tb\"}"
    );
}

#[test]
fn text_table_file_constructor_reuses_pound_header_and_pads_short_rows() {
    let dir = tempfile::tempdir().unwrap();
    let table_path = dir.path().join("table.tsv");
    fs::write(&table_path, "#old\tignored\n#main\tvalue\na\t1\nb\n").unwrap();

    let table = TextTable::from_file(table_path.to_str().unwrap()).unwrap();
    assert!(table.pound);
    assert_eq!(
        table
            .header
            .iter()
            .map(|h| h.name.as_str())
            .collect::<Vec<_>>(),
        vec!["main", "value"]
    );
    assert_eq!(
        table.rows,
        vec![
            vec!["a".to_string(), "1".to_string()],
            vec!["b".to_string(), String::new()],
        ]
    );
}

#[test]
fn text_table_file_constructor_trims_header_fields_like_string_vector() {
    let dir = tempfile::tempdir().unwrap();
    let table_path = dir.path().join("table.tsv");
    fs::write(&table_path, " gene \t value \nbla\t1\n").unwrap();

    let table = TextTable::from_file(table_path.to_str().unwrap()).unwrap();
    assert_eq!(
        table
            .header
            .iter()
            .map(|h| h.name.as_str())
            .collect::<Vec<_>>(),
        vec!["gene", "value"]
    );
    assert_eq!(table.col2num("gene").unwrap(), 0);
    assert_eq!(table.col2num("value").unwrap(), 1);
}

#[test]
fn text_table_file_constructor_headerless_mode_names_columns_from_first_row() {
    let dir = tempfile::tempdir().unwrap();
    let table_path = dir.path().join("no_header.tsv");
    fs::write(&table_path, "alpha\tbeta\ncharlie\n").unwrap();

    let table = TextTable::from_file_with_header(table_path.to_str().unwrap(), false).unwrap();
    assert!(table.pound);
    assert_eq!(
        table
            .header
            .iter()
            .map(|h| h.name.as_str())
            .collect::<Vec<_>>(),
        vec!["1", "2"]
    );
    assert_eq!(
        table.rows,
        vec![
            vec!["alpha".to_string(), "beta".to_string()],
            vec!["charlie".to_string(), String::new()],
        ]
    );
}

#[test]
fn text_table_file_constructor_applies_column_synonyms() {
    let dir = tempfile::tempdir().unwrap();
    let table_path = dir.path().join("table.tsv");
    let synonyms_path = dir.path().join("synonyms.tsv");
    fs::write(&table_path, "legacy\tvalue\nx\t1\n").unwrap();
    fs::write(&synonyms_path, "canonical\nlegacy\n\nother\nunused\n").unwrap();

    let table = TextTable::from_file_with_synonyms(
        table_path.to_str().unwrap(),
        synonyms_path.to_str().unwrap(),
    )
    .unwrap();
    assert_eq!(table.col2num("canonical").unwrap(), 0);
    assert_eq!(table.col2num_("legacy"), None);
}

#[test]
fn set_header_uses_common_nulls_and_second_numeric_pass() {
    let mut table = TextTable::with_header(vec![Header::new("n"), Header::new("v")]);
    table.header[0].len_max = 8;
    table.rows = vec![
        vec!["1".to_string(), "NA".to_string()],
        vec!["2.125".to_string(), "x".to_string()],
    ];
    table.set_header().unwrap();
    assert!(table.header[0].numeric);
    assert_eq!(table.header[0].decimals, 3);
    assert_eq!(table.header[0].len_max, 8);
    assert!(table.header[1].null);
}

#[test]
fn set_header_preserves_existing_header_metadata_like_original() {
    let mut text = Header::new("text");
    text.numeric = false;
    text.null = true;
    text.choices.insert("old".to_string());

    let mut number = Header::new("number");
    number.scientific = true;
    number.decimals = 4;

    let mut table = TextTable::with_header(vec![text, number]);
    table.rows = vec![vec!["123".to_string(), "1.2".to_string()]];
    table.set_header().unwrap();

    assert!(!table.header[0].numeric);
    assert!(table.header[0].null);
    assert_eq!(
        table.header[0].choices.iter().cloned().collect::<Vec<_>>(),
        vec!["123".to_string(), "old".to_string()]
    );
    assert!(table.header[1].scientific);
    assert_eq!(table.header[1].decimals, 4);
    assert_eq!(table.header[1].len_max, 6);
}

#[test]
fn set_header_and_sort_use_original_strtod_numeric_surface() {
    let mut table = TextTable::with_header(vec![Header::new("value")]);
    table.rows = vec![vec!["0x10".to_string()], vec!["2".to_string()]];
    table.set_header().unwrap();

    assert!(table.header[0].numeric);
    table.sort(&["value"]).unwrap();

    assert_eq!(
        table.rows,
        vec![vec!["2".to_string()], vec!["0x10".to_string()]]
    );
}

#[test]
fn sort_preserves_existing_header_metadata_like_original() {
    let mut table = TextTable::with_header(vec![Header::new("k")]);
    table.rows = vec![
        vec![String::new()],
        vec!["x".to_string()],
        vec!["12345".to_string()],
    ];
    table.set_header().unwrap();
    assert!(table.header[0].null);
    assert!(!table.header[0].numeric);
    assert_eq!(table.header[0].len_max, 5);
    assert_eq!(
        table.header[0].choices.iter().cloned().collect::<Vec<_>>(),
        vec!["12345".to_string(), "x".to_string()]
    );

    table.rows = vec![vec!["b".to_string()], vec!["a".to_string()]];
    table.sort(&["k"]).unwrap();

    assert_eq!(
        table.rows,
        vec![vec!["a".to_string()], vec!["b".to_string()]]
    );
    assert!(table.header[0].null);
    assert!(!table.header[0].numeric);
    assert_eq!(table.header[0].len_max, 5);
    assert_eq!(
        table.header[0].choices.iter().cloned().collect::<Vec<_>>(),
        vec!["12345".to_string(), "x".to_string()]
    );
}

#[test]
fn deredundify_preserves_existing_header_metadata_like_original() {
    let mut table = TextTable::with_header(vec![Header::new("k"), Header::new("score")]);
    table.rows = vec![
        vec!["x".to_string(), String::new()],
        vec!["y".to_string(), "12345".to_string()],
    ];
    table.set_header().unwrap();
    assert!(table.header[1].null);
    assert!(table.header[1].numeric);
    assert_eq!(table.header[1].len_max, 6);
    assert_eq!(
        table.header[1].choices.iter().cloned().collect::<Vec<_>>(),
        vec!["12345".to_string()]
    );

    table.rows = vec![
        vec!["b".to_string(), "2".to_string()],
        vec!["b".to_string(), "3".to_string()],
        vec!["a".to_string(), "1".to_string()],
    ];
    table
        .deredundify(&["k"], |better, worse| {
            (better[1].parse::<i32>().unwrap() > worse[1].parse::<i32>().unwrap()) as i32
        })
        .unwrap();

    assert_eq!(
        table.rows,
        vec![
            vec!["a".to_string(), "1".to_string()],
            vec!["b".to_string(), "3".to_string()],
        ]
    );
    assert!(table.header[1].null);
    assert!(table.header[1].numeric);
    assert_eq!(table.header[1].len_max, 6);
    assert_eq!(
        table.header[1].choices.iter().cloned().collect::<Vec<_>>(),
        vec!["12345".to_string()]
    );
}

#[test]
fn text_table_find_values_and_index_match_original_shape() {
    let mut table = TextTable::with_header(vec![
        Header::new("sample"),
        Header::new("gene"),
        Header::new("value"),
    ]);
    table.rows = vec![
        vec!["s1".to_string(), "bla".to_string(), "2".to_string()],
        vec!["s2".to_string(), "aac".to_string(), "1".to_string()],
        vec!["s3".to_string(), "bla".to_string(), "3".to_string()],
        vec!["s4".to_string(), "".to_string(), "4".to_string()],
    ];

    let gene_col = table.col2num("gene").unwrap();
    assert_eq!(
        table.col_nums_row2values(&[0, gene_col], 2),
        vec!["s3".to_string(), "bla".to_string()]
    );
    assert_eq!(table.find(&[gene_col], &["bla".to_string()], 1), Some(2));
    assert_eq!(
        table.col2values(gene_col),
        vec!["aac".to_string(), "bla".to_string()]
    );

    let index = table.build_index(&["gene"]).unwrap();
    assert_eq!(
        index.find(&["bla".to_string()]).unwrap(),
        &vec![0usize, 2usize]
    );
    assert!(index.find(&["missing".to_string()]).is_none());
}

#[test]
fn text_table_key_constructor_errors_keep_original_context() {
    let mut table = TextTable::with_header(vec![Header::new("sample"), Header::new("gene")]);
    table.name = "key.tsv".to_string();
    table.rows = vec![
        vec!["s1".to_string(), "bla".to_string()],
        vec!["s1".to_string(), "bla".to_string()],
    ];

    let err = match table.build_key(&["sample", "gene"]) {
        Ok(_) => panic!("duplicate key should fail"),
        Err(err) => err,
    };
    assert!(err
        .to_string()
        .contains("Duplicate key s1,bla for the key on sample,gene"));
    assert!(err.to_string().contains("In table file: key.tsv"));

    table.rows[1][1].clear();
    let err = match table.build_key(&["sample", "gene"]) {
        Ok(_) => panic!("empty key should fail"),
        Err(err) => err,
    };
    assert!(err.to_string().contains("Empty value in key, in row 2"));
    assert!(err.to_string().contains("In table file: key.tsv"));
}

#[test]
fn text_table_qc_errors_keep_original_context() {
    let mut duplicate = TextTable::with_header(vec![Header::new("gene"), Header::new("gene")]);
    duplicate.name = "qc.tsv".to_string();
    let err = duplicate.qc().unwrap_err();
    assert!(err.to_string().contains("Duplicate column name: \"gene\""));
    assert!(err.to_string().contains("In table file: qc.tsv"));

    let mut wrong_width = TextTable::with_header(vec![Header::new("gene")]);
    wrong_width.name = "qc.tsv".to_string();
    wrong_width.rows = vec![vec!["bla".to_string(), "extra".to_string()]];
    let err = wrong_width.qc().unwrap_err();
    assert!(err
        .to_string()
        .contains("Row 1 contains 2 columns whereas table has 1 columns"));
    assert!(err.to_string().contains("In table file: qc.tsv"));

    let mut tab = TextTable::with_header(vec![Header::new("gene")]);
    tab.name = "qc.tsv".to_string();
    tab.rows = vec![vec!["bla\tTEM".to_string()]];
    let err = tab.qc().unwrap_err();
    assert!(err
        .to_string()
        .contains("Field \"gene\" of row 1 contains a tab character"));
    assert!(err.to_string().contains("In table file: qc.tsv"));

    let mut eol = TextTable::with_header(vec![Header::new("gene")]);
    eol.name = "qc.tsv".to_string();
    eol.rows = vec![vec!["bla\nTEM".to_string()]];
    let err = eol.qc().unwrap_err();
    assert!(err
        .to_string()
        .contains("Field \"gene\" of row 1 contains an EOL character"));
    assert!(err.to_string().contains("In table file: qc.tsv"));
}

#[test]
fn text_table_qc_field_diagnostic_uses_original_row_index_header_name() {
    let mut table = TextTable::with_header(vec![Header::new("first"), Header::new("second")]);
    table.name = "qc.tsv".to_string();
    table.rows = vec![
        vec!["ok".to_string(), "ok".to_string()],
        vec!["bad\tfield".to_string(), "ok".to_string()],
    ];

    let err = table.qc().unwrap_err();
    assert!(err
        .to_string()
        .contains("Field \"second\" of row 2 contains a tab character"));
    assert!(err.to_string().contains("In table file: qc.tsv"));
}

#[test]
fn text_table_group_sum_min_max_and_aggregate() {
    let mut table = TextTable::with_header(vec![
        Header::new("sample"),
        Header::new("count"),
        Header::new("first"),
        Header::new("last"),
        Header::new("gene"),
        Header::new("drop"),
    ]);
    table.rows = vec![
        vec!["b", "2.50", "4", "4", "tet", "x"],
        vec!["a", "1.25", "3", "3", "bla", "x"],
        vec!["a", "2.50", "2", "5", "aac", "x"],
        vec!["b", "", "1", "7", "tet", "x"],
    ]
    .into_iter()
    .map(|row| row.into_iter().map(str::to_string).collect())
    .collect();
    table.set_header().unwrap();

    table
        .group(&["sample"], &["count"], &["first"], &["last"], &["gene"])
        .unwrap();

    assert_eq!(
        table
            .header
            .iter()
            .map(|h| h.name.as_str())
            .collect::<Vec<_>>(),
        vec!["sample", "count", "first", "last", "gene"]
    );
    assert_eq!(
        table.rows,
        vec![
            vec!["a", "3.75", "2", "5", "aac,bla"],
            vec!["b", "2.50", "1", "7", "tet"],
        ]
    );
}

#[test]
fn text_table_group_rejects_overlapping_operations_and_commas() {
    let mut table = TextTable::with_header(vec![Header::new("k"), Header::new("v")]);
    table.rows = vec![vec!["a".to_string(), "x".to_string()]];
    let err = table.group(&["k"], &[], &[], &[], &["k"]).unwrap_err();
    assert!(err
        .to_string()
        .contains("operations \"group by\" and \"aggregation\""));

    let mut table = TextTable::with_header(vec![Header::new("k"), Header::new("v")]);
    table.rows = vec![
        vec!["a".to_string(), "x".to_string()],
        vec!["a".to_string(), "y,z".to_string()],
    ];
    table.set_header().unwrap();
    let err = table.group(&["k"], &[], &[], &[], &["v"]).unwrap_err();
    assert!(err
        .to_string()
        .contains("Cannot aggregate column v for row 2"));
}

#[test]
fn text_table_is_key_and_substitue_column_match_original_names() {
    let mut table = TextTable::with_header(vec![Header::new("id"), Header::new("value")]);
    table.rows = vec![
        vec!["1".to_string(), "x".to_string()],
        vec!["2".to_string(), "x".to_string()],
    ];
    table.set_header().unwrap();
    assert!(table.is_key(table.col2num("id").unwrap()));
    assert!(!table.is_key(table.col2num("value").unwrap()));

    let mut from = "value".to_string();
    table.substitue_column(&mut from, "copy").unwrap();
    assert_eq!(from, "copy");
    assert_eq!(table.col2num("copy").unwrap(), 2);
    let err = table.duplicate_column("id", "copy").unwrap_err();
    assert!(err
        .to_string()
        .contains("Table already has column \"copy\""));
    assert_eq!(TextTable::aggr2values(""), Vec::<String>::new());
    assert_eq!(TextTable::aggr2values("z,a,z"), vec!["a", "z"]);
}
