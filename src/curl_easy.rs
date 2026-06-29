use std::env;
use std::ffi::{CStr, CString};
use std::fmt;
use std::fs::{self, File};
use std::io::Write;
use std::os::raw::{c_char, c_int, c_long, c_uint, c_void};
use std::ptr;
use std::slice;

use anyhow::{bail, Context, Result};

const CURL_ERROR_SIZE: usize = 256;
const CURLVERSION_NOW: c_int = 9;
const CURLOPT_WRITEDATA: c_int = 10001;
const CURLOPT_URL: c_int = 10002;
const CURLOPT_ERRORBUFFER: c_int = 10010;
const CURLOPT_WRITEFUNCTION: c_int = 20011;
const CURLOPT_CAINFO: c_int = 10065;
const CURLOPT_HTTP_VERSION: c_int = 84;
const CURLOPT_FTP_USE_EPSV: c_int = 85;
const CURL_HTTP_VERSION_1_0: c_long = 1;

#[allow(non_camel_case_types)]
type CURL = c_void;

type WriteCallback = unsafe extern "C" fn(*mut c_char, usize, usize, *mut c_void) -> usize;

#[repr(C)]
struct CurlVersionInfoData {
    age: c_int,
    version: *const c_char,
    version_num: c_uint,
}

#[link(name = "curl")]
unsafe extern "C" {
    fn curl_version_info(age: c_int) -> *const CurlVersionInfoData;
    fn curl_easy_init() -> *mut CURL;
    fn curl_easy_cleanup(curl: *mut CURL);
    fn curl_easy_perform(curl: *mut CURL) -> c_int;
    fn curl_easy_setopt(curl: *mut CURL, option: c_int, ...) -> c_int;
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SoftwareVersion {
    pub major: u32,
    pub minor: u32,
    pub patch: u32,
}

impl SoftwareVersion {
    pub fn str(&self) -> String {
        format!("{}.{}.{}", self.major, self.minor, self.patch)
    }
}

impl fmt::Display for SoftwareVersion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.str())
    }
}

pub fn get_lib_version() -> SoftwareVersion {
    let ver = unsafe { curl_version_info(CURLVERSION_NOW) };
    if ver.is_null() {
        return SoftwareVersion::default();
    }

    let version_num = unsafe { (*ver).version_num };
    SoftwareVersion {
        major: (version_num >> 16) & 0xff,
        minor: (version_num >> 8) & 0xff,
        patch: version_num & 0xff,
    }
}

unsafe extern "C" fn write_stream_cb(
    ptr: *mut c_char,
    size: usize,
    n_memb: usize,
    user_data: *mut c_void,
) -> usize {
    assert!(!ptr.is_null());
    assert!(size == 1);
    assert!(!user_data.is_null());

    let f = unsafe { &mut *(user_data as *mut File) };
    let bytes = unsafe { slice::from_raw_parts(ptr as *const u8, n_memb) };
    for c in bytes {
        if f.write_all(&[*c]).is_err() {
            return 0;
        }
    }

    n_memb
}

unsafe extern "C" fn write_string_cb(
    ptr: *mut c_char,
    size: usize,
    n_memb: usize,
    user_data: *mut c_void,
) -> usize {
    assert!(!ptr.is_null());
    assert!(size == 1);
    assert!(!user_data.is_null());

    let s = unsafe { &mut *(user_data as *mut Vec<u8>) };
    let bytes = unsafe { slice::from_raw_parts(ptr as *const u8, n_memb) };
    s.extend_from_slice(bytes);

    n_memb
}

#[derive(Debug)]
pub struct Curl {
    eh: *mut CURL,
    ca_info: Option<CString>,
}

impl Curl {
    pub fn new() -> Result<Self> {
        let eh = unsafe { curl_easy_init() };
        if eh.is_null() {
            bail!("Cannot initialize curl_easy");
        }

        let ca_info = env::var("CURL_CA_BUNDLE")
            .ok()
            .map(CString::new)
            .transpose()
            .context("CURL_CA_BUNDLE contains an interior NUL byte")?;

        if let Some(env_ca_bundle) = ca_info.as_ref() {
            unsafe {
                curl_easy_setopt(eh, CURLOPT_CAINFO, env_ca_bundle.as_ptr());
            }
        }
        unsafe {
            curl_easy_setopt(eh, CURLOPT_HTTP_VERSION, CURL_HTTP_VERSION_1_0);
        }

        Ok(Self { eh, ca_info })
    }

    pub fn download(&mut self, url: &str, f_name: &str) -> Result<()> {
        assert!(!f_name.is_empty());

        {
            let mut f =
                File::create(f_name).with_context(|| format!("failed to create {f_name}"))?;
            unsafe {
                curl_easy_setopt(
                    self.eh,
                    CURLOPT_WRITEFUNCTION,
                    write_stream_cb as WriteCallback,
                );
                curl_easy_setopt(self.eh, CURLOPT_WRITEDATA, &mut f as *mut File);
            }
            self.process(url, "download")?;
        }

        let bytes = fs::read(f_name).with_context(|| format!("failed to read {f_name}"))?;
        if bytes.split(|b| b.is_ascii_whitespace()).next() == Some(&b"<?xml"[..]) {
            bail!("Cannot download \"{}\"", f_name);
        }

        Ok(())
    }

    pub fn read(&mut self, url: &str) -> Result<String> {
        let mut s = Vec::with_capacity(1024);
        unsafe {
            curl_easy_setopt(
                self.eh,
                CURLOPT_WRITEFUNCTION,
                write_string_cb as WriteCallback,
            );
            curl_easy_setopt(self.eh, CURLOPT_WRITEDATA, &mut s as *mut Vec<u8>);
        }
        self.process(url, "read")?;

        String::from_utf8(s).context("CURL: response is not valid UTF-8")
    }

    fn process(&mut self, url: &str, error_msg_action: &str) -> Result<()> {
        assert!(!url.is_empty());

        let url_c = CString::new(url).context("CURL URL contains an interior NUL byte")?;
        let mut err = [0 as c_char; CURL_ERROR_SIZE + 1];

        unsafe {
            curl_easy_setopt(self.eh, CURLOPT_ERRORBUFFER, err.as_mut_ptr());
            curl_easy_setopt(self.eh, CURLOPT_URL, url_c.as_ptr());
            if url.starts_with("ftp://") {
                curl_easy_setopt(self.eh, CURLOPT_FTP_USE_EPSV, 0 as c_long);
            }
        }

        let cc = unsafe { curl_easy_perform(self.eh) };
        if cc != 0 {
            let err = unsafe { CStr::from_ptr(err.as_ptr()) }.to_string_lossy();
            bail!(
                "CURL: Cannot {error_msg_action}\n  from {url}\n  code={cc}\n  error: {err}\n  version: {}",
                get_lib_version()
            );
        }

        let _ = self.ca_info.as_ref();
        Ok(())
    }
}

impl Drop for Curl {
    fn drop(&mut self) {
        if !self.eh.is_null() {
            unsafe {
                curl_easy_cleanup(self.eh);
            }
            self.eh = ptr::null_mut();
        }
    }
}
