use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use crate::utils::DynError;
use roxmltree::Document;

#[derive(Debug, Clone)]
pub struct IFileData {
    pub ra: String,
    pub dec: String,
    pub epoch: Option<String>,
    pub source: Option<String>,
    pub fft: Option<usize>,
    pub sampling_hz: Option<f64>,
    pub bit: Option<usize>,
    pub level: Option<String>,
    pub shuffle: Option<String>,
    pub obsfreq_mhz: Option<f64>,
    pub clock_delay_s: Option<f64>,
    pub clock_rate_sps: Option<f64>,
    pub sideband: Option<String>,
    pub ant1_rotation_hz: Option<f64>,
    pub ant2_rotation_hz: Option<f64>,
    pub ant1_center_mhz: Option<f64>,
    pub ant2_center_mhz: Option<f64>,
    pub ant1_bw_mhz: Option<f64>,
    pub ant2_bw_mhz: Option<f64>,
    pub ant1_station_name: Option<String>,
    pub ant2_station_name: Option<String>,
    pub ant1_station_key: Option<String>,
    pub ant2_station_key: Option<String>,
    pub process_epochs: Vec<String>,
}

fn parse_optional_f64(params: &HashMap<String, String>, keys: &[&str]) -> Result<Option<f64>, DynError> {
    for key in keys {
        if let Some(value) = params.get(*key) {
            return Ok(Some(value.trim().parse::<f64>()?));
        }
    }
    Ok(None)
}

fn parse_clock_pair(clock_raw: &str) -> Result<(Option<f64>, Option<f64>), DynError> {
    let parts: Vec<&str> = clock_raw
        .split(|c: char| c == ',' || c == ';' || c.is_whitespace())
        .filter(|s| !s.is_empty())
        .collect();
    if parts.is_empty() {
        return Ok((None, None));
    }
    let delay = Some(parts[0].parse::<f64>()?);
    let rate = if parts.len() >= 2 {
        Some(parts[1].parse::<f64>()?)
    } else {
        None
    };
    Ok((delay, rate))
}

fn parse_ifile_kv(path: &PathBuf) -> Result<IFileData, DynError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut params = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        // Accept both full-line and trailing '#' comments in .ifile.
        let line = line.splitn(2, '#').next().unwrap_or("").trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with(';') {
            continue;
        }
        if let Some(index) = line.find('=') {
            let (key, value) = line.split_at(index);
            let key = key.trim().to_ascii_lowercase().replace('_', "");
            let value = value
                .trim_start_matches('=')
                .trim()
                .trim_matches('"')
                .trim_matches('\'')
                .to_string();
            params.insert(key, value);
        }
    }

    let mut clock_delay_s = parse_optional_f64(&params, &["clockdelay", "clockoffset"])?;
    let mut clock_rate_sps = parse_optional_f64(&params, &["clockrate", "clockdrift"])?;
    let obsfreq_mhz = parse_optional_f64(
        &params,
        &[
            "obsfreq",
            "obsfreqmhz",
            "skyfreq",
            "observingfrequencymhz",
            "observingfrequency",
            "freq",
        ],
    )?;
    if let Some(clock_raw) = params.get("clock") {
        let (delay_from_clock, rate_from_clock) = parse_clock_pair(clock_raw)?;
        if clock_delay_s.is_none() {
            clock_delay_s = delay_from_clock;
        }
        if clock_rate_sps.is_none() {
            clock_rate_sps = rate_from_clock;
        }
    }

    let ra = params
        .get("ra")
        .or_else(|| params.get("srcra"))
        .or_else(|| params.get("rightascension"))
        .ok_or("ra not found in ifile")?
        .clone();
    let dec = params
        .get("dec")
        .or_else(|| params.get("srcdec"))
        .or_else(|| params.get("declination"))
        .ok_or("dec not found in ifile")?
        .clone();

    Ok(IFileData {
        ra,
        dec,
        epoch: params.get("epoch").or_else(|| params.get("scanepoch")).cloned(),
        source: params
            .get("source")
            .or_else(|| params.get("sourcename"))
            .or_else(|| params.get("object"))
            .cloned(),
        fft: params.get("fft").map(|v| v.trim().parse()).transpose()?,
        sampling_hz: params
            .get("samplinghz")
            .or_else(|| params.get("sampling"))
            .or_else(|| params.get("fs"))
            .map(|v| v.trim().parse())
            .transpose()?,
        bit: params.get("bit").or_else(|| params.get("bits")).map(|v| v.trim().parse()).transpose()?,
        level: params.get("level").or_else(|| params.get("levels")).cloned(),
        shuffle: params.get("shuffle").cloned(),
        obsfreq_mhz,
        clock_delay_s,
        clock_rate_sps,
        sideband: params.get("sideband").cloned(),
        ant1_rotation_hz: parse_optional_f64(
            &params,
            &["ant1rotationhz", "rotation1hz", "rotationhz", "rotation"],
        )?,
        ant2_rotation_hz: parse_optional_f64(
            &params,
            &["ant2rotationhz", "rotation2hz", "rotationhz", "rotation"],
        )?,
        ant1_center_mhz: parse_optional_f64(
            &params,
            &["ant1centermhz", "center1mhz", "obsfreq1mhz", "freq1mhz"],
        )?,
        ant2_center_mhz: parse_optional_f64(
            &params,
            &["ant2centermhz", "center2mhz", "obsfreq2mhz", "freq2mhz"],
        )?,
        ant1_bw_mhz: parse_optional_f64(&params, &["ant1bwmhz", "bw1mhz"])?,
        ant2_bw_mhz: parse_optional_f64(&params, &["ant2bwmhz", "bw2mhz"])?,
        ant1_station_name: None,
        ant2_station_name: None,
        ant1_station_key: None,
        ant2_station_key: None,
        process_epochs: Vec::new(),
    })
}

fn text_of_child<'a>(node: roxmltree::Node<'a, 'a>, tag: &str) -> Option<&'a str> {
    node.children()
        .find(|n| n.tag_name().name().eq_ignore_ascii_case(tag))
        .and_then(|n| n.text())
        .map(|s| s.trim())
}

fn parse_f64_text(node: roxmltree::Node<'_, '_>, tag: &str) -> Result<Option<f64>, DynError> {
    if let Some(text) = text_of_child(node, tag) {
        return Ok(Some(text.parse::<f64>()?));
    }
    Ok(None)
}

fn normalize_epoch_text(input: &str) -> String {
    let trimmed = input.trim();
    // XML schedule format: YYYY/DDD HH:MM:SS -> YYYY-MM-DDTHH:MM:SSZ
    let parts: Vec<&str> = trimmed.split_whitespace().collect();
    if parts.len() >= 2 {
        if let Some((year, doy)) = parts[0].split_once('/') {
            if year.len() == 4
                && doy.len() == 3
                && year.chars().all(|c| c.is_ascii_digit())
                && doy.chars().all(|c| c.is_ascii_digit())
            {
                let year_i = year.parse::<i32>().ok();
                let doy_i = doy.parse::<u32>().ok();
                if let (Some(year_i), Some(mut day_of_year)) = (year_i, doy_i) {
                    let leap = (year_i % 4 == 0 && year_i % 100 != 0) || (year_i % 400 == 0);
                    let mut month_days = [31u32, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
                    if leap {
                        month_days[1] = 29;
                    }
                    let mut month = 1u32;
                    for &days in month_days.iter() {
                        if day_of_year <= days {
                            break;
                        }
                        day_of_year -= days;
                        month += 1;
                    }
                    if (1..=12).contains(&month) && day_of_year >= 1 {
                        let time = parts[1].trim();
                        if time.len() == 8
                            && time.chars().enumerate().all(|(idx, c)| {
                                if idx == 2 || idx == 5 {
                                    c == ':'
                                } else {
                                    c.is_ascii_digit()
                                }
                            })
                        {
                            return format!(
                                "{:04}-{:02}-{:02}T{}Z",
                                year_i, month, day_of_year, time
                            );
                        }
                    }
                }
            }
        }
    }
    trimmed.to_string()
}

fn parse_xml_schedule(path: &PathBuf) -> Result<IFileData, DynError> {
    let xml = std::fs::read_to_string(path)?;
    let doc = Document::parse(&xml)?;

    let process_nodes: Vec<roxmltree::Node<'_, '_>> = doc
        .descendants()
        .filter(|n| n.tag_name().name().eq_ignore_ascii_case("process"))
        .collect();
    if process_nodes.is_empty() {
        return Err("process not found in XML schedule".into());
    }
    let process = process_nodes[0];
    let source_name = text_of_child(process, "object")
        .or_else(|| {
            doc.descendants()
                .find(|n| n.tag_name().name().eq_ignore_ascii_case("source"))
                .and_then(|n| n.attribute("name"))
        })
        .ok_or("source/object not found in XML schedule")?
        .to_string();
    let epoch = text_of_child(process, "epoch").map(normalize_epoch_text);
    let process_epochs: Vec<String> = process_nodes
        .iter()
        .filter_map(|p| text_of_child(*p, "epoch"))
        .map(normalize_epoch_text)
        .collect();
    let stations_text = text_of_child(process, "stations").unwrap_or("");
    let station_keys: Vec<char> = stations_text
        .chars()
        .filter(|c| !c.is_whitespace())
        .take(2)
        .collect();
    let station1_key = station_keys
        .first()
        .ok_or("process stations must contain at least 2 keys")?;
    let station2_key = station_keys
        .get(1)
        .ok_or("process stations must contain at least 2 keys")?;

    let source_node = doc
        .descendants()
        .find(|n| {
            n.tag_name().name().eq_ignore_ascii_case("source")
                && n.attribute("name")
                    .map(|s| s.eq_ignore_ascii_case(&source_name))
                    .unwrap_or(false)
        })
        .or_else(|| {
            doc.descendants()
                .find(|n| n.tag_name().name().eq_ignore_ascii_case("source"))
        })
        .ok_or("source entry not found in XML schedule")?;
    let ra = text_of_child(source_node, "ra")
        .ok_or("ra not found in source entry")?
        .to_string();
    let dec = text_of_child(source_node, "dec")
        .ok_or("dec not found in source entry")?
        .to_string();

    let station1 = doc
        .descendants()
        .find(|n| {
            n.tag_name().name().eq_ignore_ascii_case("station")
                && n.attribute("key")
                    .map(|s| s.eq_ignore_ascii_case(&station1_key.to_string()))
                    .unwrap_or(false)
        })
        .ok_or("station 1 key not found in XML schedule")?;
    let station2 = doc
        .descendants()
        .find(|n| {
            n.tag_name().name().eq_ignore_ascii_case("station")
                && n.attribute("key")
                    .map(|s| s.eq_ignore_ascii_case(&station2_key.to_string()))
                    .unwrap_or(false)
        })
        .ok_or("station 2 key not found in XML schedule")?;
    let station1_name = text_of_child(station1, "name").map(|s| s.to_string());
    let station2_name = text_of_child(station2, "name").map(|s| s.to_string());

    let terminal1_name = text_of_child(station1, "terminal").ok_or("station1 terminal missing")?;
    let terminal2_name = text_of_child(station2, "terminal").ok_or("station2 terminal missing")?;
    let terminal1 = doc
        .descendants()
        .find(|n| {
            n.tag_name().name().eq_ignore_ascii_case("terminal")
                && n.attribute("name")
                    .map(|s| s.eq_ignore_ascii_case(terminal1_name))
                    .unwrap_or(false)
        })
        .ok_or("terminal for station1 not found")?;
    let terminal2 = doc
        .descendants()
        .find(|n| {
            n.tag_name().name().eq_ignore_ascii_case("terminal")
                && n.attribute("name")
                    .map(|s| s.eq_ignore_ascii_case(terminal2_name))
                    .unwrap_or(false)
        })
        .ok_or("terminal for station2 not found")?;

    let speed1_hz = parse_f64_text(terminal1, "speed")?;
    let speed2_hz = parse_f64_text(terminal2, "speed")?;
    let sampling_hz = match (speed1_hz, speed2_hz) {
        (Some(s1), Some(s2)) => {
            if (s1 - s2).abs() > 1e-3 {
                return Err("XML terminals have different sampling speeds; this build requires a common sampling rate".into());
            }
            Some(s1)
        }
        (Some(s), None) | (None, Some(s)) => Some(s),
        (None, None) => None,
    };
    let bit = text_of_child(terminal1, "bit")
        .map(|v| v.parse::<usize>())
        .transpose()?;
    let level = text_of_child(terminal1, "level").map(|s| s.to_string());

    let stream = process
        .children()
        .find(|n| n.tag_name().name().eq_ignore_ascii_case("stream"))
        .or_else(|| {
            doc.descendants()
                .find(|n| n.tag_name().name().eq_ignore_ascii_case("stream"))
        })
        .ok_or("stream not found in XML schedule")?;
    let frequency_hz = text_of_child(stream, "frequency")
        .map(|v| v.parse::<f64>())
        .transpose()?;
    let fft = text_of_child(stream, "fft")
        .map(|v| v.parse::<usize>())
        .transpose()?;
    let obsfreq_mhz = frequency_hz.map(|v| v / 1_000_000.0);

    let special1 = stream.children().find(|n| {
        n.tag_name().name().eq_ignore_ascii_case("special")
            && n.attribute("key")
                .map(|s| s.eq_ignore_ascii_case(&station1_key.to_string()))
                .unwrap_or(false)
    });
    let special2 = stream.children().find(|n| {
        n.tag_name().name().eq_ignore_ascii_case("special")
            && n.attribute("key")
                .map(|s| s.eq_ignore_ascii_case(&station2_key.to_string()))
                .unwrap_or(false)
    });
    let rotation1_hz = special1
        .and_then(|s| text_of_child(s, "rotation"))
        .map(|v| v.parse::<f64>())
        .transpose()?
        .unwrap_or(0.0);
    let rotation2_hz = special2
        .and_then(|s| text_of_child(s, "rotation"))
        .map(|v| v.parse::<f64>())
        .transpose()?
        .unwrap_or(0.0);
    let sideband = special1
        .and_then(|s| text_of_child(s, "sideband"))
        .or_else(|| special2.and_then(|s| text_of_child(s, "sideband")))
        .map(|s| s.to_ascii_lowercase());

    // stream frequency is treated as base observing frequency; rotation is per-station offset in Hz.
    let ant1_center_mhz = None;
    let ant2_center_mhz = None;
    let ant1_bw_mhz = speed1_hz.map(|s| s / 2.0 / 1_000_000.0);
    let ant2_bw_mhz = speed2_hz.map(|s| s / 2.0 / 1_000_000.0);

    let clock1 = doc.descendants().find(|n| {
        n.tag_name().name().eq_ignore_ascii_case("clock")
            && n.attribute("key")
                .map(|s| s.eq_ignore_ascii_case(&station1_key.to_string()))
                .unwrap_or(false)
    });
    let clock2 = doc.descendants().find(|n| {
        n.tag_name().name().eq_ignore_ascii_case("clock")
            && n.attribute("key")
                .map(|s| s.eq_ignore_ascii_case(&station2_key.to_string()))
                .unwrap_or(false)
    });
    let clock1_delay = clock1
        .and_then(|c| text_of_child(c, "delay"))
        .map(|v| v.parse::<f64>())
        .transpose()?;
    let clock2_delay = clock2
        .and_then(|c| text_of_child(c, "delay"))
        .map(|v| v.parse::<f64>())
        .transpose()?;
    let clock1_rate = clock1
        .and_then(|c| text_of_child(c, "rate"))
        .map(|v| v.parse::<f64>())
        .transpose()?;
    let clock2_rate = clock2
        .and_then(|c| text_of_child(c, "rate"))
        .map(|v| v.parse::<f64>())
        .transpose()?;
    let clock_delay_s = match (clock1_delay, clock2_delay) {
        (Some(d1), Some(d2)) => Some(d2 - d1),
        _ => None,
    };
    let clock_rate_sps = match (clock1_rate, clock2_rate) {
        (Some(r1), Some(r2)) => Some(r2 - r1),
        _ => None,
    };

    Ok(IFileData {
        ra,
        dec,
        epoch,
        source: Some(source_name),
        fft,
        sampling_hz,
        bit,
        level,
        shuffle: None,
        obsfreq_mhz,
        clock_delay_s,
        clock_rate_sps,
        sideband,
        ant1_rotation_hz: Some(rotation1_hz),
        ant2_rotation_hz: Some(rotation2_hz),
        ant1_center_mhz,
        ant2_center_mhz,
        ant1_bw_mhz,
        ant2_bw_mhz,
        ant1_station_name: station1_name,
        ant2_station_name: station2_name,
        ant1_station_key: Some(station1_key.to_string()),
        ant2_station_key: Some(station2_key.to_string()),
        process_epochs,
    })
}

pub fn parse_ifile(path: &PathBuf) -> Result<IFileData, DynError> {
    let is_xml = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.eq_ignore_ascii_case("xml"))
        .unwrap_or(false);
    if is_xml {
        return parse_xml_schedule(path);
    }
    parse_ifile_kv(path)
}
