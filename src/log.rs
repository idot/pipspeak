use std::{
    fs::File,
    io::{BufWriter, Write},
};

use anyhow::Result;
use hashbrown::HashSet;
use serde::Serialize;

use crate::config::Config;

#[derive(Debug, Default, Serialize, Clone)]
pub struct Statistics {
    pub total_reads: usize,
    pub passing_reads: usize,
    pub fraction_passing: f64,
    pub whitelist_size: usize,
    pub num_filtered_1: usize,
    pub num_filtered_2: usize,
    pub num_filtered_3: usize,
    pub num_filtered_4: usize,
    pub num_filtered_umi: usize,
    #[serde(skip)]
    pub whitelist: HashSet<Vec<u8>>,
    #[serde(skip)]
    pub counter_maps: BarcodePartCounterMaps,
}
impl Statistics {
    pub fn new() -> Self {
        Self {
            counter_maps: BarcodePartCounterMaps::new(),
            ..Self::default()
        }
    }
    pub fn calculate_metrics(&mut self) {
        self.fraction_passing = self.passing_reads as f64 / self.total_reads as f64;
        self.whitelist_size = self.whitelist.len();
    }
    pub fn whitelist_to_file(&self, file: &str) -> Result<()> {
        let mut writer = File::create(file).map(BufWriter::new)?;
        for seq in &self.whitelist {
            writer.write(seq)?;
            writer.write(b"\n")?;
        }
        Ok(())
    }
    pub fn counter_maps_to_file(&self, file: &str, config: &Config) -> Result<()> {
        let mut writer = File::create(file).map(BufWriter::new)?;
        let _ = writer.write(b"position\tbarcode\tcount\n");
        for (position, map) in self.counter_maps.maps.iter().enumerate() {
            let map = map.lock().unwrap();
            for (k, v) in map.iter() {
                let bc = config.get_barcode(*k, position)
                    .and_then(|bc| String::from_utf8(bc.to_vec()).ok())
                    .unwrap_or_else(|| "unknown".to_string());
                writer.write_all(format!("{}\t{}\t{}\n", position, bc, v).as_bytes())?;
            }
        }
        Ok(())
    }
}

#[derive(Debug, Serialize)]
pub struct Timing {
    pub timestamp: String,
    pub elapsed_time: f64,
}

#[derive(Debug, Serialize)]
pub struct FileIO {
    pub readpath_r1: String,
    pub readpath_r2: String,
    pub writepath_r1: String,
    pub writepath_r2: String,
    pub whitelist_path: String,
}

#[derive(Debug, Serialize)]
pub struct Parameters {
    pub offset: usize,
    pub umi_len: usize,
    pub exact_matching: bool,
    pub write_linkers: bool,
    pub pipspeak_version: String,
}

#[derive(Debug, Serialize)]
/// A struct to hold the information about the run
pub struct Log {
    pub parameters: Parameters,
    pub file_io: FileIO,
    pub statistics: Statistics,
    pub timing: Timing,
}
impl Log {
    pub fn stderr(&self) -> Result<()> {
        let yaml = serde_yaml::to_string(&self)?;
        eprint!("{}", yaml);
        Ok(())
    }

    pub fn to_file(&self, path: &str) -> Result<()> {
        let yaml = serde_yaml::to_string(&self)?;
        std::fs::write(path, yaml)?;
        Ok(())
    }
}


use std::collections::HashMap;
use std::sync::Mutex;

#[derive(Debug, Default, Serialize)]
pub struct BarcodePartCounterMaps {
    maps: Vec<Mutex<HashMap<usize, usize>>>,
}

impl Clone for BarcodePartCounterMaps {
    fn clone(&self) -> Self {
        let maps = self.maps.iter().map(|m| {
            let map = m.lock().unwrap();
            Mutex::new(map.clone())
        }).collect();
        Self { maps }
    }
}

impl BarcodePartCounterMaps {
    // Initialize the counter maps
    pub fn new() -> Self {
        let maps = vec![
            Mutex::new(HashMap::new()),
            Mutex::new(HashMap::new()),
            Mutex::new(HashMap::new()),
            Mutex::new(HashMap::new()),
        ];
        Self { maps }
    }

    /// Add to the respective map
    pub fn add(&self, index: usize, position: usize) {
        let mut map = self.maps[position].lock().unwrap();
        *map.entry(index).or_insert(0) += 1;
    }
}