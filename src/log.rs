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
    #[serde(skip)]
    pub barcode_umi_counter: BarcodeUmiCounter,
    #[serde(skip)]
    pub umi_base_composition: UMIBaseComposition,
}
impl Statistics {
    pub fn new() -> Self {
        Self {
            counter_maps: BarcodePartCounterMaps::new(),
            barcode_umi_counter: BarcodeUmiCounter::new(),
            umi_base_composition: UMIBaseComposition::new(16),
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
    pub fn barcode_umi_stats_to_file(&self, file: &str) -> std::io::Result<()> {
        self.barcode_umi_counter.write_barcode_stats(file)
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

/// A struct to hold the UMI counts
/// encodes the UMI as a u32
/// and the count as a u32
#[derive(Debug, Default, Serialize)]
pub struct UmiCounter {
    map: Mutex<HashMap<u32, u32>>,
}

impl Clone for UmiCounter {
    fn clone(&self) -> Self {
        let map = self.map.lock().unwrap().clone();
        Self {
            map: Mutex::new(map),
        }
    }
}

impl UmiCounter {
    pub fn new() -> Self {
        Self {
            map: Mutex::new(HashMap::new()),
        }
    }

    pub fn umi2u32(umi: &Vec<u8>) -> u32 {
        if umi.len() > 16 {
            panic!("UMI length is greater than 16")
        }
        let mut res = 0;
        for &b in umi {
            res <<= 2;
            match b {
                b'A' => res |= 0,
                b'C' => res |= 1,
                b'G' => res |= 2,
                b'T' => res |= 3,
                _ => panic!("Invalid UMI base"),
            }
        }
        res
    }

    pub fn add(&self, umi: &Vec<u8>) {
        let umi_e = Self::umi2u32(umi);
        let mut map = self.map.lock().unwrap();
        *map.entry(umi_e).or_insert(0) += 1;
    }
}

/// Holds the barcode and UMI counts
/// encodes the barcode as a u32
/// and the UMI as a u32
#[derive(Debug, Default, Serialize)]
pub struct BarcodeUmiCounter {
    map: Mutex<HashMap<u32, UmiCounter>>,
}

impl Clone for BarcodeUmiCounter {
    fn clone(&self) -> Self {
        let map = self.map.lock().unwrap().clone();
        Self {
            map: Mutex::new(map)
        }
    }
}


impl BarcodeUmiCounter {
    pub fn new() -> Self {
        Self {
            map: Mutex::new(HashMap::new()),
        }
    }

    pub fn barcodes2u32(
        b1_idx: usize,
        b2_idx: usize,
        b3_idx: usize,
        b4_idx: usize,) -> u32 {
        
        let b1 = b1_idx as u8;
        let b2 = b2_idx as u8;
        let b3 = b3_idx as u8;
        let b4 = b4_idx as u8;

        ((b1 as u32) << 24) | ((b2 as u32) << 16) | ((b3 as u32) << 8) | (b4 as u32)
    }

    pub fn add(&self, 
        b1_idx: usize,
        b2_idx: usize,
        b3_idx: usize,
        b4_idx: usize,
        umi: &Vec<u8>) {
        let barcode_e = Self::barcodes2u32(b1_idx, b2_idx, b3_idx, b4_idx);
        let mut map = self.map.lock().unwrap();
        map.entry(barcode_e).or_insert_with(UmiCounter::new).add(umi);
    }

    pub fn write_barcode_stats(&self, filename: &str) -> std::io::Result<()> {
        let mut writer = File::create(filename).map(BufWriter::new)?;
        writer.write(b"barcode,total_umi,unique_umi,median_umi\n")?;
        for (barcode, umi_counter) in self.map.lock().unwrap().iter() {
            let umi_counts: Vec<u32> = umi_counter.map.lock().unwrap().values().cloned().collect();
            let total_umis = umi_counts.iter().sum::<u32>();
            let unique_umis = umi_counts.len() as u32;

            let mut sorted_counts = umi_counts;
            sorted_counts.sort_unstable();
            let median_umi = sorted_counts[sorted_counts.len() / 2];

            writeln!(writer, "{},{},{},{}", barcode, total_umis, unique_umis, median_umi)?;
        }
        Ok(())
    }

}

/// Holds the base composition
/// counts for each base
#[derive(Debug, Default, Serialize, Clone)]
pub struct BaseComposition {
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    pub n: u64,
}

impl BaseComposition {
    pub fn empty(&self) -> bool {
        self.a == 0 && self.c == 0 && self.g == 0 && self.t == 0 && self.n == 0
    }
}

impl BaseComposition {
    pub fn add_base(&mut self, base: u8) {
        match base {
            b'A' => self.a += 1,
            b'C' => self.c += 1,
            b'G' => self.g += 1,
            b'T' => self.t += 1,
            b'N' => self.n += 1,
            _ => panic!("Invalid base"),
        }
    }
}

/// Holds the UMI base composition
/// for each position in the UMI
/// 
#[derive(Debug, Default, Serialize, Clone)]
pub struct UMIBaseComposition {
    pub bases: Vec<BaseComposition>,
}

impl UMIBaseComposition {
    pub fn new(umi_len: usize) -> Self {
        let mut bases = Vec::with_capacity(umi_len);
        for _ in 0..umi_len {
            bases.push(BaseComposition::default());
        }
        Self { bases }
    }

    pub fn add(&mut self, umi: &Vec<u8>) {
        for (i, &base) in umi.iter().enumerate() {
            self.bases[i].add_base(base);
        }
    }

    pub fn write_umi_base_composition(&self, filename: &str) -> std::io::Result<()> {
        let mut writer = File::create(filename).map(BufWriter::new)?;
        writer.write(b"position,a,c,g,t,n\n")?;

        for (i, base) in self.bases.iter().enumerate() {
            if ! base.empty() {
                writeln!( writer, "{},{},{},{},{},{}", i, base.a, base.c, base.g, base.t, base.n)?;
            }
        }

        Ok(())
    }
    
}