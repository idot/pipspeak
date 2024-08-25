use crate::barcodes::{Barcodes, Spacer};
use anyhow::Result;
use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct ConfigYamlRead {
    barcodes: std::collections::HashMap<String, String>,
    spacers: std::collections::HashMap<String, String>,
    parameters: Option<ConfigParameters>,
}


#[derive(Debug, Deserialize)]
pub struct ConfigYaml {
    barcodes: Vec<String>,
    spacers: Vec<String>,
    parameters: Option<ConfigParameters>,
}


#[derive(Debug, Deserialize)]
pub struct ConfigParameters {
    umi_len: usize,
}

pub struct Config {
    barcodes: Vec<Barcodes>,
    linkers: bool,
    umi_len: usize,
}
impl Config {
    pub fn from_file(path: &str, exact: bool, linkers: bool) -> Result<Self> {
        let contents = std::fs::read_to_string(path)?;
        let read_yaml = serde_yaml::from_str::<ConfigYamlRead>(&contents)?;
        let yaml = ConfigYaml {
            barcodes: read_yaml.barcodes.values().cloned().collect(),
            spacers: read_yaml.spacers.values().cloned().collect(),
            parameters: read_yaml.parameters,
        };
        Self::from_yaml(yaml, exact, linkers)
    }

    pub fn from_yaml(yaml: ConfigYaml, exact: bool, linkers: bool) -> Result<Self> {
        let mut barcodes = Vec::new();
        for (idx, (barcode_path, spacer)) in yaml.barcodes.iter().zip(yaml.spacers.iter()).enumerate() {
            let barcode = if idx < yaml.spacers.len() {
                let spacer = Spacer::from_str(spacer);
                Self::load_barcode(barcode_path, Some(&spacer), exact)?
            } else {
                Self::load_barcode(barcode_path, None, exact)?
            };
            barcodes.push(barcode);
        }

        let umi_len = yaml.parameters.map(|p| p.umi_len).unwrap_or(0);

        Ok(Self {
            barcodes,
            linkers,
            umi_len,
        })
    }

    pub fn build_barcode(&self, indices: &[usize]) -> Vec<u8> {
        let mut bc = Vec::new();
        for (idx, &barcode_idx) in indices.iter().enumerate() {
            bc.extend_from_slice(
                self.barcodes[idx]
                    .get_barcode(barcode_idx, self.linkers)
                    .expect(&format!("Invalid barcode index in bc{}", idx + 1)),
            );
        }
        bc
    }

    pub fn barcode_count(&self) -> usize {
        self.barcodes.len()
    }

   
    fn load_barcode(path: &str, spacer: Option<&Spacer>, exact: bool) -> Result<Barcodes> {
        if let Some(spacer) = spacer {
            Barcodes::from_file_with_spacer(path, spacer, exact)
        } else {
            Barcodes::from_file(path, exact)
        }
    }

    /// Matches a subsequence starting from `pos` against one of the barcode sets.
    /// Returns the end nucleotide position of the match and the within-set barcode index
    pub fn match_subsequence(
        &self,
        seq: &[u8],
        set_idx: usize,
        pos: usize,
        offset: Option<usize>,
    ) -> Option<(usize, usize)> {
        let bc = match self.barcodes.get(set_idx){
            Some(bc) => bc,
            None => panic!("Invalid set index: {}", set_idx),
        };
        if let Some(off) = offset {
            bc.match_subsequence(seq, pos, pos + bc.len() + off)
        } else {
            bc.match_subsequence(seq, pos, pos + bc.len())
        }
    }


    /// Returns the length of the UMI
    pub fn umi_len(&self) -> usize {
        self.umi_len
    }

    /// Returns the barcode based on index
    pub fn get_barcode(&self, b_index: usize, position: usize) -> Option<&[u8]> {
        self.barcodes.get(b_index).and_then(|bc| bc.get_barcode(position, self.linkers))
    }

}

#[cfg(test)]
mod testing {

    use super::*;

    const TEST_PATH: &str = "data/config_v3.yaml";

    #[test]
    fn load_yaml() {
        let config = Config::from_file(TEST_PATH, false, false);
        assert!(config.is_ok());
    }

    #[test]
    fn load_yaml_umi_len() {
        let config = Config::from_file("data/config_v3_umi_len.yaml", false, false);
        assert!(config.is_ok());
        assert!(config.unwrap().umi_len == 8)
    }

    #[test]
    fn load_yaml_exact() {
        let config = Config::from_file(TEST_PATH, true, false);
        assert!(config.is_ok());
    }
    /*
    #[test]
    fn barcode_lengths() {
        let config = Config::from_file(TEST_PATH, false, false).unwrap();
        assert_eq!(config.bc1.len(), 8 + 3);
        assert_eq!(config.bc2.len(), 6 + 3);
        assert_eq!(config.bc3.len(), 6 + 5);
        assert_eq!(config.bc4.len(), 8);
    }

    #[test]
    fn barcode_lengths_exact() {
        let config = Config::from_file(TEST_PATH, true, false).unwrap();
        assert_eq!(config.bc1.len(), 8 + 3);
        assert_eq!(config.bc2.len(), 6 + 3);
        assert_eq!(config.bc3.len(), 6 + 5);
        assert_eq!(config.bc4.len(), 8);
    }

    #[test]
    fn barcode_sequences() {
        let config = Config::from_file(TEST_PATH, false, false).unwrap();

        assert_eq!(config.bc1.get_barcode(0, true).unwrap(), b"AGAAACCAATG");
        assert_eq!(config.bc1.get_barcode(95, true).unwrap(), b"TCTTTGACATG");
        assert_eq!(config.bc1.get_barcode(96, true), None);

        assert_eq!(config.bc1.get_barcode(0, false).unwrap(), b"AGAAACCA");
        assert_eq!(config.bc1.get_barcode(95, false).unwrap(), b"TCTTTGAC");
        assert_eq!(config.bc1.get_barcode(96, false), None);

        assert_eq!(config.bc2.get_barcode(0, true).unwrap(), b"TCTGTGGAG");
        assert_eq!(config.bc2.get_barcode(95, true).unwrap(), b"GTAATCGAG");
        assert_eq!(config.bc2.get_barcode(96, true), None);

        assert_eq!(config.bc2.get_barcode(0, false).unwrap(), b"TCTGTG");
        assert_eq!(config.bc2.get_barcode(95, false).unwrap(), b"GTAATC");
        assert_eq!(config.bc2.get_barcode(96, false), None);

        assert_eq!(config.bc3.get_barcode(0, true).unwrap(), b"AAAGTGTCGAG");
        assert_eq!(config.bc3.get_barcode(95, true).unwrap(), b"CTGAAGTCGAG");
        assert_eq!(config.bc3.get_barcode(96, false), None);

        assert_eq!(config.bc3.get_barcode(0, false).unwrap(), b"AAAGTG");
        assert_eq!(config.bc3.get_barcode(95, false).unwrap(), b"CTGAAG");
        assert_eq!(config.bc3.get_barcode(96, false), None);

        assert_eq!(config.bc4.get_barcode(0, true).unwrap(), b"CTGGGTAT");
        assert_eq!(config.bc4.get_barcode(95, true).unwrap(), b"AAACTACA");
        assert_eq!(config.bc4.get_barcode(96, true), None);

        assert_eq!(config.bc4.get_barcode(0, false).unwrap(), b"CTGGGTAT");
        assert_eq!(config.bc4.get_barcode(95, false).unwrap(), b"AAACTACA");
        assert_eq!(config.bc4.get_barcode(96, false), None);
    }

    #[test]
    fn barcode_sequences_exact() {
        let config = Config::from_file(TEST_PATH, true, false).unwrap();

        assert_eq!(config.bc1.get_barcode(0, true).unwrap(), b"AGAAACCAATG");
        assert_eq!(config.bc1.get_barcode(95, true).unwrap(), b"TCTTTGACATG");
        assert_eq!(config.bc1.get_barcode(96, true), None);

        assert_eq!(config.bc2.get_barcode(0, true).unwrap(), b"TCTGTGGAG");
        assert_eq!(config.bc2.get_barcode(95, true).unwrap(), b"GTAATCGAG");
        assert_eq!(config.bc2.get_barcode(96, true), None);

        assert_eq!(config.bc3.get_barcode(0, true).unwrap(), b"AAAGTGTCGAG");
        assert_eq!(config.bc3.get_barcode(95, true).unwrap(), b"CTGAAGTCGAG");
        assert_eq!(config.bc3.get_barcode(96, true), None);

        assert_eq!(config.bc4.get_barcode(0, true).unwrap(), b"CTGGGTAT");
        assert_eq!(config.bc4.get_barcode(95, true).unwrap(), b"AAACTACA");
        assert_eq!(config.bc4.get_barcode(96, true), None);
    }

    #[test]
    fn construct_building_a() {
        let config = Config::from_file(TEST_PATH, false, false).unwrap();
        let bc = config.build_barcode(0, 0, 0, 0);
        let exp = [
            "AGAAACCA".as_bytes(),
            "TCTGTG".as_bytes(),
            "AAAGTG".as_bytes(),
            "CTGGGTAT".as_bytes(),
        ]
        .concat();
        assert_eq!(bc, exp);
    }

    #[test]
    fn construct_building_b() {
        let config = Config::from_file(TEST_PATH, false, false).unwrap();
        let bc = config.build_barcode(0, 95, 0, 95);
        let exp = [
            "AGAAACCA".as_bytes(),
            "GTAATC".as_bytes(),
            "AAAGTG".as_bytes(),
            "AAACTACA".as_bytes(),
        ]
        .concat();
        assert_eq!(bc, exp);
    }

    #[test]
    fn construct_building_a_exact() {
        let config = Config::from_file(TEST_PATH, true, false).unwrap();
        let bc = config.build_barcode(0, 0, 0, 0);
        let exp = [
            "AGAAACCA".as_bytes(),
            "TCTGTG".as_bytes(),
            "AAAGTG".as_bytes(),
            "CTGGGTAT".as_bytes(),
        ]
        .concat();
        assert_eq!(bc, exp);
    }

    #[test]
    fn construct_building_b_exact() {
        let config = Config::from_file(TEST_PATH, true, false).unwrap();
        let bc = config.build_barcode(0, 95, 0, 95);
        let exp = [
            "AGAAACCA".as_bytes(),
            "GTAATC".as_bytes(),
            "AAAGTG".as_bytes(),
            "AAACTACA".as_bytes(),
        ]
        .concat();
        assert_eq!(bc, exp);
    }
  */

}
