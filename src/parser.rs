use std::{
    default, io::Write, time::Duration
};
use anyhow::Result;
use psutil::process::Process;
use fxread::{FastxRead, Record};
use indicatif::ProgressBar;
use gzp::{
    deflate::Gzip,
    par::compress::{ParCompress},
};
use serde::de;

use crate::log::Statistics;
use crate::config::Config;

fn match_records(rec1: &Record, offset: usize, config: &Config, statistics: &mut Statistics) -> Option<(usize, Vec<usize>)> {
    let mut pos = 0;
    let mut barcode_indices = Vec::new();
    let default_offset = Some(2);

    for i in 0..config.barcode_count() {
        if let Some((new_pos, bc_idx)) = config.match_subsequence(rec1.seq(), i, pos, if i == 0 { Some(offset) } else { default_offset }) {
            pos = pos+new_pos;
            barcode_indices.push(bc_idx);
        } else {
            statistics.num_filtered[i] += 1;
            return None;
        }
    }
    
    statistics.passing_reads += 1;
    Some((pos, barcode_indices))
}

fn match_umi(rec1: &Record, pos: usize, umi_len: usize, statistics: &mut Statistics) -> Option<(usize, Vec<u8>)> {
    if rec1.seq().len() < pos + umi_len {
        statistics.num_filtered_umi += 1;
        None
    } else {
        let umi = rec1.seq()[pos..pos + umi_len].to_vec();
        let contains_n = umi.iter().any(|&base| base == b'N');
        if contains_n {
            statistics.num_filtered_umi += 1;
            None
        } else {
            Some((pos + umi_len, umi))
        }
    }
}

fn construct_match(rec1: &Record, pos: usize, barcode_indices: &[usize], umi: &Vec<u8>, config: &Config, statistics: &mut Statistics) -> (Vec<u8>, Vec<u8>) {
    let mut construct_seq = config.build_barcode(barcode_indices);
    for (i, &idx) in barcode_indices.iter().enumerate() {
        statistics.counter_maps.add(idx, i);
    }
    statistics.barcode_umi_counter.add(barcode_indices, umi);
    statistics.umi_base_composition.add(umi);
    construct_seq.extend_from_slice(umi);
    
    let construct_qual = rec1.qual().unwrap()[pos - construct_seq.len()..pos].to_vec();
    (construct_seq, construct_qual)
}

fn processed_message(idx: usize) -> String {
    let process = Process::current().unwrap();
    let mem_info = process.memory_info().unwrap();
    let used_mem = mem_info.rss();
    let used_gb = used_mem as f64 / 1024.0 / 1024.0 / 1024.0;
    let msg = format!("Processed {} reads, used memory: {:.2}Gb\n", idx, used_gb);
    msg
}

pub fn parse_records(
    r1: Box<dyn FastxRead<Item = Record>>,
    r2: Box<dyn FastxRead<Item = Record>>,
    r1_out: &mut ParCompress<Gzip>,
    r2_out: &mut ParCompress<Gzip>,
    config: &Config,
    offset: usize,
    umi_len: usize,
) -> Result<Statistics> {
    let pb = ProgressBar::new_spinner();
    pb.enable_steady_tick(Duration::from_millis(100));
    let mut statistics = Statistics::new(config.barcode_count());

    let record_iter = r1.zip(r2).enumerate();

    for (idx, (rec1, rec2)) in record_iter {
        statistics.total_reads += 1;

        if idx % 1000000 == 0 || (idx < 1000 && idx % 100 == 0) {
            let msg = processed_message(idx);
            print!("{}", msg);
            pb.set_message(msg);
        }

        if let Some((pos, barcode_indices)) = match_records(&rec1, offset, config, &mut statistics) {
            if let Some((pos, umi)) = match_umi(&rec1, pos, umi_len, &mut statistics) {
                let (c_seq, c_qual) = construct_match(&rec1, pos, &barcode_indices, &umi, config, &mut statistics);
                
                statistics.whitelist.insert(c_seq.clone());
                write_to_fastq(r1_out, rec1.id(), &c_seq, &c_qual)?;
                write_to_fastq(r2_out, rec2.id(), rec2.seq(), rec2.qual().unwrap())?;
            }
        }
    }

    statistics.calculate_metrics();
    pb.finish_with_message(format!(
        "Processed {} reads, {} passed filters ({:.4}%)",
        statistics.total_reads,
        statistics.passing_reads,
        statistics.fraction_passing * 100.0
    ));
    Ok(statistics)


}


/// Writes a record to a gzip fastq file
fn write_to_fastq<W: Write>(writer: &mut W, id: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}


#[cfg(test)]
mod testing {

    use super::*;

    const TEST_PATH: &str = "data/config_v3.yaml";

    #[test]
    fn parse_v3() {
        let config = Config::from_file(TEST_PATH, false, false).unwrap();
        let mut statistics = Statistics::new(config.barcode_count());
        let seq = b"NATACTGAATATGGTAATCGAGATCTGATCGAGGAAAGACAGTACACTTCGAGTGTGATATCTGTCTCTCTC".to_vec();
        let qual = b"1".repeat(72).to_vec();
        let fastq = fxread::Record::new_fastq_from_parts(b"id", &seq, &qual).unwrap();
        let result_record = match_records(&fastq, 5, &config, &mut statistics);
        assert_eq!(result_record, Some((41, vec![41, 95, 70, 18])));
        assert_eq!(statistics.passing_reads, 1);
        let result_umi = match_umi(&fastq, 41, 12, &mut statistics);
        assert_eq!(result_umi, Some((53, b"GTACACTTCGAG".to_vec())));
        assert_eq!(statistics.num_filtered_umi, 0);
        let result_seq = b"TACTGAATGTAATCATCTGAGAAAGACAGTACACTTCGAG".to_vec();
        let (seq, qual) = construct_match(&fastq, 55, &result_record.unwrap().1, &result_umi.unwrap().1, &config, &mut statistics);
        assert_eq!(seq, result_seq);
        assert_eq!(qual, b"1".repeat(40).to_vec())
    }
}