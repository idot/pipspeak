mod barcodes;
mod cli;
mod config;
mod log;

use anyhow::Result;
use chrono::Local;
use clap::Parser;
use cli::Cli;
use config::Config;
use fxread::{initialize_reader, FastxRead, Record};
use gzp::{
    deflate::Gzip,
    par::compress::{ParCompress, ParCompressBuilder},
};
use indexmap::map::MutableEntryKey;
use indicatif::ProgressBar;
use log::{FileIO, Log, Parameters, Statistics, Timing};
use std::{
    fs::File,
    io::Write,
    time::{Duration, Instant},
};
use psutil::process::Process;
use std::rc::Rc;
use std::cell::RefCell;


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

fn match_records<'a>(rec1: &'a Record, rec2: &'a Record, offset: usize, config: &Config,  statistics: &mut Statistics) -> Option<(&'a Record, &'a Record, usize, Vec<usize> )> {
    let mut pos = 0;
    let mut barcode_indices = Vec::new();
    
    for i in 0..config.barcode_count() {
        if let Some((new_pos, bc_idx)) = config.match_subsequence(&rec1.seq(), i, pos, if i == 0 { Some(offset) } else { None }) {
            pos = new_pos;
            barcode_indices.push(bc_idx);
        } else {
            statistics.num_filtered[i] += 1;
            return None;
        }
    }
    
    statistics.passing_reads += 1;
    Some((rec1, rec2, pos, barcode_indices))
}



fn match_umi<'a>(rec1: &'a Record, rec2: &'a Record, pos: usize, barcode_indices: Vec<usize>, umi_len: usize, statistics: &mut Statistics) -> Option<(&'a Record, &'a Record, usize, Vec<usize>, Vec<u8>)> {
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
            Some((rec1, rec2, pos + umi_len, barcode_indices, umi))
        }
    }
}

fn construct_match<'a>(rec1: &'a Record, rec2: &'a Record, pos: usize, barcode_indices: &[usize], umi: &Vec<u8>, config: &Config, statistics: &mut Statistics) -> (Vec<u8>, Vec<u8>, &'a Record, &'a Record) {
    let mut construct_seq = config.build_barcode(barcode_indices);
    for (i, &idx) in barcode_indices.iter().enumerate() {
        statistics.counter_maps.add(idx, i);
    }
    statistics.barcode_umi_counter.add(barcode_indices, umi);
    statistics.umi_base_composition.add(umi);
    construct_seq.extend_from_slice(umi);
    
    let construct_qual = rec1.qual().unwrap()[pos - construct_seq.len()..pos].to_vec();
    (construct_seq, construct_qual, rec1, rec2)
}


fn parse_records(
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

    let record_iter = r1
        .zip(r2)
        .inspect(|_| statistics.total_reads += 1)
        .enumerate()
        .map(|(idx, pair)| {
            if idx % 1000000 == 0 || (idx < 1000 && idx % 100 == 0) {
                let process = Process::current().unwrap();
                let mem_info = process.memory_info().unwrap();
                let used_mem = mem_info.rss();
                let used_gb = used_mem as f64 / 1024.0 / 1024.0 / 1024.0;
                let msg = format!("Processed {} reads, used memory: {:.2}Gb\n", idx, used_gb);
                print!("{}", msg);
                pb.set_message(msg);
            }
            pair
        })
        .filter_map(|(rec1, rec2)| {
            match_records(&rec1, &rec2, offset, config, &mut statistics)
        })
        .filter_map(|(rec1, rec2, pos, barcode_indices)| {
            match_umi(rec1, rec2, pos, barcode_indices, umi_len, &mut statistics)
        })
        .map(|(rec1, rec2, pos, barcode_indices, umi)| {
            construct_match(rec1, rec2, pos, &barcode_indices, &umi, config, &mut statistics)
        });
    

    for (c_seq, c_qual, rec1, rec2) in record_iter {
        statistics.whitelist.insert(c_seq.clone());
        write_to_fastq(r1_out, rec1.id(), &c_seq, &c_qual)?;
        write_to_fastq(r2_out, rec2.id(), rec2.seq(), rec2.qual().unwrap())?;
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

/// Sets the number of threads to use for writing R1 and R2 files
fn set_threads(num_threads: usize) -> (usize, usize) {
    if num_threads == 0 {
        set_threads(num_cpus::get())
    } else if num_threads == 1 {
        (1, 1)
    } else {
        if num_threads % 2 == 0 {
            (num_threads / 2, num_threads / 2)
        } else {
            (num_threads / 2, num_threads / 2 + 1)
        }
    }
}

fn main() -> Result<()> {
    let args = Cli::parse();
    let config = Config::from_file(&args.config, args.exact, args.linkers)?;
    let r1 = initialize_reader(&args.r1)?;
    let r2 = initialize_reader(&args.r2)?;

    let r1_filename = args.prefix.clone() + "_R1.fq.gz";
    let r2_filename = args.prefix.clone() + "_R2.fq.gz";
    let log_filename = args.prefix.clone() + "_log.yaml";
    let whitelist_filename = args.prefix.clone() + "_whitelist.txt";
    let countermaps_filename = args.prefix.clone() + "_barcode_position_counts.tsv";
    let barcodes_umi_filename = args.prefix.clone() + "_barcode_umi_stats.tsv";
    let umi_stats_filename = args.prefix.clone() + "_umi_composition_stats.tsv";

    let (r1_threads, r2_threads) = set_threads(args.threads);
    let mut r1_writer: ParCompress<Gzip> = ParCompressBuilder::new()
        .num_threads(r1_threads)?
        .from_writer(File::create(&r1_filename)?);
    let mut r2_writer: ParCompress<Gzip> = ParCompressBuilder::new()
        .num_threads(r2_threads)?
        .from_writer(File::create(&r2_filename)?);

    let timestamp = Local::now().to_string();
    let start_time = Instant::now();

    let umi_len = if config.umi_len() == 0 { args.umi_len }else{ config.umi_len() };

    let statistics = parse_records(
        r1,
        r2,
        &mut r1_writer,
        &mut r2_writer,
        &config,
        args.offset,
        umi_len,
    )?;
    statistics.whitelist_to_file(&whitelist_filename)?;
    statistics.counter_maps_to_file(&countermaps_filename, &config)?;
    statistics.barcode_umi_stats_to_file(&barcodes_umi_filename)?;
    statistics.umi_base_composition.write_umi_base_composition(&umi_stats_filename)?;

    let elapsed_time = start_time.elapsed().as_secs_f64();
    let timing = Timing {
        timestamp,
        elapsed_time,
    };

    let parameters = Parameters {
        offset: args.offset,
        umi_len: args.umi_len,
        exact_matching: args.exact,
        write_linkers: args.linkers,
        pipspeak_version: env!("CARGO_PKG_VERSION").to_string(),
    };

    let file_io = FileIO {
        readpath_r1: args.r1,
        readpath_r2: args.r2,
        writepath_r1: r1_filename,
        writepath_r2: r2_filename,
        whitelist_path: whitelist_filename,
    };

    let log = Log {
        parameters,
        timing,
        statistics,
        file_io,
    };

    if !args.quiet {
        log.stderr()?;
    }
    log.to_file(&log_filename)?;

    Ok(())
}
