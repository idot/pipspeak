mod barcodes;
mod cli;
mod config;
mod log;
mod parser;

use anyhow::Result;
use chrono::Local;
use clap::Parser;
use cli::Cli;
use config::Config;
use fxread::initialize_reader;
use gzp::{
    deflate::Gzip,
    par::compress::{ParCompress, ParCompressBuilder},
};


use ::log::{LevelFilter, set_max_level};
use log::{FileIO, Log, Parameters, Timing};
use std::{
    fs::File,
    time::Instant,
};


use crate::parser::parse_records;


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

    env_logger::init();
    let log_level = match args.loglevel.to_lowercase().as_str() {
        "error" => LevelFilter::Error,
        "warn" => LevelFilter::Warn,
        "info" => LevelFilter::Info,
        "debug" => LevelFilter::Debug,
        "trace" => LevelFilter::Trace,
        _ => LevelFilter::Info, // Default to Info if the log level is not recognized
    };
    set_max_level(log_level);

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
