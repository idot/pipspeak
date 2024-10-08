use clap::Parser;

#[derive(Parser, Debug)]
#[clap(author, version, about)]
pub struct Cli {
    /// Input file for R1
    #[clap(short = 'i', long, value_parser)]
    pub r1: String,

    /// Input file for R2
    #[clap(short = 'I', long, value_parser)]
    pub r2: String,

    /// Output file prefix (output files will be named <prefix>_R[12].fq.gz)
    #[clap(short = 'p', long, value_parser, default_value = "pipspeak")]
    pub prefix: String,

    /// Number of threads to use in gzip compression (0 = all threads)
    #[clap(short = 't', long, default_value = "1")]
    pub threads: usize,

    /// The amount of nucleotides away from the start of R1 to accept a barcode
    #[clap(short = 's', long, default_value = "5")]
    pub offset: usize,

    /// The yaml config file describing the file paths of the 4 barcodes and the spacers
    #[clap(short = 'c', long, value_parser)]
    pub config: String,

    /// The length of the UMI, overriden in the config file if present
    #[clap(short = 'u', long, default_value = "12")]
    pub umi_len: usize,

    /// Offset of the UMI from the last base of the barcode (0 = no offset)
    #[clap(long, default_value = "0")]
    pub umi_offset: usize,

    /// Use exact matching instead of one mismatch
    #[clap(short = 'x', long)]
    pub exact: bool,

    /// Include linkers in the output
    #[clap(short = 'l', long)]
    pub linkers: bool,

    /// Do not write anything to stderr
    #[clap(short = 'q', long)]
    pub quiet: bool,

    /// Log level
    #[clap(short = 'e', long, default_value = "info")]
    pub loglevel: String,
}
