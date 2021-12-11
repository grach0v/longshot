//! Miscellaneous utility functions for printing, parsing, type conversion, etc.

use bio::stats::{LogProb, Prob};
use chrono::prelude::*;
use clap::ArgMatches;
use errors::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use hashbrown::HashMap;

pub static INDEX_FREQ: usize = 1000;
pub static MAX_VCF_QUAL: f64 = 500.0;

pub fn print_time() -> String {
    Local::now().format("%Y-%m-%d %H:%M:%S").to_string()
}

// use this spacer instead of calling print_time() to have the spaces match up with
// lines that document the time
pub static SPACER: &str = "                   ";

pub fn parse_u8(argmatch: &ArgMatches, arg_name: &str) -> Result<u8> {
    let parse_result: u8 = argmatch
        .value_of(arg_name)
        .chain_err(|| format!("{} not defined.", arg_name))?
        .parse::<u8>()
        .chain_err(|| format!("{} must be an integer between 0 and 255!", arg_name))?;
    Ok(parse_result)
}

pub fn parse_u32(argmatch: &ArgMatches, arg_name: &str) -> Result<u32> {
    let parse_result: u32 = argmatch
        .value_of(arg_name)
        .chain_err(|| format!("{} not defined.", arg_name))?
        .parse::<u32>()
        .chain_err(|| format!("{} must be a positive integer!", arg_name))?;
    Ok(parse_result)
}

pub fn parse_usize(argmatch: &ArgMatches, arg_name: &str) -> Result<usize> {
    let parse_result: usize = argmatch
        .value_of(arg_name)
        .chain_err(|| format!("{} not defined.", arg_name))?
        .parse::<usize>()
        .chain_err(|| format!("{} must be a positive integer!", arg_name))?;
    Ok(parse_result)
}

pub fn parse_positive_f64(argmatch: &ArgMatches, arg_name: &str) -> Result<f64> {
    let parse_result: f64 = argmatch
        .value_of(arg_name)
        .chain_err(|| format!("{} not defined.", arg_name))?
        .parse::<f64>()
        .chain_err(|| format!("{} must be a positive float!", arg_name))?;
    ensure!(
        parse_result > 0.0,
        format!("{} must be a positive float!", arg_name)
    );
    Ok(parse_result)
}

pub fn parse_nonnegative_f64(argmatch: &ArgMatches, arg_name: &str) -> Result<f64> {
    let parse_result: f64 = argmatch
        .value_of(arg_name)
        .chain_err(|| format!("{} not defined.", arg_name))?
        .parse::<f64>()
        .chain_err(|| format!("{} must be a nonnegative float!", arg_name))?;
    ensure!(
        parse_result >= 0.0,
        format!("{} must be a positive float!", arg_name)
    );
    Ok(parse_result)
}

pub fn parse_prob_into_logprob(argmatch: &ArgMatches, arg_name: &str) -> Result<LogProb> {
    let parse_result: f64 = argmatch
        .value_of(arg_name)
        .chain_err(|| format!("{} not defined.", arg_name))?
        .parse::<f64>()
        .chain_err(|| format!("{} must be a float between 0.0 and 1.0!", arg_name))?;
    ensure!(
        parse_result >= 0.0 && parse_result <= 1.0,
        format!("{} must be a float between 0.0 and 1.0!", arg_name)
    );
    Ok(LogProb::from(Prob(parse_result)))
}

pub fn parse_flag(argmatch: &ArgMatches, arg_name: &str) -> Result<bool> {
    let is_flag_set = match argmatch.occurrences_of(arg_name) {
        0 => false,
        1 => true,
        _ => {
            bail!(format!("{} cannot be specified multiple times.", arg_name));
        }
    };
    Ok(is_flag_set)
}

// this is really ugly. TODO a less verbose implementation
pub fn parse_region_string(
    region_string: &str,
    chrom_lengths: &HashMap<String, u32>
) -> Result<GenomicInterval> {

    match region_string {
        r if r.contains(":") && r.contains("-") => {
            let split1: Vec<&str> = r.split(":").collect();
            if split1.len() != 2 {
                bail!("Invalid format for region. Please use <chrom> or <chrom:start-stop>");
            }
            let split2: Vec<&str> = split1[1].split("-").collect();
            if split2.len() != 2 {
                bail!("Invalid format for region. Please use <chrom> or <chrom:start-stop>");
            }
            let iv_chrom = split1[0].to_string();
            let iv_start = split2[0]
                .parse::<u32>()
                .chain_err(|| "Invalid position value specified in region string.")?; // read in as 1-based inclusive range
            let iv_end = split2[1]
                .parse::<u32>()
                .chain_err(|| "Invalid position value specified in region string.")?; // read in as 1-based inclusive range


            let tlen = chrom_lengths[&iv_chrom];

            ensure!(
                iv_start > 0,
                "--region start position must be greater than 0."
            );

            ensure!(
                iv_start < tlen,
                "--region start position exceeds the length of the contig."
            );
            ensure!(
                iv_end <= tlen,
                "--region end position exceeds the length of the contig."
            );

            Ok(GenomicInterval {
                chrom: iv_chrom,
                start_pos: iv_start - 1, // convert to 0-based inclusive range
                end_pos: iv_end - 1,     // convert to 0-based inclusive range
            })
        }
        r => {
            let r_str = r.to_string();
            let mut tid: u32 = 0;

            let tlen = chrom_lengths[&r_str];

            Ok(GenomicInterval {
                chrom: r_str,
                start_pos: 0,
                end_pos: tlen - 1,
            })
        }
    }
}

#[derive(Clone)]

// Enum ?
// Save max number and then correct
pub struct GenomicInterval {
    pub chrom: String,
    // chromosome name
    pub start_pos: u32,
    // start of interval (0-indexed)
    pub end_pos: u32,
    // end of interval (0-indexed, inclusive)
}

#[derive(Clone)]
pub struct DensityParameters {
    pub n: usize,
    pub len: usize,
    pub gq: f64,
}

pub struct BamFileInteraction {
    pub file_names: Vec<String>,
    pub open_files: Vec<rust_htslib::bam::Reader>,
    pub chrom_to_tid: HashMap<String, Vec<i32>>,
    pub chrom_to_len: HashMap<String, u32>,
    pub file_index_tid_to_chrom: Vec<Vec<String>>,
    pub target_names: Vec<String>
}

impl BamFileInteraction {
    pub fn new(names: Vec<String>) -> Result<BamFileInteraction> {
        // Create open_files
        let bam_files: Vec<rust_htslib::bam::Reader> = names
            .iter()
            .map(|name| bam::Reader::from_path(name)
                .chain_err(|| ErrorKind::BamOpenError))
            .collect::<_>().ok();

        // Create chrom_to_tid
        // And chrom_to_len
        // and file_index_tid_to_chrom

        let mut chrom_to_tid: HashMap<String, Vec<i32>>;
        let mut chrom_to_len: HashMap<String, u32>;
        let mut file_index_tid_to_chrom: Vec<Vec<String>>;

        for (bam, bam_index) in bam_files.into_iter().zip(0 .. bam_files.len()) {
            let cur_names = parse_target_names_opened(&bam)?;
            file_index_tid_to_chrom.push(vec![]);

            for chrom in chrom_to_tid.keys() {
                let Some(vec) = chrom_to_tid.get(chrom);
                vec.push(-1);
            }

            for (chrom, chrom_index) in cur_names.into_iter().zip(0 .. cur_names.len()) {
                file_index_tid_to_chrom[file_index_tid_to_chrom.len() - 1].push(chrom);

                let tlen = bam
                    .header()
                    .target_len(chrom_index as u32)
                    .chain_err(|| ErrorKind::BamHeaderTargetLenAccessError)?;

                if let (Some(vec), Some(&len)) = (chrom_to_tid.get(&chrom), chrom_to_len.get(&chrom)) {
                    vec[vec.len() - 1] = chrom_index as i32;

                    ensure! (
                        tlen == len,
                        "--different chrom lengths of one chromosome."
                    );

                } else {
                    let mut new_vec = vec![-1; bam_index];
                    new_vec.push(chrom_index as i32);
                    chrom_to_tid.insert(chrom, new_vec);

                    chrom_to_len.insert(chrom, tlen);
                }
            }
        }

        Ok(
            BamFileInteraction {
                file_names: names,
                open_files: bam_files,
                chrom_to_tid: chrom_to_tid,
                chrom_to_len: chrom_to_len,
                file_index_tid_to_chrom: file_index_tid_to_chrom,
                target_names: chrom_to_len.keys().cloned().collect()
            }
        )
    }
}

pub fn u8_to_string(u: &[u8]) -> Result<String> {
    Ok(String::from_utf8(u.to_vec()).chain_err(|| "Error converting u8 to String.")?)
}

//
pub fn dna_vec(u: &[u8]) -> Vec<char> {
    let mut v: Vec<char> = Vec::with_capacity(u.len());
    for cu in u.to_ascii_uppercase() {
        let c = cu as char;
        //assert!(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
        if c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N' {
            v.push(c);
        } else {
            eprintln!(
                "Warning: Unexpected base \"{}\" encountered. Replaced with \"N\".",
                c
            );
            v.push('N');
        }
    }
    v
}

pub fn has_non_acgt(s: &String) -> bool {
    for c in s.chars() {
        if !(c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            return true;
        }
    }
    false
}


pub fn parse_target_names(bam_file: &String) -> Result<Vec<String>> {
    // TODO replace with std::str_from_utf8 

    let bam = bam::Reader::from_path(bam_file).chain_err(|| ErrorKind::BamOpenError)?;
    let header_view = bam.header();
    let target_names_dec: Vec<&[u8]> = header_view.target_names();
    let mut target_names: Vec<String> = vec![];

    for t_name_dec in target_names_dec {
        let mut name_vec: Vec<char> = vec![];
        for decr in t_name_dec {
            let dec: u8 = *decr;
            name_vec.push(dec as char);
        }
        let name_string: String = name_vec.into_iter().collect();
        target_names.push(name_string);
    }

    Ok(target_names)
}

pub fn parse_target_names_opened(bam: &rust_htslib::bam::Reader) -> Result<Vec<String>> {
    // TODO replace with std::str_from_utf8 

    // let bam = bam::Reader::from_path(bam_file).chain_err(|| ErrorKind::BamOpenError)?;
    let header_view = bam.header();
    let target_names_dec: Vec<&[u8]> = header_view.target_names();
    let mut target_names: Vec<String> = vec![];

    for t_name_dec in target_names_dec {
        let mut name_vec: Vec<char> = vec![];
        for decr in t_name_dec {
            let dec: u8 = *decr;
            name_vec.push(dec as char);
        }
        let name_string: String = name_vec.into_iter().collect();
        target_names.push(name_string);
    }

    Ok(target_names)
}


pub fn get_whole_genome_intervals(bam_files_iteraction: &BamFileInteraction) -> Result<Vec<GenomicInterval>> {

    let mut intervals: Vec<GenomicInterval> = vec![];

    // for (tid, t_name_dec) in target_names_dec.iter().enumerate() {
    //     let mut name_vec: Vec<char> = vec![];
    //     for decr in t_name_dec.iter() {
    //         let dec: u8 = *decr;
    //         name_vec.push(dec as char);
    //     }
    //     let name_string: String = name_vec.into_iter().collect();
    //     intervals.push(GenomicInterval {
    //         // tid: tid as u32,
    //         chrom: name_string,
    //         start_pos: 0,
    //         end_pos: header_view
    //             .target_len(tid as u32)
    //             .chain_err(|| format!("Error accessing target len for tid {}", tid))?
    //             - 1,
    //     });
    // }

    for chrom in bam_files_iteraction.target_names {
        intervals.push(GenomicInterval {
            chrom: chrom,
            start_pos: 0,
            end_pos: bam_files_iteraction.chrom_to_len[&chrom] - 1
        });
    }

    Ok(intervals)
}

// the strategy used for iterating over the BAM entries is as follows:
// if a genomic region was specified, then ```get_interval_lst``` puts that genomic region into a vector (interval_lst)
// containing only that one region. Then we iterate over the region list and seek every region,
// which effectively means we just seek that single genomic interval.
//
// if a genomic region was not specified, then we want to iterate over the whole BAM file.
// so get_interval_lst returns a list of genomic intervals that contain each whole chromosome
// as described by the BAM header SQ lines. Then we iterate over the interval lst and seek each
// interval separately, which effectively just iterates over all the BAM entries.
//
// the reason for this strange design (instead of either fetching a region beforehand or not and
// then just iterating over all of it is the following:
// if an indexed reader is used, and fetch is never called, pileup() hangs and accessing BAM records
// doesn't work.

//remove result 
pub fn get_interval_lst(
    bam_files_iteraction: &BamFileInteraction,
    interval: &Option<Vec<GenomicInterval>>,
) -> Result<Vec<GenomicInterval>> {
    match interval {
        &Some(ref iv) => Ok(iv.clone()),
        &None => get_whole_genome_intervals(bam_files_iteraction),
    }
}
