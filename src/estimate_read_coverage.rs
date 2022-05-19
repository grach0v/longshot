//! This module contains a function for estimating the mean read coverage from a BAM file.
//!
//! This is used in Longshot if the user specifies the ```-A``` option, which calculates the mean read
//! coverage and then estimates a maximum read coverage based on this mean coverage.

// extern crates
extern crate rust_htslib;

// use declarations
use errors::*;
// use rust_htslib::bam;
use rust_htslib::bam::Read;
use util::{get_interval_lst, print_time, GenomicInterval, OpenedBamFiles};

/// Calculates the mean coverage of the input BAM file over the input region
///
/// #Arguments
/// -```bam_file```: the input BAM file name
/// -```interval```: the (optional) GenomicInterval within which variants should be called
///                  the reads that are used for estimating the alignment parameters are also
///                  limited to this region.
///
/// #Returns
/// Returns a result containing the mean coverage.
///
/// #Errors
/// - ```BamOpenError```: error opening the BAM file (without using index)
/// - ```IndexedBamOpenError```: error opening the indexed BAM file
/// - ```IndexedBamFetchError```: error fetching region from the indexed BAM file
/// - ```IndexedBamOpenError```: error opening the indexed BAM file
/// - ```IndexedBamPileupReadError```: error reading a pileup from the indexed BAM
/// - ```BamHeaderTargetLenAccessError```: error accessing a target (contig) len from BAM header

pub fn calculate_mean_coverage(
    bam_files_iteraction: &mut OpenedBamFiles,
    intervals: &Option<Vec<GenomicInterval>>,
) -> Result<f64> {

    let interval_lst: Vec<GenomicInterval> = get_interval_lst(bam_files_iteraction, intervals)?;

    let mut numerator: usize = 0;
    let mut denumerator: usize = 0;

    let iterable_intervals: Vec<Option<GenomicInterval>> = match intervals {
        Some(values) => {
            values.iter().map(|x| Some(x.clone())).collect()
        }
        None => {
            vec![None]
        }
    };

    for (mut bam_ix, bam) in bam_files_iteraction.open_indexed_files.iter_mut().zip(bam_files_iteraction.open_files.iter()) {
        for interval in iterable_intervals.iter() {
            // count variables and etc
            let mut prev_tid = u32::MAX;
            let mut bam_covered_positions = 0; // total positions
            let mut total_read_bases = 0;
            let mut total_bam_ref_positions = 0;

            // interval_lst has either the single specified genomic region, or list of regions covering all chromosomes
            // for more information about this design decision, see get_interval_lst implementation in util.rs

            for iv in interval_lst.iter() {
                bam_ix
                    .fetch((iv.tid, iv.start_pos, iv.end_pos + 1))
                    .chain_err(|| ErrorKind::IndexedBamFetchError)?;

                // iterate over the BAM pileups for every position in this interval
                for p in bam_ix.pileup() {
                    let pileup = p.chain_err(|| ErrorKind::IndexedBamPileupReadError)?;

                    let tid: u32 = pileup.tid();

                    // if we're on a different contig/chrom than the previous BAM record, we need to read
                    // in the sequence for that contig/chrom from the FASTA into the ref_seq vector
                    if tid != prev_tid {
                        total_bam_ref_positions += bam
                            .header()
                            .target_len(tid)
                            .chain_err(|| ErrorKind::BamHeaderTargetLenAccessError)?;
                    }

                    // if we've over-stepped the bounds of the region then skip to next interval
                    if tid != iv.tid || pileup.pos() < iv.start_pos || pileup.pos() > iv.end_pos {
                        prev_tid = tid;
                        continue;
                    }

                    let mut depth: usize = 0;

                    // pileup the bases for a single position and count number of each base
                    for alignment in pileup.alignments() {
                        let record = alignment.record();

                        // filter out some read QC failures
                        if record.is_unmapped()
                            || record.is_secondary()
                            || record.is_quality_check_failed()
                            || record.is_duplicate()
                            || record.is_supplementary()
                        {
                            continue;
                        }

                        // increment depth count
                        depth += 1;
                    }

                    bam_covered_positions += 1;
                    total_read_bases += depth;
                    prev_tid = tid;
                }
            }

            // here we possibly print a warning
            // this warning is to catch the case where the -A option is used but the BAM file alignments
            // are limited to a small region and that region wasn't specified with --region argument.
            // This causes the estimated mean coverage to be way off.
            let total_ref_positions: usize = match interval {
                Some(ref iv) => (iv.end_pos - iv.start_pos + 1) as usize,
                None => {
                    // output a warning if the number of covered bases is significantly less than the ref positions
                    if total_bam_ref_positions / 2 > bam_covered_positions {
                        eprintln!("{} WARNING: Max coverage calculation is highly likely to be incorrect. The number of reference \
                                    bases covered by the bam file ({}) differs significantly from the expected number of positions in the \
                                    reference ({}). If you are using a bam file that only covers part of the genome, please specify \
                                    this region exactly with the --region argument so the number of reference bases is known. \
                                    Alternatively, disable maximum coverage filtering by setting -C to a large number.",
                                print_time(), bam_covered_positions, total_bam_ref_positions);
                    }

                    total_bam_ref_positions as usize
                }
            };

            // print the total reference positions and number of observed bases in BAM file
            eprintln!(
                "{} Total reference positions: {}",
                print_time(),
                total_ref_positions
            );
            eprintln!("{} Total bases in bam: {}", print_time(), total_read_bases);

            numerator += total_read_bases;
            denumerator += total_ref_positions;
        }
    }

    Ok(numerator as f64 / denumerator as f64)
}

