use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::record::CigarStringView;
use rust_htslib::bam::record::Cigar;
use std::error::Error;
use util::*;
use bio::stats::{LogProb, Prob, PHREDProb};
use bio::io::fasta;
use bio::pattern_matching::bndm;
use realignment;

static VERBOSE: bool = false;

// an extension of the rust-htslib cigar representation
// has the cigar operation as well as the position on the reference of the operation start,
// and the position on the read of the operation start
pub struct CigarPos {
    pub cig: Cigar,
    pub ref_pos: u32,
    pub read_pos: u32,
}

//************************************************************************************************
// BEGINNING OF RUST-HTSLIB BASED CODE *****************************************************************
//************************************************************************************************

// This function (create_augmented_cigarlist) is a modification of Rust-Htslib function rust_htslib::bam::record::CigarStringView::read_pos
// https://github.com/rust-bio/rust-htslib/blob/master/src/bam/record.rs
// Rust-Htslib license is copied here as per its terms:
//The MIT License (MIT)
//
//Copyright (c) 2016 Johannes Köster, the Rust-Htslib team.

//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

quick_error! {
    #[derive(Debug)]
    pub enum CigarOrAnchorError {
        UnsupportedOperation(msg: String) {
            description("Unsupported CIGAR operation")
            display(x) -> ("{}: {}", x.description(), msg)
        }
        UnexpectedOperation(msg: String) {
            description("CIGAR operation not allowed at this point")
            display(x) -> ("{}: {}", x.description(), msg)
        }
        AnchorRangeOutsideRead(msg: String) {
            description("Attempted to find sequence anchors for a range completely outside of the sequence.")
            display(x) -> ("{}: {}", x.description(), msg)
        }
    }
}

/// Creates a vector of CigarPos, which contain each cigar operation and the reference position
/// and read position at the start of the operation
///
/// # Arguments
///
/// * `cigar_string_view` - cigar string to create cigar position list for
///
pub fn create_augmented_cigarlist(refpos: u32,
                                  cigar_string_view: &CigarStringView)
                                  -> Result<Vec<CigarPos>, CigarOrAnchorError> {
    let mut rpos = refpos;
    let mut qpos = 0u32; // position within read
    let mut j = 0; // index into cigar operation vector

    let mut cigar_list: Vec<CigarPos> = Vec::with_capacity(cigar_string_view.len());
    // find first cigar operation referring to qpos = 0 (and thus bases in record.seq()),
    // because all augmentations of qpos and rpos before that are invalid
    for (i, c) in cigar_string_view.iter().enumerate() {
        match c {
            &Cigar::Match(_) |
            &Cigar::Diff(_) |
            &Cigar::Equal(_) |
            // this is unexpected, but bwa + GATK indel realignment can produce insertions
            // before matching positions
            &Cigar::Ins(_) => {
                j = i;
                break;
            }
            &Cigar::SoftClip(_) => {
                j = i;
                break;
            }
            &Cigar::Del(_) => {
                return Err(CigarOrAnchorError::UnexpectedOperation(
                    "'deletion' (D) found before any operation describing read sequence".to_owned()
                ));
            }
            &Cigar::Back(_) => {
                return Err(CigarOrAnchorError::UnsupportedOperation(
                    "'back' (B) operation is deprecated according to htslib/bam_plcmd.c and is not in SAMv1 spec".to_owned()
                ));
            }
            &Cigar::RefSkip(_) => {
                return Err(CigarOrAnchorError::UnexpectedOperation(
                    "'reference skip' (N) found before any operation describing read sequence".to_owned()
                ));
            }
            &Cigar::HardClip(_) if i > 0 && i < cigar_string_view.len() - 1 => {
                return Err(CigarOrAnchorError::UnexpectedOperation(
                    "'hard clip' (H) found in between operations, contradicting SAMv1 spec that hard clips can only be at the ends of reads".to_owned()
                ));
            }
            // if we have reached the end of the CigarString with only pads and hard clips, we have no read position matching the variant
            &Cigar::Pad(_) | &Cigar::HardClip(_) if i == cigar_string_view.len() - 1 => return Ok(cigar_list),
            // skip leading HardClips and Pads, as they consume neither read sequence nor reference sequence
            &Cigar::Pad(_) | &Cigar::HardClip(_) => ()
        }
    }

    while j < cigar_string_view.len() {
        // append CigarPos to list (contains cigar op, ref position, read position)
        match &cigar_string_view[j] {
            &Cigar::Match(_) |
            &Cigar::Diff(_) |
            &Cigar::Equal(_) |
            &Cigar::Ins(_) |
            &Cigar::Del(_) |
            &Cigar::RefSkip(_) => {
                // push the current cigar operation, ref position, and query position onto the list
                cigar_list.push(CigarPos {
                    cig: cigar_string_view[j].clone(),
                    ref_pos: rpos,
                    read_pos: qpos,
                });
            }
            &Cigar::SoftClip(_) |
            &Cigar::Pad(_) |
            &Cigar::HardClip(_) |
            &Cigar::Back(_) => {}
        }

        // move the reference and query positions forward
        match &cigar_string_view[j] {
            &Cigar::Match(l) | &Cigar::Diff(l) | &Cigar::Equal(l) => {
                rpos += l;
                qpos += l;
                j += 1;
            }
            &Cigar::SoftClip(l) => {
                qpos += l;
                j += 1;
            }
            &Cigar::Ins(l) => {
                qpos += l;
                j += 1;
            }
            &Cigar::RefSkip(l) |
            &Cigar::Del(l) => {
                rpos += l;
                j += 1;
            }
            &Cigar::Pad(_) => {
                j += 1;
            }
            &Cigar::HardClip(_) if j < cigar_string_view.len() - 1 => {
                return Err(CigarOrAnchorError::UnexpectedOperation(
                    "'hard clip' (H) found in between operations, contradicting SAMv1 spec that hard clips can only be at the ends of reads".to_owned()
                ));
            }
            &Cigar::HardClip(_) => {
                j += 1;
            }
            &Cigar::Back(_) => {
                return Err(CigarOrAnchorError::UnsupportedOperation(
                    "'back' (B) operation is deprecated according to htslib/bam_plcmd.c and is not in SAMv1 spec".to_owned()
                ));
            }
        }
    }

    Ok(cigar_list)
}

pub struct AnchorPositions {
    pub left_anchor_ref: u32,
    pub right_anchor_ref: u32,
    pub left_anchor_read: u32,
    pub right_anchor_read: u32,
}

pub fn find_anchors(bam_record: &Record,
                    cigarpos_list: &Vec<CigarPos>,
                    var_interval: GenomicInterval,
                    ref_seq: &Vec<char>,
                    read_seq: &Vec<char>,
                    target_names: &Vec<String>,
                    extract_params: ExtractFragmentParameters)
                    -> Result<Option<AnchorPositions>, CigarOrAnchorError> {

    //let min_window_length = extract_params.min_window_length;
    let max_window_padding = extract_params.max_window_padding;

    let anchor_length = extract_params.anchor_length as u32;

    let l_max = var_interval.start_pos as usize - max_window_padding;
    let r_max = var_interval.end_pos as usize + max_window_padding;
    let mut ref_seq_max_window: Vec<u8> = vec![];

    for c in ref_seq[l_max..r_max+1].iter() {
        ref_seq_max_window.push(*c as u8);
    }

    if var_interval.chrom != target_names[bam_record.tid() as usize] ||
        (var_interval.start_pos as i32) >=
            bam_record.cigar().end_pos().expect("Error while accessing CIGAR end position") ||
        (var_interval.end_pos as i32) < bam_record.pos() {
        return Err(CigarOrAnchorError::AnchorRangeOutsideRead(
            "Attempted to find sequence anchors for a range completely outside of the sequence.".to_owned()
        ));
    }

    let mut left_anchor_ref: u32 = 0;
    let mut right_anchor_ref: u32 = 0;
    let mut left_anchor_read: u32 = 0;
    let mut right_anchor_read: u32 = 0;

    let mut left_ix = 0;
    let mut right_ix = 0;

    let contains_pos = |pos: u32, cigar_op_start: u32, cigar_op_length: u32| {
        cigar_op_start <= pos && cigar_op_start + cigar_op_length > pos
    };

    // find left_ix and right_ix.
    // these are the indexes into the CIGAR of the CIGAR operation holding the
    // left side of the range and right side of the range

    // how slow is this? could be sped up with binary search in the future.
    for (i, cigarpos) in cigarpos_list.iter().enumerate() {
        match cigarpos.cig {
            Cigar::Match(l) |
            Cigar::Diff(l) |
            Cigar::Equal(l) |
            Cigar::Ins(l) |
            Cigar::Del(l) |
            Cigar::RefSkip(l) => {
                if contains_pos(var_interval.start_pos, cigarpos.ref_pos, l) {
                    left_ix = i;
                }
                if contains_pos(var_interval.end_pos, cigarpos.ref_pos, l) {
                    right_ix = i;
                    break;
                }
            }
            Cigar::Pad(_) |
            Cigar::Back(_) |
            Cigar::SoftClip(_) |
            Cigar::HardClip(_) => {
                return Err(CigarOrAnchorError::UnexpectedOperation(
                    "CIGAR operation found in cigarpos_list that should have been removed already.".to_owned()
                ));
            }
        }
    }

    // step backwards to find left anchor
    let mut match_len_left = 0;
    let mut seen_indel_left = false;
    let mut found_anchor_left = false;
    //println!("**************************************************");
    //println!("searching for left anchor...");

    for i in (0..left_ix + 1).rev() {
        match cigarpos_list[i].cig {
            Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) => {
                match_len_left += l;
                let mut potential_anchor = false;
                if seen_indel_left && match_len_left > anchor_length {
                    // we have found a sufficiently long match and it's NOT the cigar op containing var_interval.start_pos
                    left_anchor_ref = cigarpos_list[i].ref_pos + match_len_left - anchor_length;
                    left_anchor_read = cigarpos_list[i].read_pos + match_len_left - anchor_length;
                    potential_anchor = true;
                } else if !seen_indel_left &&
                    cigarpos_list[i].ref_pos < var_interval.start_pos - anchor_length {
                    // we have found a sufficiently long match but it's the cigar ops surrounding var_interval.start_pos
                    // the var_interval.start_pos we are trying to anchor is just inside one huge match
                    left_anchor_ref = var_interval.start_pos - anchor_length;
                    left_anchor_read = cigarpos_list[i].read_pos +
                        (var_interval.start_pos - anchor_length -
                            cigarpos_list[i].ref_pos);
                    potential_anchor = true;
                }

                if potential_anchor {
                    while left_anchor_ref >= cigarpos_list[i].ref_pos {

                        let l_anc: usize = left_anchor_ref as usize;
                        let r_anc: usize = (left_anchor_ref + anchor_length) as usize;
                        if l_anc <= 0 || r_anc >= ref_seq.len() {
                            break;
                        }

                        // check if there is an exact match between anchor sequence on read and ref
                        let anchor_on_read: Vec<char> = read_seq[(left_anchor_read as usize)..
                            (left_anchor_read + anchor_length) as usize]
                            .to_vec();
                        assert_eq!(anchor_on_read.len(),anchor_length as usize);
                        let anchor_on_ref: Vec<char> = ref_seq[l_anc..r_anc].to_vec();
                        let anchor_match = anchor_on_read == anchor_on_ref;

                        // check that the anchor sequence is unique in the region
                        let mut pattern: Vec<u8>= vec![];
                        for c in anchor_on_ref.iter() {
                            pattern.push(*c as u8);
                        };
                        assert_eq!(pattern.len(), anchor_length as usize);
                        let bndm = bndm::BNDM::new(&pattern);
                        let occ: Vec<usize> = bndm.find_all(&ref_seq_max_window).collect();

                        //println!("ref_seq_max_window: {}", String::from_utf8(ref_seq_max_window.clone()).unwrap());
                        //println!("anchor_on_ref: {}", String::from_utf8(pattern.clone()).unwrap());
                        //let s: String = anchor_on_read.into_iter().collect();
                        //println!("anchor_on_read: {}", s);
                        //println!("pattern: {}", String::from_utf8(pattern.clone()).unwrap());
                        //println!("anchor_match: {} occ.len(): {} l_anc <= l_max: {}", anchor_match, occ.len(), l_anc <= l_max);
                        //println!("***************");

                        if (anchor_match && occ.len() == 1) || l_anc <= l_max {
                            found_anchor_left = true;
                            break;
                        }

                        left_anchor_ref -= 1;
                        left_anchor_read -= 1;
                    }
                    if found_anchor_left {
                        break;
                    }
                }
            }
            Cigar::Ins(_) |
            Cigar::Del(_) |
            Cigar::RefSkip(_) => {
                match_len_left = 0;
                seen_indel_left = true;
            }
            Cigar::Pad(_) |
            Cigar::Back(_) |
            Cigar::SoftClip(_) |
            Cigar::HardClip(_) => {
                return Err(CigarOrAnchorError::UnexpectedOperation(
                    "CIGAR operation found in cigarpos_list that should have been removed already.".to_owned()
                ));
            }
        }
    }

    if !found_anchor_left {
        return Ok(None); // failed to find a left anchor
    }

    //println!("**************************************************");
    //println!("searching for right anchor...");
    // step forwards to find right anchor
    let mut match_len_right = 0;
    let mut seen_indel_right = false;
    let mut found_anchor_right = false;
    for i in right_ix..cigarpos_list.len() {
        match cigarpos_list[i].cig {
            Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) => {
                match_len_right += l;
                let mut potential_anchor = false;
                if seen_indel_right && match_len_right > anchor_length {
                    // we have found a sufficiently long match and it's NOT the cigar op containing var_interval.end_pos
                    right_anchor_ref = cigarpos_list[i].ref_pos + l - match_len_right +
                        anchor_length;
                    right_anchor_read = cigarpos_list[i].read_pos + l - match_len_right +
                        anchor_length;
                    potential_anchor = true;
                } else if !seen_indel_right &&
                    cigarpos_list[i].ref_pos + l > var_interval.end_pos + anchor_length {
                    // we have found a sufficiently long match but it's the cigar ops surrounding var_interval.end_pos
                    // the var_interval.end_pos we are trying to anchor is just inside one huge match
                    right_anchor_ref = var_interval.end_pos + anchor_length;
                    right_anchor_read = cigarpos_list[i].read_pos +
                        (var_interval.end_pos + anchor_length -
                            cigarpos_list[i].ref_pos);
                    potential_anchor = true;
                }

                if potential_anchor {
                    while right_anchor_ref < cigarpos_list[i].ref_pos + l {

                        let l_anc: usize = (right_anchor_ref - anchor_length) as usize;
                        let r_anc: usize = right_anchor_ref as usize;
                        if l_anc <= 0 || r_anc >= ref_seq.len() {
                            break;
                        }

                        // check if there is an exact match between anchor sequence on read and ref
                        let anchor_on_read: Vec<char> = read_seq[(right_anchor_read - anchor_length) as usize..
                            right_anchor_read as usize]
                            .to_vec();
                        assert_eq!(anchor_on_read.len(),anchor_length as usize);
                        let anchor_on_ref: Vec<char> = ref_seq[l_anc..r_anc].to_vec();
                        let anchor_match = anchor_on_read == anchor_on_ref;

                        // check that the anchor sequence is unique in the region
                        let mut pattern: Vec<u8> = vec![];
                        for c in anchor_on_ref.iter() {
                            pattern.push(*c as u8);
                        };
                        assert_eq!(pattern.len(), anchor_length as usize);
                        let bndm = bndm::BNDM::new(&pattern);
                        let occ: Vec<usize> = bndm.find_all(&ref_seq_max_window).collect();

                        //println!("ref_seq_max_window: {}", String::from_utf8(ref_seq_max_window.clone()).unwrap());
                        //println!("anchor_on_ref: {}", String::from_utf8(pattern.clone()).unwrap());
                        //let s: String = anchor_on_read.into_iter().collect();
                        //println!("anchor_on_read: {}", s);
                        //println!("pattern: {}", String::from_utf8(pattern.clone()).unwrap());
                        //println!("anchor_match: {} occ.len(): {} r_anc >= r_max: {}", anchor_match, occ.len(), r_anc >= r_max);
                        //println!("***************");

                        if (anchor_match && occ.len() == 1) || r_anc >= r_max {
                            found_anchor_right = true;
                            break;
                        }

                        right_anchor_ref += 1;
                        right_anchor_read += 1;
                    }
                    if found_anchor_right {
                        break;
                    }
                }
            }
            Cigar::Ins(_) |
            Cigar::Del(_) |
            Cigar::RefSkip(_) => {
                match_len_right = 0;
                seen_indel_right = true;
            }
            Cigar::Pad(_) |
            Cigar::Back(_) |
            Cigar::SoftClip(_) |
            Cigar::HardClip(_) => {
                return Err(CigarOrAnchorError::UnexpectedOperation(
                    "CIGAR operation found in cigarpos_list that should have been removed already.".to_owned()
                ));
            }
        }
    }

    if !found_anchor_right {
        return Ok(None); // failed to find a right anchor
    }

    // commented out 1/16/18 -- we no longer care if we found anchors or not,
    // we're going to use them anyway.

    // return none if read window or ref window is larger or smaller than the allowed lengths
    //let ref_window_len = (right_anchor_ref - left_anchor_ref) as usize;
    //let read_window_len = (right_anchor_read - left_anchor_read) as usize;
    //if ref_window_len < min_window_length || read_window_len < min_window_length ||
    //    ref_window_len > max_window_length || read_window_len > max_window_length {
    //    return Ok(None);
    //}

    // return none if any of the anchors are out of bounds
    if right_anchor_ref as usize >= ref_seq.len() ||
        right_anchor_read as usize >= bam_record.seq().len() {
        return Ok(None);
    }

    // return anchor pos
    Ok(Some(AnchorPositions {
        left_anchor_ref: left_anchor_ref,
        right_anchor_ref: right_anchor_ref,
        left_anchor_read: left_anchor_read,
        right_anchor_read: right_anchor_read,
    }))
}

fn extract_var_cluster(read_seq: &Vec<char>,
                       ref_seq: &Vec<char>,
                       var_cluster: Vec<Var>,
                       anchors: AnchorPositions,
                       extract_params: ExtractFragmentParameters,
                       align_params: AlignmentParameters)
                       -> Vec<FragCall> {

    let mut calls: Vec<FragCall> = vec![];

    //let ref_window = ref_seq[(anchors.left_anchor_ref as usize)..
    //(anchors.right_anchor_ref as usize) + 1]
    //        .to_vec();
    let window_capacity = (anchors.right_anchor_ref - anchors.left_anchor_ref + 10) as usize;

    let read_window: Vec<char> = read_seq[(anchors.left_anchor_read as usize)..
        (anchors.right_anchor_read as usize) + 1]
        .to_vec();

    let mut max_score: LogProb = LogProb::ln_zero();
    let mut max_hap: usize = 0;
    let n_vars: usize = var_cluster.len() as usize; // number of variants in cluster
    let n_haps: usize = 2usize.pow(n_vars as u32); // number of possible haplotypes for cluster
    assert!(n_vars <= extract_params.short_hap_max_snvs);

    let in_hap = |var: usize, hap: usize| (hap & 2usize.pow(var as u32)) > 0;

    // allele_scores0[i] contains the Log sum of the probabilities of all short haplotypes
    // that had a '0' at the ith variant of the cluster.
    // similarly for allele_scores1, for '1' variants.
    let mut allele_scores0: Vec<LogProb> = vec![LogProb::ln_zero(); n_vars];
    let mut allele_scores1: Vec<LogProb> = vec![LogProb::ln_zero(); n_vars];
    let mut score_total: LogProb = LogProb::ln_zero();

    if VERBOSE {
        for var in var_cluster.clone() {
            println!("{} {} {} {}",
                     var.chrom,
                     var.pos0,
                     var.ref_allele,
                     var.var_allele);
        }
        let read_seq_str: String = read_window.clone().into_iter().collect();
        println!("read: {}", read_seq_str);
    }


    for hap in 0..n_haps {
        let mut hap_window: Vec<char> = Vec::with_capacity(window_capacity);
        let mut i: usize = anchors.left_anchor_ref as usize;
        for var in 0..n_vars {
            while i < var_cluster[var].pos0 {
                hap_window.push(ref_seq[i]);
                i += 1;
            }
            if in_hap(var, hap) {
                for c in var_cluster[var].var_allele.chars() {
                    hap_window.push(c);
                }
            } else {
                for c in var_cluster[var].ref_allele.chars() {
                    hap_window.push(c);
                }
            }
            i += var_cluster[var].ref_allele.len();
        }

        while i <= anchors.right_anchor_ref as usize {
            hap_window.push(ref_seq[i]);
            i += 1;
        }

        // we now want to score hap_window
        let score: LogProb = match extract_params.alignment_type{
            AlignmentType::NumericallyStableAllAlignment => {
                realignment::sum_all_alignments_numerically_stable(&read_window,
                                                                   &hap_window,
                                                                   align_params.ln(),
                                                                   extract_params.band_width)
            }
            AlignmentType::FastAllAlignment => {
                realignment::sum_all_alignments(&read_window,
                                                &hap_window,
                                                align_params,
                                                extract_params.band_width)
            }
            AlignmentType::MaxAlignment => {
                realignment::max_alignment(&read_window,
                                           &hap_window,
                                           align_params.ln(),
                                           extract_params.band_width)
            }
        };

        for var in 0..n_vars {
            if in_hap(var, hap) {
                allele_scores1[var] = LogProb::ln_add_exp(allele_scores1[var], score);
            } else {
                allele_scores0[var] = LogProb::ln_add_exp(allele_scores0[var], score);
            }
        }
        if VERBOSE {
            let hap_seq_str: String = hap_window.into_iter().collect();
            println!("hap:{} {} PHRED: {}", hap, hap_seq_str, *PHREDProb::from(score));
        }
        // add current alignment score to the total score sum
        score_total = LogProb::ln_add_exp(score_total, score);

        if score > max_score {
            max_score = score;
            max_hap = hap;
        }
    }


    for v in 0..n_vars {
        if in_hap(v, max_hap) {
            // the best haplotype has a '1' variant allele at this position

            // quality score is the probability call is wrong, in other words the ratio of the
            // sum of haplotype scores that had a '0' here over the total sum of scores
            let qual = allele_scores0[v] - score_total;
            if VERBOSE {
                println!("adding call: {} {} {} {}; allele = 1; qual = {};", var_cluster[v].chrom, var_cluster[v].pos0, var_cluster[v].ref_allele, var_cluster[v].var_allele, *Prob::from(qual));
            }
            calls.push(FragCall {
                frag_ix: None,
                var_ix: var_cluster[v].ix,
                allele: '1',
                qual: qual,
                one_minus_qual: LogProb::ln_one_minus_exp(&qual)
            });
        } else {
            // the best haplotype has a '0' ref allele at this position

            // quality score is the probability call is wrong, in other words the ratio of the
            // sum of haplotype scores that had a '1' here over the total sum of scores
            let qual = allele_scores1[v] - score_total;
            if VERBOSE {
                println!("adding call: {} {} {} {}; allele = 0; qual = {};", var_cluster[v].chrom, var_cluster[v].pos0, var_cluster[v].ref_allele, var_cluster[v].var_allele, *Prob::from(qual));
            }
            calls.push(FragCall {
                frag_ix: None,
                var_ix: var_cluster[v].ix,
                allele: '0',
                qual: qual,
                one_minus_qual: LogProb::ln_one_minus_exp(&qual)
            });
        }
    }

    if VERBOSE {
        println!("--------------------------------------");
    }

    calls
}


pub fn extract_fragment(bam_record: &Record,
                        cigarpos_list: &Vec<CigarPos>,
                        vars: Vec<Var>,
                        ref_seq: &Vec<char>,
                        target_names: &Vec<String>,
                        extract_params: ExtractFragmentParameters,
                        align_params: AlignmentParameters)
                        -> Option<Fragment> {


    if bam_record.is_quality_check_failed() || bam_record.is_duplicate() ||
        bam_record.is_secondary() || bam_record.is_unmapped() || bam_record.mapq() < extract_params.min_mapq
        {
            return None;
        }

    // TODO assert that every single variant in vars is on the same chromosome
    let id: String = u8_to_string(bam_record.qname());

    if VERBOSE {println!("Extracting fragment for read {}...", id);}

    let mut fragment = Fragment {
        id: id,
        calls: vec![],
        p_read_hap: [LogProb::from(Prob(0.5)),LogProb::from(Prob(0.5))]
    };
    let read_seq: Vec<char> = dna_vec(&bam_record.seq().as_bytes());
    let mut cluster_lst: Vec<(AnchorPositions, Vec<Var>)> = vec![];
    let mut var_anchor_lst: Vec<(Var, AnchorPositions)> = vec![];

    // populate a list with tuples of each variant, and anchor sequences for its alignment
    for ref var in vars {
        let var_interval = GenomicInterval{
            chrom: var.chrom.clone(),
            start_pos: var.pos0 as u32,
            end_pos: var.pos0 as u32
        };
        match find_anchors(&bam_record,
                                   &cigarpos_list,
                                   var_interval,
                                   &ref_seq,
                                   &read_seq,
                                   &target_names,
                                   extract_params).expect("CIGAR or Anchor Error while finding anchor sequences.") {
            Some(anchors) => {var_anchor_lst.push((var.clone(), anchors));},
            _ => {}
        };
    }

    // now that we have anchors for each var the read covers,
    // group the variants into clusters to align together if adjacent anchors overlap
    // we generate anchors for the whole cluster by taking the first-left and last-right anchor pos of the cluster
    {
        // populate cluster_lst with variant clusters
        let mut var_cluster: Vec<Var> = vec![];
        let mut var_anchors: Vec<AnchorPositions> = vec![];
        let mut l;

        // generate clusters of SNVs that should be considered
        for (var, anc) in var_anchor_lst {
            l = var_cluster.len();
            if l == 0 ||
                (anc.left_anchor_ref < var_anchors[l - 1].right_anchor_ref
                    && l < extract_params.short_hap_max_snvs) {
                var_cluster.push(var);
                var_anchors.push(anc);
            } else {

                // sequence anchor that covers the whole cluster of variants
                let combined_anchor = AnchorPositions{
                    left_anchor_ref: var_anchors[0].left_anchor_ref,
                    right_anchor_ref: var_anchors[l - 1].right_anchor_ref,
                    left_anchor_read: var_anchors[0].left_anchor_read,
                    right_anchor_read: var_anchors[l - 1].right_anchor_read
                };

                cluster_lst.push((combined_anchor, var_cluster.clone()));

                var_cluster.clear();
                var_anchors.clear();
                var_cluster.push(var);
                var_anchors.push(anc);
            }
        }

        l = var_cluster.len();

        if l > 0 {
            // sequence anchor that covers the whole cluster of variants
            let combined_anchor = AnchorPositions{
                left_anchor_ref: var_anchors[0].left_anchor_ref,
                right_anchor_ref: var_anchors[l - 1].right_anchor_ref,
                left_anchor_read: var_anchors[0].left_anchor_read,
                right_anchor_read: var_anchors[l - 1].right_anchor_read
            };

            cluster_lst.push((combined_anchor, var_cluster.clone()));
        }
    }

    // now extract alleles for the variant cluster
    for (anchors, var_cluster) in cluster_lst {

        for call in extract_var_cluster(&read_seq,
                                        ref_seq,
                                        var_cluster,
                                        anchors,
                                        extract_params,
                                        align_params) {
            fragment.calls.push(call);
        }
    }

    if fragment.calls.len() > 0 {
        Some(fragment)
    } else {
        None
    }
}

pub fn extract_fragments(bamfile_name: &String,
                         fastafile_name: &String,
                         varlist: &VarList,
                         interval: &Option<GenomicInterval>,
                         extract_params: ExtractFragmentParameters,
                         align_params: AlignmentParameters)
                         -> Vec<Fragment> {

    let t_names = parse_target_names(&bamfile_name);

    let mut prev_tid = 4294967295; // huge value so that tid != prev_tid on first iter
    let mut fasta = fasta::IndexedReader::from_file(fastafile_name).unwrap();
    let mut ref_seq: Vec<char> = vec![];

    let mut flist: Vec<Fragment> = vec![];

    // TODO: this uses a lot of duplicate code, need to figure out a better solution.
    let mut complete = 0;

    match interval {
        &Some(ref iv) => {
            let mut bam = bam::IndexedReader::from_path(bamfile_name).unwrap();
            let iv_tid = bam.header().tid(iv.chrom.as_bytes()).unwrap();
            bam.fetch(iv_tid, iv.start_pos, iv.end_pos + 1).ok().expect("Error seeking BAM file while extracting fragments.");
            for r in bam.records() {
                let record = r.unwrap();

                let tid: usize = record.tid() as usize;
                let chrom: String = t_names[record.tid() as usize].clone();

                if tid != prev_tid {
                    let mut ref_seq_u8: Vec<u8> = vec![];
                    fasta.read_all(&chrom, &mut ref_seq_u8).expect("Failed to read fasta sequence record.");
                    ref_seq = dna_vec(&ref_seq_u8);
                }

                let start_pos = record.pos();
                let end_pos =
                    record.cigar().end_pos().expect("Error while accessing CIGAR end position") - 1;

                let bam_cig: CigarStringView = record.cigar();
                let cigarpos_list: Vec<CigarPos> =
                    create_augmented_cigarlist(start_pos as u32, &bam_cig).expect("Error creating augmented cigarlist.");


                let interval = GenomicInterval {
                    chrom: chrom,
                    start_pos: start_pos as u32,
                    end_pos: end_pos as u32,
                };

                // get the list of variants that overlap this read
                let read_vars = varlist.get_variants_range(interval);

                // print the percentage of variants processed every 10%
                if read_vars.len() > 0 && ((read_vars[0].ix as f64 / varlist.lst.len() as f64) * 10.0) as usize > complete {
                    complete = ((read_vars[0].ix as f64 / varlist.lst.len() as f64) * 10.0) as usize;
                    if complete < 10 {
                        eprintln!("{}    {}% of variants processed...", print_time(), complete * 10);
                    }
                }

                let frag = extract_fragment(&record,
                                            &cigarpos_list,
                                            read_vars,
                                            &ref_seq,
                                            &t_names,
                                            extract_params,
                                            align_params);

                match frag {
                    Some(f) => {flist.push(f);}
                    None => {}
                }

                prev_tid = tid;
            }

            eprintln!("{}    100% of variants processed.",print_time());

        }
        &None => {
            let mut bam = bam::Reader::from_path(bamfile_name).unwrap();
            for r in bam.records() {
                let record = r.unwrap();

                let tid: usize = record.tid() as usize;
                let chrom: String = t_names[record.tid() as usize].clone();

                if tid != prev_tid {
                    let mut ref_seq_u8: Vec<u8> = vec![];
                    fasta.read_all(&chrom, &mut ref_seq_u8).expect("Failed to read fasta sequence record.");
                    ref_seq = dna_vec(&ref_seq_u8);
                }

                let start_pos = record.pos();
                let end_pos =
                    record.cigar().end_pos().expect("Error while accessing CIGAR end position") - 1;

                let interval = GenomicInterval {
                    chrom: chrom,
                    start_pos: start_pos as u32,
                    end_pos: end_pos as u32,
                };

                let bam_cig: CigarStringView = record.cigar();
                let cigarpos_list: Vec<CigarPos> =
                    create_augmented_cigarlist(start_pos as u32, &bam_cig).expect("Error creating augmented cigarlist.");

                // get the list of variants that overlap this read
                let read_vars = varlist.get_variants_range(interval);

                // print the percentage of variants processed every 10%
                if read_vars.len() > 0 && ((read_vars[0].ix as f64 / varlist.lst.len() as f64) * 10.0) as usize > complete {
                    complete = ((read_vars[0].ix as f64 / varlist.lst.len() as f64) * 10.0) as usize;
                    if complete < 10 {
                        eprintln!("{}    {}% of variants processed...",print_time(), complete*10);
                    }
                }

                let frag = extract_fragment(&record,
                                            &cigarpos_list,
                                            read_vars,
                                            &ref_seq,
                                            &t_names,
                                            extract_params,
                                            align_params);

                match frag {
                    Some(f) => {flist.push(f);}
                    None => {}
                }

                prev_tid = tid;
            }

            eprintln!("{}    100% of variants processed.",print_time());

        }
    };

    // label every fragment call with its index in the fragment list.
    for i in 0..flist.len() {

        for j in 0..flist[i].calls.len() {
            flist[i].calls[j].frag_ix = Some(i);
        }
    }

    flist
}

//************************************************************************************************
// END OF RUST-HTSLIB BASED CODE *****************************************************************
//************************************************************************************************


#[cfg(test)]
mod tests {
    use super::*;

}
