#!/usr/bin/env python3
#Description: ENCODE overlap_peaks 0.0.1 modified by WKL
#Author: WKL
#Environment: python3.6
#Date: 2018/04/09
#Example: python3 chip_overlap_peaks.py --in1="R1.tagAlign.gz" --peak1="r1pr1.gappedPeak.gz" --peak2="r1pr2.gappedPeak.gz" \
#   --peakpool="r1.gappedPeak.gz" --xcor1="R1.cc.qc" --pair=False --chromsizes="male.hg19.chrom.sizes" \
#    --peaktype="gappedPeak" --peakas="gappedPeak.as" --signal1="r1.fc_signal.bw" --prefix="final" \

import os
import re
import getopt
import logging
import subprocess
import sys
sys.path.insert(1,"/zs32_2/klwang/ChIP_seq/test_data/ENCODE_code")
import common
sys.path.remove("/zs32_2/klwang/ChIP_seq/test_data/ENCODE_code")

logger = logging.getLogger(__name__)
#logger.addHandler(dxpy.DXLogHandler())
logger.propagate = False
logger.setLevel(logging.INFO)

def usage():
    print("")
    print("Final peaks by overlap between replicates or pseudoreplicates.")
    print("usage: python %s -option <argument>" %sys.argv[0])
    print("  -h/--help                   ")
    print("  --in1=<STRING>              [required] Rep1 tagAlign input file.")
    print("  --in2=<STRING>                         Rep2 tagAlign input file.")    
    print("  --peak1=<STRING>                       The 1st peak file.")
    print("  --peak2=<STRING>                       The 2nd peak file.") 
    print("  --peakpool=<STRING>                    The pooled peak file.") 
    print("  --peakpoolpr1=<STRING>                 The pr1 of pooled peak file.")
    print("  --peakpoolpr2=<STRING>                 The pr2 of pooled peak file.")      
    print("  --xcor1=<STRING>            [required] Rep1 xcor input file.")
    print("  --xcor2=<STRING>                       Rep2 xcor input file.")   
    print("  --pair=<True or False>                 Paired_end or single_end.")     
    print("  --chromsizes=<STRING>                  The chromosome sizes file.")
    print("  --peaktype=<STRING>                    The type of peak.")
    print("  --peakas=<STRING>                      The peak_as file.")
    print("  --prefix=<STRING>                      The prefix for final peak file.")
    print("  --signal1=<STRING>                     The signal file of rep1.")
    print("  --signal2=<STRING>                     The signal file of rep2")
    print("  --signalpool=<STRING>                  The signal file of pool")
    print("  --fragment=<INT>                       The fragment length.")

## deal with options
try:
    opts, args = getopt.getopt( sys.argv[1:], "h", ["help", "in1=", "in2=", "peak1=", "peak2=", "peakpool=", "peakpoolpr1=", "peakpoolpr2=", "xcor1=", "xcor2=", "pair=", "chromsizes=", "peaktype=", "peakas=", "prefix=", "signal1=", "signal2=", "signalpool=", "fragment="] )
except getopt.GetoptError:
    print("ERROR: Get option error.")
    usage()
    sys.exit(2)

for opt, val in opts:
    if opt in ( "-h", "--help" ):
        usage()
        sys.exit(1)
    else:
        if opt in ( "--in1", ):
            rep1_ta = val
        if opt in ( "--in2", ):
            rep2_ta = val
        if opt in ( "--peak1", ):
            rep1_peaks = val
        if opt in ( "--peak2", ):
            rep2_peaks = val
        if opt in ( "--peakpool", ):
            pooled_peaks = val
        if opt in ( "--peakpoolpr1", ):
            pooledpr1_peaks = val
        if opt in ( "--peakpoolpr2", ):
            pooledpr2_peaks = val
        if opt in ( "--xcor1", ):
            rep1_xcor = val
        if opt in ( "--xcor2", ):
            rep2_xcor = val
        if opt in ( "--pair", ):
            paired_end = val
        if opt in ( "--chromsizes", ):
            chrom_sizes = val
        if opt in ( "--peaktype", ):
            peak_type = val
        if opt in ( "--peakas", ):
            as_file = val
        if opt in ( "--prefix", ):
            prefix = val
        if opt in ( "--signal1", ):
            rep1_signal = val
        if opt in ( "--signal2", ):
            rep2_signal = val
        if opt in ( "--signalpool", ):
            pooled_signal = val
        if opt in ( "--fragment", ):
            fragment_length = int(val)
if not "rep2_ta" in dir():
    rep2_ta = None
if not "pooledpr1_peaks" in dir():
    pooledpr1_peaks = None
if not "pooledpr2_peaks" in dir():
    pooledpr2_peaks = None
if not "rep2_xcor" in dir():
    rep2_xcor = None
if not "prefix" in dir():
    prefix = None
if not "rep1_signal" in dir():
    rep1_signal = None
if not "rep2_signal" in dir():
    rep2_signal = None
if not "pooled_signal" in dir():
    pooled_signal = None
if not "paired_end" in dir():
    paired_end = None
if not "fragment_length" in dir():
    fragment_length = None
'''
try:
    debug = False if not debug else debug
except:
    debug = False
'''
if not paired_end or paired_end=="False":
    paired_end=False
else:
    paired_end=True

def xcor_parse(fname):
    with open(fname, 'r') as xcor_file:
        if not xcor_file:
            return None
        lines = xcor_file.read().splitlines()
        line = lines[0].rstrip('\n')
        # CC_SCORE FILE format:
        #   Filename <tab>
        #   numReads <tab>
        #   estFragLen <tab>
        #   corr_estFragLen <tab>
        #   PhantomPeak <tab>
        #   corr_phantomPeak <tab>
        #   argmin_corr <tab>
        #   min_corr <tab>
        #   phantomPeakCoef <tab>
        #   relPhantomPeakCoef <tab>
        #   QualityTag
        headers = ['Filename',
                   'numReads',
                   'estFragLen',
                   'corr_estFragLen',
                   'PhantomPeak',
                   'corr_phantomPeak',
                   'argmin_corr',
                   'min_corr',
                   'phantomPeakCoef',
                   'relPhantomPeakCoef',
                   'QualityTag']
        metrics = line.split('\t')
        headers.pop(0)
        metrics.pop(0)
        xcor_qc = dict(zip(headers, metrics))
    return xcor_qc

def xcor_only(input_tagAlign, paired_end):
    input_tagAlign_file = input_tagAlign
    input_tagAlign_filename = input_tagAlign_file
    input_tagAlign_basename = input_tagAlign_file.rstrip('.gz')
    uncompressed_TA_filename = input_tagAlign_basename
    out, err = common.run_pipe(['gzip -d %s' % (input_tagAlign_filename)])
    # =================================
    # Subsample tagAlign file
    # ================================
    NREADS = 15000000
    if paired_end:
        end_infix = 'MATE1'
    else:
        end_infix = 'SE'
    subsampled_TA_filename = \
        input_tagAlign_basename + \
        ".sample.%d.%s.tagAlign.gz" % (NREADS/1000000, end_infix)
    steps = [
        'grep -v "chrM" %s' % (uncompressed_TA_filename),
        'shuf -n %d --random-source=%s' % (NREADS, uncompressed_TA_filename)]
    if paired_end:
        steps.extend([r"""awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}'"""])
    steps.extend(['gzip -cn'])
    out, err = common.run_pipe(steps, outfile=subsampled_TA_filename)
    # Calculate Cross-correlation QC scores
    CC_scores_filename = subsampled_TA_filename + ".cc.qc"
    CC_plot_filename = subsampled_TA_filename + ".cc.plot.pdf"
    # CC_SCORE FILE format
    # Filename <tab>
    # numReads <tab>
    # estFragLen <tab>
    # corr_estFragLen <tab>
    # PhantomPeak <tab>
    # corr_phantomPeak <tab>
    # argmin_corr <tab>
    # min_corr <tab>
    # phantomPeakCoef <tab>
    # relPhantomPeakCoef <tab>
    # QualityTag
    #
    # spp_tarball = SPP_VERSION_MAP.get(spp_version)
    # assert spp_tarball, "spp version %s is not supported" % (spp_version)
    # # install spp
    # subprocess.check_output(shlex.split('R CMD INSTALL %s' % (spp_tarball)))
    # run spp
    run_spp_command = '/zs32_2/klwang/ChIP_seq/test_data/ENCODE_code/phantompeakqualtools/run_spp.R'
    out, err = common.run_pipe([
        "Rscript3.3.3 %s -c=%s -p=%d -filtchr=chrM -savp=%s -out=%s"
        % (run_spp_command, subsampled_TA_filename, cpu_count(),
           CC_plot_filename, CC_scores_filename)])
    out, err = common.run_pipe([
        r"""sed -r  's/,[^\t]+//g' %s""" % (CC_scores_filename)],
        outfile="temp")
    out, err = common.run_pipe([
        "mv temp %s" % (CC_scores_filename)])
    CC_scores_file = CC_scores_filename
    CC_plot_file = CC_plot_filename
    xcor_qc = xcor_parse(CC_scores_filename)
    # Return the outputs
    output = {
        "CC_scores_file": CC_scores_file,
        "CC_plot_file": CC_plot_file,
        "paired_end": paired_end,
        "RSC": float(xcor_qc.get('relPhantomPeakCoef')),
        "NSC": float(xcor_qc.get('phantomPeakCoef')),
        "est_frag_len": float(xcor_qc.get('estFragLen'))
    }
    return output

def pool(inputs, prefix=None):
    input_filenames = []
    for input_file in inputs:
        dxf = input_file
        input_filenames.append(dxf)
    # uses last extension - presumably they are all the same
    extension = splitext(splitext(input_filenames[-1])[0])[1]
    if prefix:
        pooled_filename = prefix + "_pooled%s.gz" % (extension)
    else:
        pooled_filename = \
            '-'.join([splitext(splitext(fn)[0])[0] for fn in input_filenames]) + "_pooled%s.gz" % (extension)
    out, err = common.run_pipe([
        'gzip -dc %s' % (' '.join(input_filenames)),
        'gzip -cn'],
        outfile=pooled_filename)
    pooled = pooled_filename
    output = {
        "pooled": pooled
    }
    return output

def internal_pseudoreplicate_overlap(rep1_peaks, rep2_peaks, pooled_peaks,
                                     rep1_ta, rep1_xcor,
                                     paired_end, chrom_sizes, as_file,
                                     peak_type, prefix, fragment_length=None):
    #
    rep1_peaks_file      = rep1_peaks
    rep2_peaks_file      = rep2_peaks
    pooled_peaks_file    = pooled_peaks
    rep1_ta_file         = rep1_ta
    rep1_xcor_file       = rep1_xcor
    chrom_sizes_file     = chrom_sizes
    as_file_file         = as_file
    # Input filenames - necessary to define each explicitly because input files
    # could have the same name, in which case subsequent
    # file would overwrite previous file
    '''
    if rep1_peaks_file.rfind("/") != -1:
        rep1_peaks_fn      = '%s/rep1-%s' % (rep1_peaks_file[0:rep1_peaks_file.rfind("/")],rep1_peaks_file.split("/")[-1])
        rep2_peaks_fn      = '%s/rep2-%s' % (rep2_peaks_file[0:rep2_peaks_file.rfind("/")],rep2_peaks_file.split("/")[-1])
        pooled_peaks_fn    = '%s/pooled-%s' % (pooled_peaks_file[0:pooled_peaks_file.rfind("/")],pooled_peaks_file.split("/")[-1])
        rep1_ta_fn         = '%s/r1ta_%s' % (rep1_ta_file[0:rep1_ta_file.rfind("/")],rep1_ta_file.split("/")[-1])
        rep1_xcor_fn       = '%s/r1xc_%s' % (rep1_xcor_file[0:rep1_xcor_file.rfind("/")],rep1_xcor_file.split("/")[-1])
        chrom_sizes_fn     = 'chrom.sizes'
        as_file_fn         = '%s.as' % (peak_type)
    else:
        rep1_peaks_fn      = 'rep1-%s' % (rep1_peaks_file)
        rep2_peaks_fn      = 'rep2-%s' % (rep2_peaks_file)
        pooled_peaks_fn    = 'pooled-%s' % (pooled_peaks_file)
        rep1_ta_fn         = 'r1ta_%s' % (rep1_ta_file)
        rep1_xcor_fn       = 'r1xc_%s' % (rep1_xcor_file)
        chrom_sizes_fn     = 'chrom.sizes'
        as_file_fn         = '%s.as' % (peak_type)
    '''
    rep1_peaks_fn      = rep1_peaks_file
    rep2_peaks_fn      = rep2_peaks_file
    pooled_peaks_fn    = pooled_peaks_file
    rep1_ta_fn         = rep1_ta_file
    rep1_xcor_fn       = rep1_xcor_file
    chrom_sizes_fn     = chrom_sizes_file
    as_file_fn         = as_file_file
    # Output filenames
    if prefix:
        basename = prefix
    else:
        # strip off the peak and compression extensions
        m = re.match(
            '(.*)(\.%s)+(\.((gz)|(Z)|(bz)|(bz2)))' % (peak_type),
            pooled_peaks)
        if m:
            basename = m.group(1)
        else:
            basename = pooled_peaks
    #
    overlapping_peaks_fn    = '%s.replicated.%s' % (basename, peak_type)
    overlapping_peaks_bb_fn = overlapping_peaks_fn + '.bb'
    rejected_peaks_fn       = '%s.rejected.%s' % (basename, peak_type)
    rejected_peaks_bb_fn    = rejected_peaks_fn + '.bb'
    #
    # Intermediate filenames
    overlap_tr_fn = 'replicated_tr.%s' % (peak_type)
    overlap_pr_fn = 'replicated_pr.%s' % (peak_type)
    #
    logger.info(subprocess.check_output('set -x; ls -l', shell=True))
    # the only difference between the peak_types is how the extra columns are
    # handled
    if peak_type == "narrowPeak":
        awk_command = r"""awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}'"""
        cut_command = 'cut -f 1-10'
        bed_type = 'bed6+4'
    elif peak_type == "gappedPeak":
        awk_command = r"""awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$18-$17; if (($31/s1 >= 0.5) || ($31/s2 >= 0.5)) {print $0}}'"""
        cut_command = 'cut -f 1-15'
        bed_type = 'bed12+3'
    elif peak_type == "broadPeak":
        awk_command = r"""awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}'"""
        cut_command = 'cut -f 1-9'
        bed_type = 'bed6+3'
    else:
        assert peak_type in ['narrowPeak', 'gappedPeak', 'broadPeak'], "%s is unrecognized.  peak_type should be narrowPeak, gappedPeak or broadPeak." % (peak_type)
    #
    # Find pooled peaks that overlap Rep1 and Rep2 where overlap is defined as
    # the fractional overlap wrt any one of the overlapping peak pairs  > 0.5
    out, err = common.run_pipe([
        'intersectBed -wo -a %s -b %s' % (pooled_peaks_fn, rep1_peaks_fn),
        awk_command,
        cut_command,
        'sort -u',
        'intersectBed -wo -a stdin -b %s' % (rep2_peaks_fn),
        awk_command,
        cut_command,
        'sort -u'
        ], overlap_tr_fn)
    print(
        "%d peaks overlap with both true replicates"
        % (common.count_lines(overlap_tr_fn)))
    # this is a simplicate analysis
    # overlapping peaks are just based on pseudoreps of the one pool
    out, err = common.run_pipe([
        'cat %s' % (overlap_tr_fn),
        'sort -u'
        ], overlapping_peaks_fn)
    print(
        "%d peaks overlap"
        % (common.count_lines(overlapping_peaks_fn)))
    #
    # rejected peaks
    out, err = common.run_pipe([
        'intersectBed -wa -v -a %s -b %s' % (pooled_peaks_fn, overlapping_peaks_fn)
        ], rejected_peaks_fn)
    print("%d peaks were rejected" % (common.count_lines(rejected_peaks_fn)))
    #
    # calculate FRiP (Fraction of Reads in Peaks)
    #
    # Extract the fragment length estimate from column 3 of the
    # cross-correlation scores file or use the user-defined
    # fragment_length if given.
    if fragment_length is not None:
        fraglen = fragment_length
        fragment_length_given_by_user = True
    else:
        fraglen = common.xcor_fraglen(rep1_xcor_fn)
        fragment_length_given_by_user = False
    #
    # FRiP
    reads_in_peaks_fn = 'reads_in_%s.ta' % (peak_type)
    n_reads, n_reads_in_peaks, frip_score = common.frip(
        rep1_ta_fn, rep1_xcor_fn, overlapping_peaks_fn,
        chrom_sizes_fn, fraglen,
        reads_in_peaks_fn=reads_in_peaks_fn)
    #
    # count peaks
    npeaks_in = common.count_lines(common.uncompress(pooled_peaks_fn))
    npeaks_out = common.count_lines(overlapping_peaks_fn)
    npeaks_rejected = common.count_lines(rejected_peaks_fn)
    #
    # make bigBed files for visualization
    overlapping_peaks_bb_fn = common.bed2bb(
        overlapping_peaks_fn, chrom_sizes_fn, as_file_fn, bed_type=bed_type)
    rejected_peaks_bb_fn = common.bed2bb(
        rejected_peaks_fn, chrom_sizes_fn, as_file_fn, bed_type=bed_type)
    #
    # Upload file outputs from the local file system.
    #
    overlapping_peaks = common.compress(overlapping_peaks_fn)
    overlapping_peaks_bb = overlapping_peaks_bb_fn
    rejected_peaks = common.compress(rejected_peaks_fn)
    rejected_peaks_bb = rejected_peaks_bb_fn
    #
    output = {
        "overlapping_peaks"     : overlapping_peaks,
        "overlapping_peaks_bb"  : overlapping_peaks_bb,
        "rejected_peaks"        : rejected_peaks,
        "rejected_peaks_bb"     : rejected_peaks_bb,
        "npeaks_in"             : npeaks_in,
        "npeaks_out"            : npeaks_out,
        "npeaks_rejected"       : npeaks_rejected,
        "frip_nreads"           : n_reads,
        "frip_nreads_in_peaks"  : n_reads_in_peaks,
        "frip_score"            : frip_score,
        "fragment_length_used"  : fraglen,
        "fragment_length_given_by_user": fragment_length_given_by_user
    }
    #
    return output


def replicated_overlap(rep1_peaks, rep2_peaks, pooled_peaks,
                       pooledpr1_peaks, pooledpr2_peaks,
                       rep1_ta, rep1_xcor, rep2_ta, rep2_xcor,
                       paired_end, chrom_sizes, as_file, peak_type, prefix,
                       fragment_length=None):
    #
    rep1_peaks_file      = rep1_peaks
    rep2_peaks_file      = rep2_peaks
    pooled_peaks_file    = pooled_peaks
    pooledpr1_peaks_file = pooledpr1_peaks
    pooledpr2_peaks_file = pooledpr2_peaks
    rep1_ta_file         = rep1_ta
    rep2_ta_file         = rep2_ta
    rep1_xcor_file       = rep1_xcor
    rep2_xcor_file       = rep2_xcor
    chrom_sizes_file     = chrom_sizes
    as_file_file         = as_file
    #
    # Input filenames - necessary to define each explicitly because input files
    # could have the same name, in which case subsequent
    # file would overwrite previous file
    '''
    rep1_peaks_fn      = 'rep1-%s' % (rep1_peaks_file)
    rep2_peaks_fn      = 'rep2-%s' % (rep2_peaks_file)
    pooled_peaks_fn    = 'pooled-%s' % (pooled_peaks_file)
    pooledpr1_peaks_fn = 'pooledpr1-%s' % (pooledpr1_peaks_file)
    pooledpr2_peaks_fn = 'pooledpr2-%s' % (pooledpr2_peaks_file)
    rep1_ta_fn         = 'r1ta_%s' % (rep1_ta_file)
    rep2_ta_fn         = 'r2ta_%s' % (rep2_ta_file)
    rep1_xcor_fn       = 'r1cc_%s' % (rep1_xcor_file)
    rep2_xcor_fn       = 'r2cc_%s' % (rep2_xcor_file)
    chrom_sizes_fn     = 'chrom.sizes'
    as_file_fn         = '%s.as' % (peak_type)
    '''
    rep1_peaks_fn      = rep1_peaks_file
    rep2_peaks_fn      = rep2_peaks_file
    pooled_peaks_fn    = pooled_peaks_file
    pooledpr1_peaks_fn = pooledpr1_peaks_file
    pooledpr2_peaks_fn = pooledpr2_peaks_file
    rep1_ta_fn         = rep1_ta_file
    rep2_ta_fn         = rep2_ta_file
    rep1_xcor_fn       = rep1_xcor_file
    rep2_xcor_fn       = rep2_xcor_file
    chrom_sizes_fn     = chrom.sizes
    as_file_fn         = peak_type
    #
    # Output filenames
    if prefix:
        basename = prefix
    else:
        # strip off the peak and compression extensions
        m = re.match(
            '(.*)(\.%s)+(\.((gz)|(Z)|(bz)|(bz2)))' % (peak_type),
            pooled_peaks)
        if m:
            basename = m.group(1)
        else:
            basename = pooled_peaks
    #
    overlapping_peaks_fn    = '%s.replicated.%s' % (basename, peak_type)
    overlapping_peaks_bb_fn = overlapping_peaks_fn + '.bb'
    rejected_peaks_fn       = '%s.rejected.%s' % (basename, peak_type)
    rejected_peaks_bb_fn    = rejected_peaks_fn + '.bb'
    #
    # Intermediate filenames
    overlap_tr_fn = 'replicated_tr.%s' % (peak_type)
    overlap_pr_fn = 'replicated_pr.%s' % (peak_type)
    #
    pool_replicates_subjob = pool(inputs=[rep1_ta, rep2_ta], prefix='pooled_reps')
    # If fragment length was given by user we skip pooled_replicates
    # _xcor_subjob, set the pool_xcor_filename to None, and update
    # the flag fragment_length_given_by_user. Otherwise, run the subjob
    # to be able to extract the fragment length fron cross-correlations.
    if fragment_length is not None:
        pool_xcor_filename = None
        fraglen = fragment_length
        fragment_length_given_by_user = True
    else:
        pooled_replicates_xcor_subjob = xcor_only(pool_replicates_subjob.get("pooled"), paired_end)
        pooled_replicates_xcor_subjob.wait_on_done()
        pool_xcor_link = pooled_replicates_xcor_subjob.describe()['output'].get("CC_scores_file")
        pool_xcor_file = pool_xcor_link
        #pool_xcor_filename = 'poolcc_%s' % (pool_xcor_file)
        pool_xcor_filename = pool_xcor_file
        fraglen = common.xcor_fraglen(pool_xcor_filename)
        fragment_length_given_by_user = False
    #
    pool_replicates_subjob.wait_on_done()
    pool_ta_link = pool_replicates_subjob.describe()['output'].get("pooled")
    pool_ta_file = pool_ta_link
    #pool_ta_filename = 'poolta_%s' % (pool_ta_file)
    pool_ta_filename = pool_ta_file
    #
    logger.info(subprocess.check_output('set -x; ls -l', shell=True))
    #
    # the only difference between the peak_types is how the extra columns are
    # handled
    if peak_type == "narrowPeak":
        awk_command = r"""awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}'"""
        cut_command = 'cut -f 1-10'
        bed_type = 'bed6+4'
    elif peak_type == "gappedPeak":
        awk_command = r"""awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$18-$17; if (($31/s1 >= 0.5) || ($31/s2 >= 0.5)) {print $0}}'"""
        cut_command = 'cut -f 1-15'
        bed_type = 'bed12+3'
    elif peak_type == "broadPeak":
        awk_command = r"""awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}'"""
        cut_command = 'cut -f 1-9'
        bed_type = 'bed6+3'
    else:
        assert peak_type in ['narrowPeak', 'gappedPeak', 'broadPeak'], "%s is unrecognized.  peak_type should be narrowPeak, gappedPeak or broadPeak." % (peak_type)
    #
    # Find pooled peaks that overlap Rep1 and Rep2 where overlap is defined as
    # the fractional overlap wrt any one of the overlapping peak pairs  > 0.5
    out, err = common.run_pipe([
        'intersectBed -wo -a %s -b %s' % (pooled_peaks_fn, rep1_peaks_fn),
        awk_command,
        cut_command,
        'sort -u',
        'intersectBed -wo -a stdin -b %s' % (rep2_peaks_fn),
        awk_command,
        cut_command,
        'sort -u'
        ], overlap_tr_fn)
    print(
        "%d peaks overlap with both true replicates"
        % (common.count_lines(overlap_tr_fn)))
    # Find pooled peaks that overlap PseudoRep1 and PseudoRep2 where
    # overlap is defined as the fractional overlap wrt any one of the
    # overlapping peak pairs  > 0.5
    out, err = common.run_pipe([
        'intersectBed -wo -a %s -b %s' % (pooled_peaks_fn, pooledpr1_peaks_fn),
        awk_command,
        cut_command,
        'sort -u',
        'intersectBed -wo -a stdin -b %s' % (pooledpr2_peaks_fn),
        awk_command,
        cut_command,
        'sort -u'
        ], overlap_pr_fn)
    print(
        "%d peaks overlap with both pooled pseudoreplicates"
        % (common.count_lines(overlap_pr_fn)))
    # Combine peak lists
    out, err = common.run_pipe([
        'cat %s %s' % (overlap_tr_fn, overlap_pr_fn),
        'sort -u'
        ], overlapping_peaks_fn)
    print(
        "%d peaks overlap with true replicates or with pooled pseudoreplicates"
        % (common.count_lines(overlapping_peaks_fn)))
    # rejected peaks
    out, err = common.run_pipe([
        'intersectBed -wa -v -a %s -b %s' % (pooled_peaks_fn, overlapping_peaks_fn)
        ], rejected_peaks_fn)
    print("%d peaks were rejected" % (common.count_lines(rejected_peaks_fn)))
    # calculate FRiP (Fraction of Reads in Peaks)
    reads_in_peaks_fn = 'reads_in_%s.ta' % (peak_type)
    n_reads, n_reads_in_peaks, frip_score = common.frip(
        pool_ta_filename, pool_xcor_filename, overlapping_peaks_fn,
        chrom_sizes_fn, fraglen, reads_in_peaks_fn=reads_in_peaks_fn)
    # count peaks
    npeaks_in        = common.count_lines(common.uncompress(pooled_peaks_fn))
    npeaks_out       = common.count_lines(overlapping_peaks_fn)
    npeaks_rejected  = common.count_lines(rejected_peaks_fn)
    # make bigBed files for visualization
    overlapping_peaks_bb_fn = common.bed2bb(
        overlapping_peaks_fn, chrom_sizes_fn, as_file_fn, bed_type=bed_type)
    rejected_peaks_bb_fn = common.bed2bb(
        rejected_peaks_fn, chrom_sizes_fn, as_file_fn, bed_type=bed_type)
    # Upload file outputs from the local file system.
    overlapping_peaks       = common.compress(overlapping_peaks_fn)
    overlapping_peaks_bb    = overlapping_peaks_bb_fn
    rejected_peaks          = common.compress(rejected_peaks_fn)
    rejected_peaks_bb       = rejected_peaks_bb_fn
    #
    output = {
        "overlapping_peaks"     : overlapping_peaks,
        "overlapping_peaks_bb"  : overlapping_peaks_bb,
        "rejected_peaks"        : rejected_peaks,
        "rejected_peaks_bb"     : rejected_peaks_bb,
        "npeaks_in"             : npeaks_in,
        "npeaks_out"            : npeaks_out,
        "npeaks_rejected"       : npeaks_rejected,
        "frip_nreads"           : n_reads,
        "frip_nreads_in_peaks"  : n_reads_in_peaks,
        "frip_score"            : frip_score,
        "fragment_length_used"  : fraglen,
        "fragment_length_given_by_user": fragment_length_given_by_user
    }
    return output

def main(rep1_peaks, rep2_peaks, pooled_peaks,
         rep1_ta, rep1_xcor,
         paired_end, chrom_sizes, as_file, peak_type,
         pooledpr1_peaks=None, pooledpr2_peaks=None,
         rep2_ta=None, rep2_xcor=None,
         prefix=None,
         rep1_signal=None, rep2_signal=None, pooled_signal=None,
         fragment_length=None):

    replicate_inputs = [pooledpr1_peaks, pooledpr2_peaks]
    simplicate_experiment = \
        all([replicate_input is None for replicate_input in replicate_inputs])
    # that means pooledpr1_peaks is None and pooledpr2_peaks is None

    if simplicate_experiment:
        output = internal_pseudoreplicate_overlap(
            rep1_peaks, rep2_peaks, pooled_peaks,
            rep1_ta, rep1_xcor,
            paired_end, chrom_sizes, as_file, peak_type, prefix,
            fragment_length)
    else:
        output = replicated_overlap(
            rep1_peaks, rep2_peaks, pooled_peaks,
            pooledpr1_peaks, pooledpr2_peaks,
            rep1_ta, rep1_xcor, rep2_ta, rep2_xcor,
            paired_end, chrom_sizes, as_file, peak_type, prefix,
            fragment_length)

    # These are just passed through for convenience so that signals and tracks
    # are available in one place.  Both input and output are optional.
    if rep1_signal:
        output.update({"rep1_signal": rep1_signal})
    if rep2_signal:
        output.update({"rep2_signal": rep2_signal})
    if pooled_signal:
        output.update({"pooled_signal": pooled_signal})

    finalpeaks_dirname = rep1_peaks.split("/")[0]
    prefix = "final"
    peak_file = "%s/%s.%s.peaksfile" % (finalpeaks_dirname, prefix, peak_type,) 
    with open(peak_file,"w") as fh:
        for key, val in output.items():
            if isinstance(val, list):
                fh.write(": ".join([key, ", ".join(val)])+"\n")
            else:
                fh.write(": ".join([key, str(val)])+"\n")    

    logging.info("Exiting with output: %s", output)
    return output


if __name__ == '__main__':
    main(rep1_peaks, rep2_peaks, pooled_peaks,
         rep1_ta, rep1_xcor,
         paired_end, chrom_sizes, as_file, peak_type,
         pooledpr1_peaks, pooledpr2_peaks,
         rep2_ta, rep2_xcor,
         prefix,
         rep1_signal, rep2_signal, pooled_signal,
         fragment_length)
