#!/usr/bin/env python
#Description: encode_macs2 0.0.1 modified by WKL
#Author: WKL
#Environment: python3.6
#Date: 2018/04/08
#Example: python3 chip_macs2.py --in1="R1.tagAlign.gz" --ct1="C1.tagAlign.gz" --xcor1="R1.cc.qc" --in1_PE=False \
#    --chromsizes="encode_db/male.hg19.chrom.sizes" --genomesize="hs" --narrowpeak="encode_db/narrowPeak.as" \
#    --gappedpeak="encode_db/gappedPeak.as" --broadpeak="encode_db/broadPeak.as"

import os
import re
import gzip
import subprocess
import getopt
import shlex
from multiprocessing import cpu_count
from os.path import splitext
import logging
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
    print("Call peaks using MACS2.")
    print("usage: python %s -option <argument>" %sys.argv[0])
    print("  -h/--help                   ")
    print("  --in1=<STRING>              [required] Rep1 tagAlign input file.")
    print("  --in2=<STRING>                         Rep2 tagAlign input file.")    
    print("  --ct1=<STRING>              [required] Ctl1 tagAlign input file.")
    print("  --ct2=<STRING>                         Ctl2 tagAlign input file.")    
    print("  --xcor1=<STRING>            [required] Rep1 xcor input file.")
    print("  --xcor2=<STRING>                       Rep2 xcor input file.")    
    print("  --in1_PE=<True or False>    [required] Rep1 paired_end or single_end.")
    print("  --in2_PE=<True or False>               Rep2 paired_end or single_end.")
    print("  --chromsizes=<STRING>                  The chromosome sizes file.")
    print("  --genomesize=<STRING>                  The species such as 'hs'.")
    print("  --narrowpeak=<STRING>                  The narrowpeak_as file.")
    print("  --gappedpeak=<STRING>                  The gappedpeak_as file.")
    print("  --broadpeak=<STRING>                   The broadpeak_as file.")
    print("  --fragment=<INT>                       The fragment length.")

## deal with options
try:
    opts, args = getopt.getopt( sys.argv[1:], "h", ["help", "in1=", "in2=", "ct1=", "ct2=", "xcor1=", "xcor2=", "in1_PE=", "in2_PE=", "chromsizes=", "genomesize=", "narrowpeak=", "gappedpeak=", "broadpeak=", "fragment="] )
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
        if opt in ( "--ct1", ):
            ctl1_ta = val
        if opt in ( "--ct2", ):
            ctl2_ta = val
        if opt in ( "--xcor1", ):
            rep1_xcor = val
        if opt in ( "--xcor2", ):
            rep2_xcor = val
        if opt in ( "--in1_PE", ):
            rep1_paired_end = val
        if opt in ( "--in2_PE", ):
            rep2_paired_end = val
        if opt in ( "--chromsizes", ):
            chrom_sizes = val
        if opt in ( "--genomesize", ):
            genomesize = val
        if opt in ( "--narrowpeak", ):
            narrowpeak_as = val
        if opt in ( "--gappedpeak", ):
            gappedpeak_as = val
        if opt in ( "--broadpeak", ):
            broadpeak_as = val
        if opt in ( "--fragment", ):
            fragment_length = int(val)
if not "rep2_ta" in dir():
    rep2_ta = None
if not "ctl2_ta" in dir():
    ctl2_ta = None
if not "rep2_xcor" in dir():
    rep2_xcor = None
if not "rep2_paired_end" in dir():
    rep2_paired_end = None
if not "fragment_length" in dir():
    fragment_length = None
'''
try:
    debug = False if not debug else debug
except:
    debug = False
'''
if not rep1_paired_end or rep1_paired_end=="False":
    rep1_paired_end=False
else:
    rep1_paired_end=True

def macs2(experiment, control, xcor_scores_input, chrom_sizes,
         narrowpeak_as, gappedpeak_as, broadpeak_as, genomesize, prefix=None,
         fragment_length=None):
    # Define the output filenames
    peaks_dirname = 'peaks_macs'
    if not os.path.exists(peaks_dirname):
        os.makedirs(peaks_dirname)
    if not prefix:
        prefix = experiment
    if prefix.endswith('.gz'):
        prefix = prefix[:-3]
    narrowPeak_fn    = "%s/%s.narrowPeak" % (peaks_dirname, prefix)
    gappedPeak_fn    = "%s/%s.gappedPeak" % (peaks_dirname, prefix)
    broadPeak_fn     = "%s/%s.broadPeak"  % (peaks_dirname, prefix)
    narrowPeak_gz_fn = narrowPeak_fn + ".gz"
    gappedPeak_gz_fn = gappedPeak_fn + ".gz"
    broadPeak_gz_fn  = broadPeak_fn  + ".gz"
    fc_signal_fn     = "%s/%s.fc_signal.bw"     % (peaks_dirname, prefix)
    pvalue_signal_fn = "%s/%s.pvalue_signal.bw" % (peaks_dirname, prefix)
    # Extract the fragment length estimate from column 3 of the
    # cross-correlation scores file
    # if the fragment_length argument is given, use that instead
    if fragment_length is not None:
        fraglen = str(fragment_length)
        logger.info("User given fragment length %s" % fraglen)
    else:
        with open(xcor_scores_input, 'r') as fh:
            firstline = fh.readline()
            fraglen = firstline.split()[2]  # third column
            logger.info("Fraglen %s" % (fraglen))
    # ===========================================
    # Generate narrow peaks and preliminary signal tracks
    # ============================================
    command = 'macs2 callpeak ' + \
              '-t %s -c %s ' % (experiment, control) + \
              '-f BED -n %s/%s ' % (peaks_dirname, prefix) + \
              '-g %s -p 1e-2 --nomodel --shift 0 --extsize %s --keep-dup all -B --SPMR' % (genomesize, fraglen)
    logger.info(command)
    returncode = common.block_on(command)
    logger.info("MACS2 exited with returncode %d" % (returncode))
    assert returncode == 0, "MACS2 non-zero return"
    # MACS2 sometimes calls features off the end of chromosomes.  Fix that.
    clipped_narrowpeak_fn = common.slop_clip(
        '%s/%s_peaks.narrowPeak' % (peaks_dirname, prefix),
        chrom_sizes)
    # Rescale Col5 scores to range 10-1000 to conform to narrowPeak.as format
    # (score must be <1000)
    rescaled_narrowpeak_fn = common.rescale_scores(
        clipped_narrowpeak_fn,
        scores_col=5)
    # Sort by Col8 in descending order and replace long peak names in Column 4
    # with Peak_<peakRank>
    pipe = ['sort -k 8gr,8gr %s' % (rescaled_narrowpeak_fn),
            r"""awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}'""",
            'tee %s' % (narrowPeak_fn),
            'gzip -cn']
    out, err = common.run_pipe(pipe, '%s' % (narrowPeak_gz_fn))
    # remove additional files
    #rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.xls ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.bed ${peakFile}_summits.bed
    # ===========================================
    # Generate Broad and Gapped Peaks
    # ============================================
    command = 'macs2 callpeak ' + \
              '-t %s -c %s ' % (experiment, control) + \
              '-f BED -n %s/%s ' % (peaks_dirname, prefix) + \
              '-g %s -p 1e-2 --broad --nomodel --shift 0 --extsize %s --keep-dup all' % (genomesize, fraglen)
    logger.info(command)
    returncode = common.block_on(command)
    logger.info("MACS2 exited with returncode %d" % (returncode))
    assert returncode == 0, "MACS2 non-zero return"
    # MACS2 sometimes calls features off the end of chromosomes.  Fix that.
    clipped_broadpeak_fn = common.slop_clip(
        '%s/%s_peaks.broadPeak' % (peaks_dirname, prefix),
        chrom_sizes)
    # Rescale Col5 scores to range 10-1000 to conform to narrowPeak.as format
    # (score must be <1000)
    rescaled_broadpeak_fn = common.rescale_scores(
        clipped_broadpeak_fn,
        scores_col=5)
    # Sort by Col8 (for broadPeak) or Col 14(for gappedPeak) in descending
    # order and replace long peak names in Column 4 with Peak_<peakRank>
    pipe = ['sort -k 8gr,8gr %s' % (rescaled_broadpeak_fn),
            r"""awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}'""",
            'tee %s' % (broadPeak_fn),
            'gzip -cn']
    out, err = common.run_pipe(pipe, '%s' % (broadPeak_gz_fn))
    # MACS2 sometimes calls features off the end of chromosomes.  Fix that.
    clipped_gappedpeaks_fn = common.slop_clip(
        '%s/%s_peaks.gappedPeak' % (peaks_dirname, prefix),
        chrom_sizes,
        bed_type='gappedPeak')
    # Rescale Col5 scores to range 10-1000 to conform to narrowPeak.as format
    # (score must be <1000)
    rescaled_gappedpeak_fn = common.rescale_scores(
        clipped_gappedpeaks_fn,
        scores_col=5)
    pipe = ['sort -k 14gr,14gr %s' % (rescaled_gappedpeak_fn),
            r"""awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}'""",
            'tee %s' % (gappedPeak_fn),
            'gzip -cn']
    out, err = common.run_pipe(pipe, '%s' % (gappedPeak_gz_fn))
    # remove additional files
    #rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.xls ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_peaks.bed ${peakFile}_summits.bed
    # ===========================================
    # For Fold enrichment signal tracks
    # ============================================
    # This file is a tab delimited file with 2 columns Col1 (chromosome name),
    # Col2 (chromosome size in bp).
    command = 'macs2 bdgcmp ' + \
              '-t %s/%s_treat_pileup.bdg ' % (peaks_dirname, prefix) + \
              '-c %s/%s_control_lambda.bdg ' % (peaks_dirname, prefix) + \
              '--outdir %s -o %s_FE.bdg ' % (peaks_dirname, prefix) + \
              '-m FE'
    logger.info(command)
    returncode = common.block_on(command)
    logger.info("MACS2 exited with returncode %d" % (returncode))
    assert returncode == 0, "MACS2 non-zero return"
    # Remove coordinates outside chromosome sizes (stupid MACS2 bug)
    pipe = ['slopBed -i %s/%s_FE.bdg -g %s -b 0' % (peaks_dirname, prefix, chrom_sizes),
            'bedClip stdin %s %s/%s.fc.signal.bedgraph' % (chrom_sizes, peaks_dirname, prefix)]
    out, err = common.run_pipe(pipe)
    #rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_FE.bdg
    # Convert bedgraph to bigwig
    command = 'bedGraphToBigWig ' + \
              '%s/%s.fc.signal.bedgraph ' % (peaks_dirname, prefix) + \
              '%s ' % (chrom_sizes) + \
              '%s' % (fc_signal_fn)
    logger.info(command)
    returncode = common.block_on(command)
    logger.info("bedGraphToBigWig exited with returncode %d" % (returncode))
    assert returncode == 0, "bedGraphToBigWig non-zero return"
    #rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.fc.signal.bedgraph
    # ===========================================
    # For -log10(p-value) signal tracks
    # ============================================
    # Compute sval =
    # min(no. of reads in ChIP, no. of reads in control) / 1,000,000
    out, err = common.run_pipe([
        'gzip -dc %s' % (experiment),
        'wc -l'])
    chipReads = out.strip()
    out, err = common.run_pipe([
        'gzip -dc %s' % (control),
        'wc -l'])
    controlReads = out.strip()
    sval = str(min(float(chipReads), float(controlReads))/1000000)
    logger.info(
        "chipReads = %s, controlReads = %s, sval = %s"
        % (chipReads, controlReads, sval))
    returncode = common.block_on(
        'macs2 bdgcmp ' +
        '-t %s/%s_treat_pileup.bdg ' % (peaks_dirname, prefix) +
        '-c %s/%s_control_lambda.bdg ' % (peaks_dirname, prefix) +
        '--outdir %s -o %s_ppois.bdg ' % (peaks_dirname, prefix) +
        '-m ppois -S %s' % (sval))
    logger.info("MACS2 exited with returncode %d" % (returncode))
    assert returncode == 0, "MACS2 non-zero return"
    # Remove coordinates outside chromosome sizes (stupid MACS2 bug)
    pipe = ['slopBed -i %s/%s_ppois.bdg -g %s -b 0' % (peaks_dirname, prefix, chrom_sizes),
            'bedClip stdin %s %s/%s.pval.signal.bedgraph' % (chrom_sizes, peaks_dirname, prefix)]
    out, err = common.run_pipe(pipe)
    #rm -rf ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_ppois.bdg
    # Convert bedgraph to bigwig
    command = 'bedGraphToBigWig ' + \
              '%s/%s.pval.signal.bedgraph ' % (peaks_dirname, prefix) + \
              '%s ' % (chrom_sizes) + \
              '%s' % (pvalue_signal_fn)
    logger.info(command)
    returncode = common.block_on(command)
    logger.info("bedGraphToBigWig exited with returncode %d" % (returncode))
    assert returncode == 0, "bedGraphToBigWig non-zero return"
    #rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.pval.signal.bedgraph
    #rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_treat_pileup.bdg ${peakFile}_control_lambda.bdg
    # ===========================================
    # Generate bigWigs from beds to support trackhub visualization of peak files
    # ============================================
    narrowPeak_bb_fname = common.bed2bb('%s' %(narrowPeak_fn), chrom_sizes, narrowpeak_as, bed_type='bed6+4')
    gappedPeak_bb_fname = common.bed2bb('%s' %(gappedPeak_fn), chrom_sizes, gappedpeak_as, bed_type='bed12+3')
    broadPeak_bb_fname =  common.bed2bb('%s' %(broadPeak_fn),  chrom_sizes, broadpeak_as,  bed_type='bed6+3')
    #Temporary during development to create empty files just to get the applet to exit 
    # narrowPeak_bb_fn = "%s.bb" % (narrowPeak_fn)
    # gappedPeak_bb_fn = "%s.bb" % (gappedPeak_fn)
    # broadPeak_bb_fn  = "%s.bb" % (broadPeak_fn)
    #
    # Upload the file outputs
    #
    narrowPeak    = narrowPeak_gz_fn
    gappedPeak    = gappedPeak_gz_fn
    broadPeak     = broadPeak_gz_fn
    narrowPeak_bb = narrowPeak_bb_fname
    gappedPeak_bb = gappedPeak_bb_fname
    broadPeak_bb  = broadPeak_bb_fname
    fc_signal     = fc_signal_fn
    pvalue_signal = pvalue_signal_fn
    #
    # Build the output structure.
    #
    output = {
        "narrowpeaks":    narrowPeak,
        "gappedpeaks":    gappedPeak,
        "broadpeaks":     broadPeak,
        "narrowpeaks_bb": narrowPeak_bb,
        "gappedpeaks_bb": gappedPeak_bb,
        "broadpeaks_bb":  broadPeak_bb,
        "fc_signal":     fc_signal,
        "pvalue_signal": pvalue_signal
    }
    return output


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

def pseudoreplicator(input_tags, prefix=None):
    input_tags_file = input_tags
    input_tags_filename = input_tags_file
    # introspect the file to determine tagAlign (thus SE) or BEDPE (thus PE)
    # strip extension as appropriate
    subprocess.check_output('ls', shell=True)
    with gzip.open(input_tags_filename) as f:
        firstline = f.readline()
    logger.info('First line of input_tags:\n%s' % (firstline))
    #
    se_cols = 6
    pe_cols = 10
    firstline = firstline.decode("utf-8")
    if re.match('^(\S+[\t\n]){%d}$' % (se_cols), firstline):
        paired_end = False
        input_tags_basename = prefix or input_tags_filename.rstrip('.tagAlign.gz')
        filename_infix = 'SE'
        logger.info("Detected single-end data")
    elif re.match('^(\S+[\t\n]){%d}$' % (pe_cols), firstline):
        paired_end = True
        input_tags_basename = prefix or input_tags_filename.rstrip('.bedpe.gz')
        filename_infix = 'PE2SE'
        logger.info("Detected paired-end data")
    else:
        raise IOError(
            "%s is neither a BEDPE or tagAlign file" % (input_tags_filename))
    #
    pr_ta_filenames = \
        [input_tags_basename + ".%s.pr1.tagAlign.gz" % (filename_infix),
         input_tags_basename + ".%s.pr2.tagAlign.gz" % (filename_infix)]
    #
    # count lines in the file
    out, err = common.run_pipe([
        'gzip -dc %s' % (input_tags_filename),
        'wc -l'])
    # number of lines in each split
    nlines = (int(out)+1)/2
    # Shuffle and split BEDPE file into 2 equal parts
    # by using the input to seed shuf we ensure multiple runs with the same
    # input will produce the same output
    # Produces two files named splits_prefix0n, n=1,2
    splits_prefix = 'temp_split'
    out, err = common.run_pipe([
        'gzip -dc %s' % (input_tags_filename),
        'shuf --random-source=%s' % (input_tags_filename),
        'split -a 2 -d -l %d - %s' % (nlines, splits_prefix)])
    # Convert read pairs to reads into standard tagAlign file
    for i, index in enumerate(['00', '01']):  # could be made multi-threaded
        steps = ['cat %s' % (splits_prefix+index)]
        if paired_end:
            steps.extend([r"""awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}'"""])
        steps.extend(['gzip -cn'])
        out, err = common.run_pipe(steps, outfile=pr_ta_filenames[i])
    #
    pseudoreplicate1_file = pr_ta_filenames[0]
    pseudoreplicate2_file = pr_ta_filenames[1]
    #
    output = {
        "pseudoreplicate1": pseudoreplicate1_file,
        "pseudoreplicate2": pseudoreplicate2_file
    }
    return output

def main(rep1_ta, ctl1_ta, rep1_xcor, rep1_paired_end, chrom_sizes, genomesize,
         narrowpeak_as, gappedpeak_as, broadpeak_as,
         rep2_ta=None, ctl2_ta=None, rep2_xcor=None, rep2_paired_end=None,
         fragment_length=None):
    rep1_ta_file = rep1_ta
    rep1_ta_filename = rep1_ta_file
    ntags_rep1 = common.count_lines(rep1_ta_filename)
    #
    simplicate_experiment = rep1_ta and not rep2_ta
    if simplicate_experiment:
        logger.info("No rep2 tags specified so processing as a simplicate experiment.")
    else:
        logger.info("Rep1 and rep2 tags specified so processing as a replicated experiment.")
    #
    if not simplicate_experiment:
        assert rep1_paired_end == rep2_paired_end, 'Mixed PE/SE not supported'
        rep2_ta_file = rep2_ta
        rep2_ta_filename = rep2_ta_file
        ntags_rep2 = common.count_lines(rep2_ta_filename)
    paired_end = rep1_paired_end
    #
    unary_control = (ctl1_ta == ctl2_ta) or not ctl2_ta
    ctl1_ta_file = ctl1_ta
    ctl1_ta_filename = ctl1_ta_file
    #
    if not unary_control:
        ctl2_ta_file = ctl2_ta
        ctl2_ta_filename = ctl2_ta_file
    else:
        ctl2_ta_file = ctl1_ta_file
        ctl2_ta_filename = ctl1_ta_file
    #
    ntags_ctl1 = common.count_lines(ctl1_ta_filename)
    ntags_ctl2 = common.count_lines(ctl2_ta_filename)
    rep1_control = ctl1_ta  # default.  May be changed later.
    rep1_ctl_msg = "control rep1"
    rep2_control = ctl2_ta  # default.  May be changed later.
    rep2_ctl_msg = "control rep2"
    #
    rep_info = [(ntags_rep1, 'replicate 1', rep1_ta_filename)]
    if not simplicate_experiment:
        rep_info.append((ntags_rep2, 'replicate 2', rep2_ta_filename))
    rep_info.extend(
        [(ntags_ctl1, 'control 1', ctl1_ta_filename),
         (ntags_ctl2, 'control 2', ctl2_ta_filename)])
    for n, name, filename in rep_info:
        logger.info("Found %d tags in %s file %s" % (n, name, filename))
    #
    subprocess.check_output('ls -l', shell=True, stderr=subprocess.STDOUT)
    #
    if not simplicate_experiment:
        #Pool replicates
        pool_replicates_subjob = pool(inputs=[rep1_ta, rep2_ta], prefix='pooled_reps')
        pooled_replicates = pool_replicates_subjob.get("pooled")
        #Pool cross-correlation
        pooled_replicates_xcor_subjob = xcor_only(pooled_replicates, paired_end)
    #
    if unary_control:
        logger.info("Only one control supplied.")
        if not simplicate_experiment:
            logger.info("Using one control for both replicate 1 and 2 and for the pool.")
        rep2_control = rep1_control
        control_for_pool = rep1_control
        pool_ctl_msg = "one control"
    else:
        #Pool controls
        pool_controls_subjob = pool(inputs=[ctl1_ta, ctl2_ta],prefix="PL_ctls")
        pooled_controls = pool_controls_subjob.get("pooled")
        # always use the pooled controls for the pool
        control_for_pool = pooled_controls
        pool_ctl_msg = "pooled controls"
        # use the pooled controls for the reps depending on the ratio of rep to
        # control reads
        ratio_ctl_reads = float(ntags_ctl1)/float(ntags_ctl2)
        if ratio_ctl_reads < 1:
                ratio_ctl_reads = 1/ratio_ctl_reads
        ratio_cutoff = 1.2
        if ratio_ctl_reads > ratio_cutoff:
                logger.info(
                    "Number of reads in controls differ by > factor of %f. Using pooled controls."
                    % (ratio_cutoff))
                rep1_control = pooled_controls
                rep2_control = pooled_controls
        else:
                if ntags_ctl1 < ntags_rep1:
                        logger.info("Fewer reads in control replicate 1 than experiment replicate 1.  Using pooled controls for replicate 1.")
                        rep1_control = pooled_controls
                        rep1_ctl_msg = "pooled controls"
                elif not simplicate_experiment and ntags_ctl2 < ntags_rep2:
                        logger.info("Fewer reads in control replicate 2 than experiment replicate 2.  Using pooled controls for replicate 2.")
                        rep2_control = pooled_controls
                        rep2_ctl_msg = "pooled controls"
                else:
                    logger.info(
                        "Using distinct controls for replicate 1 and 2.")
                    rep1_control = ctl1_ta  # default.  May be changed later.
                    rep2_control = ctl2_ta  # default.  May be changed later.
                    rep1_ctl_msg = "control rep1"
                    rep2_ctl_msg = "control rep2"
    #
    rep1_pr_subjob = pseudoreplicator(input_tags=rep1_ta)
    if not simplicate_experiment:
        rep2_pr_subjob = pseudoreplicator(input_tags=rep2_ta)
    #
        pool_pr1_subjob = pool(inputs=[rep1_pr_subjob.get("pseudoreplicate1"),rep2_pr_subjob.get("pseudoreplicate1")],prefix='PPR1')
        pool_pr2_subjob = pool(inputs=[rep1_pr_subjob.get("pseudoreplicate2"),rep2_pr_subjob.get("pseudoreplicate2")],prefix='PPR2')
    #
    common_args = {
        'chrom_sizes':      chrom_sizes,
        'genomesize':       genomesize,
        'narrowpeak_as':    narrowpeak_as,
        'gappedpeak_as':    gappedpeak_as,
        'broadpeak_as':     broadpeak_as
        }
    # if the fragment_length argument is given, update macs2 input
    if fragment_length is not None:
        common_args.update({'fragment_length' : fragment_length})
    #macs2(experiment, control, xcor_scores_input, chrom_sizes,narrowpeak_as, gappedpeak_as, broadpeak_as, genomesize, prefix=None,fragment_length=None)
    common_args.update({'prefix': 'r1'})
    rep1_peaks_subjob      = macs2( rep1_ta,
                                    rep1_control,
                                    rep1_xcor, **common_args)
    #
    common_args.update({'prefix': 'r1pr1'})
    rep1pr1_peaks_subjob   = macs2( rep1_pr_subjob.get("pseudoreplicate1"),
                                    rep1_control,
                                    rep1_xcor, **common_args)
    #
    common_args.update({'prefix': 'r1pr2'})
    rep1pr2_peaks_subjob   = macs2( rep1_pr_subjob.get("pseudoreplicate2"),
                                    rep1_control,
                                    rep1_xcor, **common_args)
    #
    if not simplicate_experiment:
        common_args.update({'prefix': 'r2'})
        rep2_peaks_subjob      = macs2( rep2_ta,
                                        rep2_control,
                                        rep2_xcor, **common_args)
    #
        common_args.update({'prefix': 'r2pr1'})
        rep2pr1_peaks_subjob   = macs2( rep2_pr_subjob.get("pseudoreplicate1"),
                                        rep2_control,
                                        rep2_xcor, **common_args)
    #
        common_args.update({'prefix': 'r2pr2'})
        rep2pr2_peaks_subjob   = macs2( rep2_pr_subjob.get("pseudoreplicate2"),
                                        rep2_control,
                                        rep2_xcor, **common_args)
    #
        common_args.update({'prefix': 'pool'})
        pooled_peaks_subjob    = macs2( pooled_replicates,
                                        control_for_pool,   
                                        pooled_replicates_xcor_subjob.get("CC_scores_file"), **common_args)
    #
        common_args.update({'prefix': 'ppr1'})
        pooledpr1_peaks_subjob = macs2( pool_pr1_subjob.get("pooled"),
                                        control_for_pool,
                                        pooled_replicates_xcor_subjob.get("CC_scores_file"), **common_args)
    #
        common_args.update({'prefix': 'ppr2'})
        pooledpr2_peaks_subjob = macs2( pool_pr2_subjob.get("pooled"),
                                        control_for_pool,
                                        pooled_replicates_xcor_subjob.get("CC_scores_file"), **common_args)
    #
    output = {
        'rep1_narrowpeaks':         rep1_peaks_subjob.get("narrowpeaks"),
        'rep1_gappedpeaks':         rep1_peaks_subjob.get("gappedpeaks"),
        'rep1_broadpeaks':          rep1_peaks_subjob.get("broadpeaks"),
        'rep1_narrowpeaks_bb':      rep1_peaks_subjob.get("narrowpeaks_bb"),
        'rep1_gappedpeaks_bb':      rep1_peaks_subjob.get("gappedpeaks_bb"),
        'rep1_broadpeaks_bb':       rep1_peaks_subjob.get("broadpeaks_bb"),
        'rep1_fc_signal':           rep1_peaks_subjob.get("fc_signal"),
        'rep1_pvalue_signal':       rep1_peaks_subjob.get("pvalue_signal"),
        # 
        'rep1pr1_narrowpeaks':      rep1pr1_peaks_subjob.get("narrowpeaks"),
        'rep1pr1_gappedpeaks':      rep1pr1_peaks_subjob.get("gappedpeaks"),
        'rep1pr1_broadpeaks':       rep1pr1_peaks_subjob.get("broadpeaks"),
        'rep1pr1_fc_signal':        rep1pr1_peaks_subjob.get("fc_signal"),
        'rep1pr1_pvalue_signal':    rep1pr1_peaks_subjob.get("pvalue_signal"),
        #
        'rep1pr2_narrowpeaks':      rep1pr2_peaks_subjob.get("narrowpeaks"),
        'rep1pr2_gappedpeaks':      rep1pr2_peaks_subjob.get("gappedpeaks"),
        'rep1pr2_broadpeaks':       rep1pr2_peaks_subjob.get("broadpeaks"),
        'rep1pr2_fc_signal':        rep1pr2_peaks_subjob.get("fc_signal"),
        'rep1pr2_pvalue_signal':    rep1pr2_peaks_subjob.get("pvalue_signal")
    }
    #
    if not simplicate_experiment:
        output.update({
            'rep2_narrowpeaks':         rep2_peaks_subjob.get("narrowpeaks"),
            'rep2_gappedpeaks':         rep2_peaks_subjob.get("gappedpeaks"),
            'rep2_broadpeaks':          rep2_peaks_subjob.get("broadpeaks"),
            'rep2_narrowpeaks_bb':      rep2_peaks_subjob.get("narrowpeaks_bb"),
            'rep2_gappedpeaks_bb':      rep2_peaks_subjob.get("gappedpeaks_bb"),
            'rep2_broadpeaks_bb':       rep2_peaks_subjob.get("broadpeaks_bb"),
            'rep2_fc_signal':           rep2_peaks_subjob.get("fc_signal"),
            'rep2_pvalue_signal':       rep2_peaks_subjob.get("pvalue_signal"),
            #
            'rep2pr1_narrowpeaks':      rep2pr1_peaks_subjob.get("narrowpeaks"),
            'rep2pr1_gappedpeaks':      rep2pr1_peaks_subjob.get("gappedpeaks"),
            'rep2pr1_broadpeaks':       rep2pr1_peaks_subjob.get("broadpeaks"),
            'rep2pr1_fc_signal':        rep2pr1_peaks_subjob.get("fc_signal"),
            'rep2pr1_pvalue_signal':    rep2pr1_peaks_subjob.get("pvalue_signal"),
            #
            'rep2pr2_narrowpeaks':      rep2pr2_peaks_subjob.get("narrowpeaks"),
            'rep2pr2_gappedpeaks':      rep2pr2_peaks_subjob.get("gappedpeaks"),
            'rep2pr2_broadpeaks':       rep2pr2_peaks_subjob.get("broadpeaks"),
            'rep2pr2_fc_signal':        rep2pr2_peaks_subjob.get("fc_signal"),
            'rep2pr2_pvalue_signal':    rep2pr2_peaks_subjob.get("pvalue_signal"),
            #
            'pooled_narrowpeaks':       pooled_peaks_subjob.get("narrowpeaks"),
            'pooled_gappedpeaks':       pooled_peaks_subjob.get("gappedpeaks"),
            'pooled_broadpeaks':        pooled_peaks_subjob.get("broadpeaks"),
            'pooled_narrowpeaks_bb':    pooled_peaks_subjob.get("narrowpeaks_bb"),
            'pooled_gappedpeaks_bb':    pooled_peaks_subjob.get("gappedpeaks_bb"),
            'pooled_broadpeaks_bb':     pooled_peaks_subjob.get("broadpeaks_bb"),
            'pooled_fc_signal':         pooled_peaks_subjob.get("fc_signal"),
            'pooled_pvalue_signal':     pooled_peaks_subjob.get("pvalue_signal"),
            #
            'pooledpr1_narrowpeaks':    pooledpr1_peaks_subjob.get("narrowpeaks"),
            'pooledpr1_gappedpeaks':    pooledpr1_peaks_subjob.get("gappedpeaks"),
            'pooledpr1_broadpeaks':     pooledpr1_peaks_subjob.get("broadpeaks"),
            'pooledpr1_fc_signal':      pooledpr1_peaks_subjob.get("fc_signal"),
            'pooledpr1_pvalue_signal':  pooledpr1_peaks_subjob.get("pvalue_signal"),
            #
            'pooledpr2_narrowpeaks':    pooledpr2_peaks_subjob.get("narrowpeaks"),
            'pooledpr2_gappedpeaks':    pooledpr2_peaks_subjob.get("gappedpeaks"),
            'pooledpr2_broadpeaks':     pooledpr2_peaks_subjob.get("broadpeaks"),
            'pooledpr2_fc_signal':      pooledpr2_peaks_subjob.get("fc_signal"),
            'pooledpr2_pvalue_signal':  pooledpr2_peaks_subjob.get("pvalue_signal")
        })
    peaks_dirname = 'peaks_macs'
    prefix = rep1_ta_filename.split("/")[-1]
    peak_file = "%s/%s.peaksfile" % (peaks_dirname, prefix) 
    with open(peak_file,"w") as fh:
        for key, val in output.items():
            if isinstance(val, list):
                fh.write(": ".join([key, ", ".join(val)])+"\n")
            else:
                fh.write(": ".join([key, str(val)])+"\n")
    return output

if __name__ == '__main__':
    main(rep1_ta, ctl1_ta, rep1_xcor, rep1_paired_end, chrom_sizes, genomesize,
         narrowpeak_as, gappedpeak_as, broadpeak_as,
         rep2_ta, ctl2_ta, rep2_xcor, rep2_paired_end,
         fragment_length)