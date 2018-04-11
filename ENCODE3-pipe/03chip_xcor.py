#!/usr/bin/env python3
#Description: ENCODE xcor 0.0.1 modified by WKL
#Environment: python3.6
#Author: WKL
#Date: 2018/04/07
#Examples: python3 chip_xcor.py --in1="R1.bam" --pair=False --spp="1.14"


import os
import getopt
import subprocess
import shlex
from multiprocessing import cpu_count
import dxpy
import logging
import sys
sys.path.insert(1,"/zs32_2/klwang/ChIP_seq/test_data/ENCODE_code")
import common
sys.path.remove("/zs32_2/klwang/ChIP_seq/test_data/ENCODE_code")


logger = logging.getLogger(__name__)
#logger.addHandler(dxpy.DXLogHandler())
logger.propagate = False
logger.setLevel(logging.INFO)


# SPP_VERSION_MAP = {
#     "1.10.1": "/phantompeakqualtools/spp_1.10.1.tar.gz",
#     "1.14":   "/phantompeakqualtools/spp-1.14.tar.gz"
# }

def usage():
    print("")
    print("Align fastq reads using bwa.")
    print("usage: python %s -option <argument>" %sys.argv[0])
    print("  -h/--help                   ")
    print("  --in1=<STRING>     [required] Bam input file.")
    print("  --pair=<True or False>        Paired_end or single_end.")
    print("  --spp=<STRING>                The spp version used.")

## deal with options
try:
    opts, args = getopt.getopt( sys.argv[1:], "h", ["help", "in1=", "pair=", "spp=" ] )
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
            input_bam = val
        if opt in ( "--pair", ):
            paired_end = val
        if opt in ( "--spp", ):
            spp_version = val
'''
try:
    debug = False if not debug else debug
except:
    debug = False
'''
if not "paired_end" in dir() or not paired_end or paired_end=="False":
    paired_end=False
else:
    paired_end=True


def strip_extensions(filename, extensions):
    basename = filename
    for extension in extensions:
        basename = basename.rpartition(extension)[0] or basename
    return basename


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


def main(input_bam, paired_end, spp_version):

    # The following line(s) initialize your data object inputs on the platform
    # into dxpy.DXDataObject instances that you can start using immediately.

    input_bam_file = input_bam

    input_bam_filename = input_bam_file
    input_bam_basename = input_bam_file.rstrip('.bam')

    intermediate_TA_filename = input_bam_basename + ".tagAlign"
    if paired_end:
        end_infix = 'PE2SE'
    else:
        end_infix = 'SE'
    final_TA_filename = input_bam_basename + '.' + end_infix + '.tagAlign.gz'

    # ===================
    # Create tagAlign file
    # ===================

    out, err = common.run_pipe([
        "bamToBed -i %s" % (input_bam_filename),
        r"""awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}'""",
        "tee %s" % (intermediate_TA_filename),
        "gzip -cn"],
        outfile=final_TA_filename)

    # ================
    # Create BEDPE file
    # ================
    if paired_end:
        final_BEDPE_filename = input_bam_basename + ".bedpe.gz"
        # need namesorted bam to make BEDPE
        final_nmsrt_bam_prefix = input_bam_basename + ".nmsrt"
        final_nmsrt_bam_filename = final_nmsrt_bam_prefix + ".bam"
        samtools_sort_command = \
            "samtools sort -n %s %s" % (input_bam_filename, final_nmsrt_bam_prefix)
        logger.info(samtools_sort_command)
        subprocess.check_output(shlex.split(samtools_sort_command))
        out, err = common.run_pipe([
            "bamToBed -bedpe -mate1 -i %s" % (final_nmsrt_bam_filename),
            "gzip -cn"],
            outfile=final_BEDPE_filename)

    # =================================
    # Subsample tagAlign file
    # ================================
    logger.info(
        "Intermediate tA md5: %s" % (common.md5(intermediate_TA_filename)))
    NREADS = 15000000
    if paired_end:
        end_infix = 'MATE1'
    else:
        end_infix = 'SE'
    subsampled_TA_filename = \
        input_bam_basename + \
        ".filt.nodup.sample.%d.%s.tagAlign.gz" % (NREADS/1000000, end_infix)
    steps = [
        'grep -v "chrM" %s' % (intermediate_TA_filename),
        'shuf -n %d --random-source=%s' % (NREADS, intermediate_TA_filename)]
    if paired_end:
        steps.extend([r"""awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}'"""])
    steps.extend(['gzip -cn'])
    out, err = common.run_pipe(steps, outfile=subsampled_TA_filename)
    logger.info(
        "Subsampled tA md5: %s" % (common.md5(subsampled_TA_filename)))

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

    tagAlign_file = final_TA_filename
    if paired_end:
        BEDPE_file = final_BEDPE_filename

    CC_scores_file = CC_scores_filename
    CC_plot_file = CC_plot_filename
    xcor_qc = xcor_parse(CC_scores_filename)

    # Return the outputs
    output = {
        "tagAlign_file": tagAlign_file,
        "CC_scores_file": CC_scores_file,
        "CC_plot_file": CC_plot_file,
        "paired_end": paired_end,
        "RSC": float(xcor_qc.get('relPhantomPeakCoef')),
        "NSC": float(xcor_qc.get('phantomPeakCoef')),
        "est_frag_len": float(xcor_qc.get('estFragLen'))
    }
    if paired_end:
        output.update({"BEDPE_file": BEDPE_file})

    parse_file = subsampled_TA_filename + ".parse"
    with open(parse_file,"w") as fh:
        for key, val in output.items():
            if isinstance(val, list):
                fh.write(": ".join([key, ", ".join(val)])+"\n")
            else:
                fh.write(": ".join([key, str(val)])+"\n")
    logger.info("Exiting with output:\n%s" % (pformat(output)))
    return output


'''
def main(input_bam, paired_end, spp_version):
    if paired_end:
        print("paired_end==True")
    else:
        print("paired_end==False")
'''

if __name__ == '__main__':
    main(input_bam, paired_end, spp_version)
