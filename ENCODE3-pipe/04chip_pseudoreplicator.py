#!/usr/bin/env python3
#Description: ENCODE pseudoreplicator 0.0.1 modified by WKL
#Author: WKL
#Environment: python3.6
#Date: 2018/04/08
#Example: python3 chip_xcor.py --in1="R1.tagAlign.gz" 

import subprocess
import re
import gzip
import os
import getopt
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
    print("Align fastq reads using bwa.")
    print("usage: python %s -option <argument>" %sys.argv[0])
    print("  -h/--help                   ")
    print("  --in1=<STRING>     [required] FASTA/FASTQ input file.")
    print("  --prefix=<STRING>             The prefix of the output files.")

## deal with options
try:
    opts, args = getopt.getopt( sys.argv[1:], "h", ["help", "in1=", "prefix=" ] )
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
            input_tags = val
        if opt in ( "--prefix", ):
            prefix = val
if not "prefix" in dir():
    prefix = None

def main(input_tags, prefix=None):
    input_tags_file = input_tags
    input_tags_filename = input_tags_file
    # introspect the file to determine tagAlign (thus SE) or BEDPE (thus PE)
    # strip extension as appropriate
    subprocess.check_output('ls', shell=True)
    with gzip.open(input_tags_filename) as f:
        firstline = f.readline()
    logger.info('First line of input_tags:\n%s' % (firstline))

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

    pr_ta_filenames = \
        [input_tags_basename + ".%s.pr1.tagAlign.gz" % (filename_infix),
         input_tags_basename + ".%s.pr2.tagAlign.gz" % (filename_infix)]

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

    pseudoreplicate1_file = pr_ta_filenames[0]
    pseudoreplicate2_file = pr_ta_filenames[1]

    output = {
        "pseudoreplicate1": pseudoreplicate1_file,
        "pseudoreplicate2": pseudoreplicate2_file
    }

    return output

if __name__ == '__main__':
    main(input_tags, prefix)