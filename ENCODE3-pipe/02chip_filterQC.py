#!/usr/bin/env python3
#Description: filter_qc 0.0.1 modified by WKL
#Author: WKL
#Environment: python3.6
#Date: 2018/04/04
#Examples: python3 chip_filterQC.py  --in="R1.bam" --pair=False, --scrub=False, --picard="1.119" --samtools_params="-q 30" --samtools="1.2" --debug=False

from __future__ import print_function
import os
import re
import subprocess
import shlex
import getopt
#import dxpy
import logging
from multiprocessing import Pool
from pprint import pprint, pformat
import sys
sys.path.insert(1,"/zs32_2/klwang/ChIP_seq/test_data/ENCODE_code")
import common
sys.path.remove("/zs32_2/klwang/ChIP_seq/test_data/ENCODE_code")


logger = logging.getLogger(__name__)
#logger.addHandler(dxpy.DXLogHandler())
logger.propagate = False
logger.setLevel(logging.INFO)

SAMTOOLS_PATH = {
    "1.3.1": "/zs32/home/klwang/softwares/samtools-1.3.1/bin/samtools",
    "1.2": "/opt/tools/seq-analysis/samtools-1.2/samtools"
}

PICARD_PATH = {
    "1.119": "/opt/tools/seq-analysis/picard-tools-1.119/MarkDuplicates.jar"
}

def usage():
    print("")
    print("Filter and QC of bam file.")
    print("usage: python %s -option <argument>" %sys.argv[0])
    print("  -h/--help                   ")
    print("  --in=<STRING>                 Bam input file.")
    print("  --pair=<True or False>        Paired_end or singled_end.")
    print("  --scrub=<True or False>       Scrub or not for the bam file.")
    print("  --picard=<STRING>             The picard version used.")
    print("  --samtools_params=<STRING>    The samtools parameters used.")
    print("  --samtools=<STRING>           The samtools version used.")
    print("  --debug=<True or False>       Whether debug or not.")

## deal with options
try:
    opts, args = getopt.getopt( sys.argv[1:], "h", ["help", "in=", "pair=", "scrub=", "picard=", "samtools_params=", "samtools=", "debug=" ] )
except getopt.GetoptError:
    print("ERROR: Get option error.")
    usage()
    sys.exit(2)

for opt, val in opts:
    if opt in ( "-h", "--help" ):
        usage()
        sys.exit(1)
    else:
        if opt in ( "--in", ):
            input_bam = val
        if opt in ( "--pair", ):
            paired_end = val
        if opt in ( "--scrub", ):
            scrub = val
        if opt in ( "--picard", ):
            picard_version = val
        if opt in ( "--samtools_params", ):
            samtools_params = val
        if opt in ( "--samtools", ):
            samtools_version = val
        if opt in ( "--debug", ):
            debug = val

if not "debug" in dir() or not debug or debug=="False":
    debug=False
else:
    debug=True
if not "paired_end" in dir() or not paired_end or paired_end=="False":
    paired_end=False
else:
    paired_end=True
if not "scrub" in dir() or not scrub or scrub=="False":
    scrub=False
else:
    scrub=True

def dup_parse(fname):
    with open(fname, 'r') as dup_file:
        if not dup_file:
            return None
        lines = iter(dup_file.read().splitlines())
        for line in lines:
            if line.startswith('## METRICS CLASS'):
                headers = lines.__next__().rstrip('\n').lower()
                metrics = lines.__next__().rstrip('\n')
                break
        headers = headers.split('\t')
        metrics = metrics.split('\t')
        headers.pop(0)
        metrics.pop(0)
        dup_qc = dict(zip(headers, metrics))
    return dup_qc


def pbc_parse(fname):
    with open(fname, 'r') as pbc_file:
        if not pbc_file:
            return None
        lines = pbc_file.read().splitlines()
        line = lines[0].rstrip('\n')
        # PBC File output:
        #   TotalReadPairs <tab>
        #   DistinctReadPairs <tab>
        #   OneReadPair <tab>
        #   TwoReadPairs <tab>
        #   NRF=Distinct/Total <tab>
        #   PBC1=OnePair/Distinct <tab>
        #   PBC2=OnePair/TwoPair
        headers = ['TotalReadPairs',
                   'DistinctReadPairs',
                   'OneReadPair',
                   'TwoReadPairs',
                   'NRF',
                   'PBC1',
                   'PBC2']
        metrics = line.split('\t')
        pbc_qc = dict(zip(headers, metrics))
    return pbc_qc


def flagstat_parse(fname):
    with open(fname, 'r') as flagstat_file:
        if not flagstat_file:
            return None
        flagstat_lines = flagstat_file.read().splitlines()
    qc_dict = {
        # values are regular expressions,
        # will be replaced with scores [hiq, lowq]
        'in_total': 'in total',
        'duplicates': 'duplicates',
        'mapped': 'mapped',
        'paired_in_sequencing': 'paired in sequencing',
        'read1': 'read1',
        'read2': 'read2',
        'properly_paired': 'properly paired',
        'with_self_mate_mapped': 'with itself and mate mapped',
        'singletons': 'singletons',
        # i.e. at the end of the line
        'mate_mapped_different_chr': 'with mate mapped to a different chr$',
        # RE so must escape
        'mate_mapped_different_chr_hiQ':
            'with mate mapped to a different chr \(mapQ>=5\)'
    }
    for (qc_key, qc_pattern) in qc_dict.items():
        qc_metrics = next(re.split(qc_pattern, line)
                          for line in flagstat_lines
                          if re.search(qc_pattern, line))
        (hiq, lowq) = qc_metrics[0].split(' + ')
        qc_dict[qc_key] = [int(hiq.rstrip()), int(lowq.rstrip())]
    return qc_dict


def shell_command(command_string):
    logger.info(command_string)
    try:
        return subprocess.check_output(shlex.split(command_string))
    except subprocess.CalledProcessError as e:
        logger.error(
            "%s exited with return code %s" % (command_string, e.returncode))
        logger.error(e.output)
        return None
    except:
        raise

def scrub_fun(in_filepath, out_filepath):
    # Check the input.
    logger.debug("Input flagstat for %s" % (in_filepath))
    logger.debug(shell_command("samtools flagstat %s" % (in_filepath)))
    # Set up the paths to inputs and outputs.
    dirname = os.path.dirname(out_filepath)
    header_path = os.path.join(dirname, "header.txt")
    sam_path = os.path.join(dirname, "scrubbed.sam")
    # Cache the header.
    shell_command("samtools view -H %s -o %s" % (in_filepath, header_path))
    # Scrub the sequence information from these fields:
    # 6 = CIGAR, 10 = query sequence, 11 = PHRED, and suppress optional tags
    # For example, unscrubbed read might look like:
    # SPADE:8:33:220:1107#0 0 chr21 8994907 37 9M1I26M * 0 0 ATTGTTGACAAAAACTCGACAAACAATTGGAGAATC    bbbR]`T`^]TTSSS^_W`BBBBBBBBBBBBBBBBB    X0:i:1  X1:i:0  MD:Z:35 PG:Z:MarkDuplicates     XG:i:1  NM:i:1  XM:i:0  XO:i:1  XT:A:U
    # Scrubbed version would look like:
    # SPADE:8:33:220:1107#0 0 chr21 8994907 37 36M     * 0 0 NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN    *
    common.run_pipe([
        'samtools view %s' % (in_filepath),
        r"""awk '{OFS="\t"} {s=""; for(i=1;i<=length($10);i++) s=(s "N"); $6=(i-1 "M"); $10=s; $11="*"; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}'"""
        ], sam_path)
    # Add back the header.
    common.run_pipe([
        'cat %s %s' % (header_path, sam_path),
        'samtools view -S -b - -o %s' % (out_filepath)])
    # Check the output.
    logger.debug("Output flagstat for %s" % (out_filepath))
    logger.debug(shell_command("samtools flagstat %s" % (out_filepath)))

def scrub_main(input_bams):
    # Initialize data object inputs on the platform
    # into dxpy.DXDataObject instances.
    input_bams = [item for item in input_bams]
    # Download each file input to a new directory in the the local file system
    # using variable names for the filenames.
    # Construct output filenames.
    # Dispatch jobs to a pool of workers.
    out_paths = []
    pool = Pool()  # default is pool of cpu_count() workers
    for i, bam in enumerate(input_bams):
        dirname = str(i)
        filename = bam.split("/")[-1]
        os.mkdir(dirname)
        #in_path = os.path.join(dirname, filename)
        in_path = os.path.join(bam)
        out_path = os.path.join(dirname, "scrub-" + filename)
        out_paths.append(out_path)
        pool.apply_async(scrub_fun, (in_path, out_path))
    # Close the worker pool and block until all jobs are complete.
    pool.close()
    pool.join()
    # Populate output fields and return.
    scrubbed_bams = [path for path in out_paths]
    output = {
        "scrubbed_bams": [output_bam for output_bam in scrubbed_bams]
    }
    return output

def main(input_bam, paired_end, samtools_version, samtools_params, picard_version, scrub, debug):
    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    samtools = SAMTOOLS_PATH.get(samtools_version)
    assert samtools, "samtools version %s is not supported" % (samtools_version)
    picard = PICARD_PATH.get(picard_version)
    assert picard, "picard version %s is not supported" % (picard_version)
    logger.info("In postprocess with samtools %s and picard %s" % (samtools, picard))

    raw_bam_file = input_bam
    raw_bam_filename = raw_bam_file
    raw_bam_basename = raw_bam_file.rstrip('.bam')
    raw_bam_file_mapstats_filename = raw_bam_basename + '.flagstat.qc'
    subprocess.check_output('set -x; ls -l', shell=True)

    # Generate initial mapping statistics
    with open(raw_bam_file_mapstats_filename, 'w') as fh:
        flagstat_command = "%s flagstat %s" % (samtools, raw_bam_filename)
        logger.info(flagstat_command)
        subprocess.check_call(shlex.split(flagstat_command), stdout=fh)

    filt_bam_prefix = raw_bam_basename + ".filt.srt"
    filt_bam_filename = filt_bam_prefix + ".bam"
    if paired_end:
        # =============================
        # Remove  unmapped, mate unmapped
        # not primary alignment, reads failing platform
        # Remove low MAPQ reads
        # Only keep properly paired reads
        # Obtain name sorted BAM file
        # ==================
        tmp_filt_bam_prefix = "%s.tmp" % (filt_bam_prefix)  # was tmp.prefix.nmsrt
        tmp_filt_bam_filename = tmp_filt_bam_prefix + ".bam"
        out, err = common.run_pipe([
            # filter: -F 1804 FlAG bits to exclude; -f 2 FLAG bits to reqire;
            # -q 30 exclude MAPQ < 30; -u uncompressed output
            # exclude FLAG 1804: unmapped, next segment unmapped, secondary
            # alignments, not passing platform q, PCR or optical duplicates
            # require FLAG 2: properly aligned
            "%s view -F 1804 -f 2 %s -u %s" % (samtools, samtools_params, raw_bam_filename),
            # sort:  -n sort by name; - take input from stdin;
            # out to specified filename
            # Will produce name sorted BAM
            "%s sort -n - %s" % (samtools, tmp_filt_bam_prefix)])
        if err:
            logger.error("samtools error: %s" % (err))
        # Remove orphan reads (pair was removed)
        # and read pairs mapping to different chromosomes
        # Obtain position sorted BAM
        subprocess.check_output('set -x; ls -l', shell=True)
        out, err = common.run_pipe([
            # fill in mate coordinates, ISIZE and mate-related flags
            # fixmate requires name-sorted alignment; -r removes secondary and
            # unmapped (redundant here because already done above?)
            # - send output to stdout
            "%s fixmate -r %s -" % (samtools, tmp_filt_bam_filename),
            # repeat filtering after mate repair
            "%s view -F 1804 -f 2 -u -" % (samtools),
            # produce the coordinate-sorted BAM
            "%s sort - %s" % (samtools, filt_bam_prefix)])
        subprocess.check_output('set -x; ls -l', shell=True)
    else:  # single-end data
        # =============================
        # Remove unmapped, mate unmapped
        # not primary alignment, reads failing platform
        # Remove low MAPQ reads
        # Obtain name sorted BAM file
        # ==================
        with open(filt_bam_filename, 'w') as fh:
            samtools_filter_command = (
                "%s view -F 1804 %s -b %s"
                % (samtools, samtools_params, raw_bam_filename)
                )
            logger.info(samtools_filter_command)
            subprocess.check_call(
                shlex.split(samtools_filter_command),
                stdout=fh)

    # ========================
    # Mark duplicates
    # ======================
    tmp_filt_bam_filename = raw_bam_basename + ".dupmark.bam"
    dup_file_qc_filename = raw_bam_basename + ".dup.qc"
    picard_string = ' '.join([
        "java -Xmx4G -jar %s" % (picard),
        "INPUT=%s" % (filt_bam_filename),
        "OUTPUT=%s" % (tmp_filt_bam_filename),
        "METRICS_FILE=%s" % (dup_file_qc_filename),
        "VALIDATION_STRINGENCY=LENIENT",
        "ASSUME_SORTED=true",
        "REMOVE_DUPLICATES=false"
        ])
    logger.info(picard_string)
    subprocess.check_output(shlex.split(picard_string))
    os.rename(tmp_filt_bam_filename, filt_bam_filename)

    if paired_end:
        final_bam_prefix = raw_bam_basename + ".filt.srt.nodup"
    else:
        final_bam_prefix = raw_bam_basename + ".filt.nodup.srt"
    final_bam_filename = final_bam_prefix + ".bam"  # To be stored
    final_bam_index_filename = final_bam_filename + ".bai"  # To be stored
    # QC file
    final_bam_file_mapstats_filename = final_bam_prefix + ".flagstat.qc"

    if paired_end:
        samtools_dedupe_command = \
            "%s view -F 1804 -f2 -b %s" % (samtools, filt_bam_filename)
    else:
        samtools_dedupe_command = \
            "%s view -F 1804 -b %s" % (samtools, filt_bam_filename)

    # ============================
    # Remove duplicates
    # Index final position sorted BAM
    # ============================
    with open(final_bam_filename, 'w') as fh:
        logger.info(samtools_dedupe_command)
        subprocess.check_call(
            shlex.split(samtools_dedupe_command),
            stdout=fh)
    # Index final bam file
    samtools2=SAMTOOLS_PATH.get("1.3.1")
    samtools_index_command = \
        "%s index %s %s" % (samtools2, final_bam_filename, final_bam_index_filename)
    logger.info(samtools_index_command)
    subprocess.check_output(shlex.split(samtools_index_command))

    # Generate mapping statistics
    with open(final_bam_file_mapstats_filename, 'w') as fh:
        flagstat_command = "%s flagstat %s" % (samtools, final_bam_filename)
        logger.info(flagstat_command)
        subprocess.check_call(shlex.split(flagstat_command), stdout=fh)

    # =============================
    # Compute library complexity
    # =============================
    # Sort by name
    # convert to bedPE and obtain fragment coordinates
    # sort by position and strand
    # Obtain unique count statistics
    pbc_file_qc_filename = final_bam_prefix + ".pbc.qc"
    # PBC File output
    # TotalReadPairs [tab]
    # DistinctReadPairs [tab]
    # OneReadPair [tab]
    # TwoReadPairs [tab]
    # NRF=Distinct/Total [tab]
    # PBC1=OnePair/Distinct [tab]
    # PBC2=OnePair/TwoPair
    if paired_end:
        steps = [
            "%s sort -no %s -" % (samtools, filt_bam_filename),
            "bamToBed -bedpe -i stdin",
            r"""awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}'"""]
    else:
        steps = [
            "bamToBed -i %s" % (filt_bam_filename),
            r"""awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}'"""]
    steps.extend([
        "grep -v 'chrM'",
        "sort",
        "uniq -c",
        r"""awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{if(m2){printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}else{printf "%d\t%d\t%d\t%d\t%f\t%f\t%s\n",mt,m0,m1,m2,m0/mt,m1/m0,"Inf"}}'"""
        ])
    out, err = common.run_pipe(steps, pbc_file_qc_filename)
    if err:
        logger.error("PBC file error: %s" % (err))

    output = {}
    logger.info("Uploading results files to the project")
    filtered_bam = final_bam_filename
    filtered_bam_index = final_bam_index_filename
    output.update({
        "filtered_bam": filtered_bam,
        "filtered_bam_index": filtered_bam_index
    })

    # If the scrub parameter is true, pass the bams to the scrub applet.
    if scrub:
        scrub_subjob = scrub_main([input_bam, filtered_bam])
        scrubbed_unfiltered_bam = scrub_subjob.get("scrubbed_bams")[0]
        scrubbed_filtered_bam = scrub_subjob.get("scrubbed_bams")[1]
        # Add the optional scrubbed outputs.
        output.update({
            "scrubbed_unfiltered_bam": scrubbed_unfiltered_bam,
            "scrubbed_filtered_bam": scrubbed_filtered_bam
        })

    # Upload or calculate the remaining outputs.
    filtered_mapstats = final_bam_file_mapstats_filename
    dup_file = dup_file_qc_filename
    pbc_file = pbc_file_qc_filename

    logger.info("Calcualting QC metrics")
    dup_qc = dup_parse(dup_file_qc_filename)
    pbc_qc = pbc_parse(pbc_file_qc_filename)
    initial_mapstats_qc = flagstat_parse(raw_bam_file_mapstats_filename)
    final_mapstats_qc = flagstat_parse(final_bam_file_mapstats_filename)
    if paired_end:
        useable_fragments = final_mapstats_qc.get('in_total')[0]/2
    else:
        useable_fragments = final_mapstats_qc.get('in_total')[0]
    logger.info("initial_mapstats_qc: %s" % (initial_mapstats_qc)),
    logger.info("final_mapstats_qc: %s" % (final_mapstats_qc)),
    logger.info("dup_qc: %s" % (dup_qc))
    logger.info("pbc_qc: %s" % (pbc_qc))

    # Return links to the output files and values.
    output.update({
        "filtered_mapstats": filtered_mapstats,
        "dup_file_qc": dup_file,
        "pbc_file_qc": pbc_file,
        "paired_end": paired_end,
        "n_reads_input": str(initial_mapstats_qc.get('in_total')[0]),
        "picard_read_pairs_examined": str(dup_qc.get('read_pairs_examined')),
        "picard_unpaired_reads_examined": str(dup_qc.get('unpaired_reads_examined')),
        "picard_read_pair_duplicates": str(dup_qc.get('read_pair_duplicates')),
        "picard_unpaired_read_duplicates": str(dup_qc.get('unpaired_read_duplicates')),
        "useable_fragments": str(useable_fragments),
        "NRF": str(pbc_qc.get('NRF')),
        "PBC1": str(pbc_qc.get('PBC1')),
        "PBC2": str(pbc_qc.get('PBC2')),
        "duplicate_fraction": str(dup_qc.get('percent_duplication'))
    })
    parse_file = final_bam_prefix + ".parse"
    with open(parse_file,"w") as fh:
        for key, val in output.items():
            if isinstance(val, list):
                fh.write(": ".join([key, ", ".join(val)])+"\n")
            else:
                fh.write(": ".join([key, str(val)])+"\n")
    logger.info("Exiting with output:\n%s" % (pformat(output)))
    return output


'''

def main(input_bam, paired_end, samtools_version, samtools_params, picard_version, scrub, debug):
    if paired_end:
        print("paired_end==True")
    else:
        print("paired_end==False")
    if scrub:
        print("scrub==True")
    else:
        print("scrub==False")
    if debug:
        print("debug==True")
    else:
        print("debug==False")
'''
if __name__ == '__main__':
    main(input_bam, paired_end, samtools_version, samtools_params, picard_version, scrub, debug)
