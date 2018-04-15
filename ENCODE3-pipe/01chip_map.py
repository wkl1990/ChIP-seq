#!/usr/bin/env python3
#Description: ENCODE_map 0.0.1 modified by WKL
#Author: WKL
#Environment: Python3.6
#Date: 2018/04/03
#Examples: python3 chip_map.py --in1="R1.fq" --in2="R2.fq" --ref="/zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta" \
#    --crop_len=native --bwa="0.7.10" --bwa_params="-q 5 -l 32 -k 2" --samtools="1.2" --debug=False

import os
import re
import getopt
import subprocess
import shlex
from multiprocessing import cpu_count
#import dxpy
import logging
import sys
sys.path.insert(1,"/zs32_2/klwang/ChIP_seq/test_data/ENCODE_code")
import common
sys.path.remove("/zs32_2/klwang/ChIP_seq/test_data/ENCODE_code")


logger = logging.getLogger(__name__)
#logger.addHandler(dxpy.DXLogHandler())
logger.propagate = False
logger.setLevel(logging.INFO)

BWA_PATH = {
    "0.7.10": "/opt/tools/seq-analysis/speedseq/bin/bwa",
    "0.7.12": "/opt/tools/seq-analysis/bwa-0.7.12/bwa"
}

SAMTOOLS_PATH = {
    "1.3.1": "/zs32/home/klwang/softwares/samtools-1.3.1/bin/samtools",
    "1.2": "/opt/tools/seq-analysis/samtools-1.2/samtools"
}

TRIMMOMATIC_PATH = "/opt/tools/seq-analysis/Trimmomatic-0.36/trimmomatic-0.36.jar"

# the order of this list is important.
# strip_extensions strips from the right inward, so
# the expected right-most extensions should appear first (like .gz)
STRIP_EXTENSIONS = ['.gz', '.fq', '.fastq', '.fa', '.fasta']

def usage():
    print("")
    print("Align fastq reads using bwa.")
    print("usage: python %s -option <argument>" %sys.argv[0])
    print("  -h/--help                   ")
    print("  --in1=<STRING>     [required] FASTA/FASTQ input file.")
    print("  --in2=<STRING>                FASTQ (PE, R2) input file.")
    print("  --ref=<STRING>                The reference file.")
    print("  --crop_len=<True or native>   True will crop fastq reads.")
    print("  --bwa=<STRING>                The bwa version used.")
    print("  --bwa_params=<STRING>         The bwa parameters used.")
    print("  --samtools=<STRING>           The samtools version used.")
    print("  --debug=<True or False>       Whether debug or not.")

## deal with options
try:
    opts, args = getopt.getopt( sys.argv[1:], "h", ["help", "in1=", "in2=", "ref=", "crop_len=", "bwa=", "bwa_params=", "samtools=", "debug=" ] )
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
            reads1 = val
        if opt in ( "--in2", ):
            reads2 = val
        if opt in ( "--ref", ):
            reference = val
        if opt in ( "--crop_len", ):
            crop_length = val
        if opt in ( "--bwa", ):
            bwa_version = val
        if opt in ( "--bwa_params", ):
            bwa_aln_params = val
        if opt in ( "--samtools", ):
            samtools_version = val        
        if opt in ( "--debug", ):
            debug = val
if not "reads2" in dir():
    reads2 = None 

#print(reads2)
'''
try:
    debug = False if not debug else debug
except:
    debug = False
'''
if not "debug" in dir() or not debug or debug=="False":
    debug=False
else:
    debug=True

def strip_extensions(filename, extensions):
    basename = filename
    for extension in extensions:
        basename = basename.rpartition(extension)[0] or basename
    return basename


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

def crop(reads1_file, reads2_file, crop_length, debug):
    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    logger.setLevel(logging.INFO)
    if crop_length == 'native':
        output = dict(zip(
            ["cropped_reads1", "cropped_reads2"], [reads1_file, reads2_file]))
    else:
        reads1_filename = reads1_file
        reads1_basename = strip_extensions(reads1_filename, STRIP_EXTENSIONS)
        if reads2_file:
            end_string = "PE"
            reads2_filename = reads2_file
            reads2_basename = \
                strip_extensions(reads2_filename, STRIP_EXTENSIONS)
            output_fwd_paired_filename = reads1_basename + '-crop-paired.fq.gz'
            output_fwd_unpaired_filename = \
                reads1_basename + '-crop-unpaired.fq.gz'
            output_rev_paired_filename = reads2_basename + '-crop-paired.fq.gz'
            output_rev_unpaired_filename = \
                reads2_basename + '-crop-unpaired.fq.gz'
            SE_output_filename = None
            adapter = "/".join([TRIMMOMATIC_PATH[0:TRIMMOMATIC_PATH.rfind("/")],"adapters/TruSeq3-PE.fa:2:30:10"])
        else:
            end_string = "SE"
            reads2_filename = None
            reads2_basename = None
            output_fwd_paired_filename = None
            output_fwd_unpaired_filename = None
            output_rev_paired_filename = None
            output_rev_unpaired_filename = None
            SE_output_filename = reads1_basename + "-crop.fq.gz"
            adapter = "/".join([TRIMMOMATIC_PATH[0:TRIMMOMATIC_PATH.rfind("/")],"adapters/TruSeq3-SE.fa:2:30:10"])
        crop_command = ' '.join([s for s in [
            'java -jar',
            TRIMMOMATIC_PATH,
            end_string,
            '-threads %d' % (cpu_count()),
            reads1_filename,
            reads2_filename,
            SE_output_filename,
            output_fwd_paired_filename,
            output_fwd_unpaired_filename,
            output_rev_paired_filename,
            output_rev_unpaired_filename,
            'ILLUMINACLIP:%s' % (adapter), 
            'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35',]
            #'MINLEN:%s' % (crop_length),
            #'CROP:%s' % (crop_length)]
            if s])
        logger.info("Cropping with: %s" % (crop_command))
        print(subprocess.check_output(shlex.split(crop_command)))
        print(subprocess.check_output(shlex.split('ls -l')))
        if SE_output_filename:
            cropped_reads = [SE_output_filename, None]
        else:
            cropped_reads = [output_fwd_paired_filename,
                output_rev_paired_filename]
        output = dict(zip(["cropped_reads1", "cropped_reads2"], cropped_reads))
    logger.info("returning from crop with output %s" % (output))
    return output


def postprocess(indexed_reads, unmapped_reads, reference,
                bwa_version, samtools_version, debug):
    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    samtools = SAMTOOLS_PATH.get(samtools_version)
    assert samtools, "samtools version %s is not supported" % (samtools_version)
    bwa = BWA_PATH.get(bwa_version)
    assert bwa, "BWA version %s is not supported" % (bwa_version)
    logger.info("In postprocess with samtools %s and bwa %s" % (samtools, bwa))
    indexed_reads_filenames = []
    unmapped_reads_filenames = []
    for i, reads in enumerate(indexed_reads):
        read_pair_number = i+1
        fn = reads
        logger.info("indexed_reads %d: %s" % (read_pair_number, fn))
        indexed_reads_filenames.append(fn)
        unmapped = unmapped_reads[i]
        fn = unmapped
        logger.info("unmapped reads %d: %s" % (read_pair_number, fn))
        unmapped_reads_filenames.append(fn)
    reference_filename = reference
    logger.info("reference: %s" % (reference_filename))
    # extract the reference files from the tar
    reference_filename = reference_filename
    logger.info("Using reference file: %s" % (reference_filename))
    paired_end = len(indexed_reads) == 2
    if paired_end:
        r1_basename = strip_extensions(
            unmapped_reads_filenames[0], STRIP_EXTENSIONS)
        r2_basename = strip_extensions(
            unmapped_reads_filenames[1], STRIP_EXTENSIONS)
        if r1_basename.rfind("/") == -1:
            reads_basename = r1_basename + r2_basename
        else:
            reads_basename = "/".join([r1_basename[0:r1_basename.rfind("/")], (r1_basename.split('/')[-1] + r2_basename.split('/')[-1])])
    else:
        reads_basename = strip_extensions(
            unmapped_reads_filenames[0], STRIP_EXTENSIONS)
    raw_bam_filename = '%s.raw.srt.bam' % (reads_basename)
    raw_bam_mapstats_filename = '%s.raw.srt.bam.flagstat.qc' % (reads_basename)
    if paired_end:
        reads1_filename = indexed_reads_filenames[0]
        reads2_filename = indexed_reads_filenames[1]
        unmapped_reads1_filename = unmapped_reads_filenames[0]
        unmapped_reads2_filename = unmapped_reads_filenames[1]
        raw_sam_filename = reads_basename + ".raw.sam"
        badcigar_filename = "badreads.tmp"
        steps = [
            "%s sampe -P %s %s %s %s %s"
            % (bwa, reference_filename, reads1_filename, reads2_filename,
               unmapped_reads1_filename, unmapped_reads2_filename),
            "tee %s" % (raw_sam_filename),
            r"""awk 'BEGIN {FS="\t" ; OFS="\t"} ! /^@/ && $6!="*" { cigar=$6; gsub("[0-9]+D","",cigar); n = split(cigar,vals,"[A-Z]"); s = 0; for (i=1;i<=n;i++) s=s+vals[i]; seqlen=length($10) ; if (s!=seqlen) print $1"\t" ; }'""",
            "sort",
            "uniq"]
        out, err = common.run_pipe(steps, badcigar_filename)
        print(out)
        if err:
            logger.error("sampe error: %s" % (err))
        steps = [
            "cat %s" % (raw_sam_filename),
            "grep -v -F -f %s" % (badcigar_filename)]
    else:  # single end
        reads_filename = indexed_reads_filenames[0]
        unmapped_reads_filename = unmapped_reads_filenames[0]
        steps = [
            "%s samse %s %s %s"
            % (bwa, reference_filename,
               reads_filename, unmapped_reads_filename)]
    if samtools_version == "0.1.9":
        steps.extend([
            "%s view -Su -" % (samtools),
            "%s sort - %s"
            % (samtools, raw_bam_filename.rstrip('.bam'))])  # samtools adds .bam
    else:
        steps.extend([
            "%s view -@%d -Su -" % (samtools, cpu_count()),
            "%s sort -@%d - %s"
            % (samtools, cpu_count(), raw_bam_filename.rstrip('.bam'))])  # samtools adds .bam
    logger.info("Running pipe: %s" % (steps))
    out, err = common.run_pipe(steps)
    if out:
        print(out)
    if err:
        logger.error("samtools error: %s" % (err))
    with open(raw_bam_mapstats_filename, 'w') as fh:
        subprocess.check_call(
            shlex.split("%s flagstat %s" % (samtools, raw_bam_filename)),
            stdout=fh)
    print(subprocess.check_output('ls -l', shell=True))
    mapped_reads = raw_bam_filename
    mapping_statistics = raw_bam_mapstats_filename
    flagstat_qc = flagstat_parse(raw_bam_mapstats_filename)
    output = {
        'mapped_reads': mapped_reads,
        'mapping_statistics': mapping_statistics,
        'n_mapped_reads': flagstat_qc.get('mapped')[0]  # 0 is hi-q reads
    }
    logger.info("Returning from postprocess with output: %s" % (output))
    return output


def process(reads_file, reference, bwa_aln_params, bwa_version, debug):
    # reads_file, reference should be links to file objects.
    # reference should be a tar of files generated by bwa index and
    # the tar should be uncompressed to avoid repeating the decompression.
    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    bwa = BWA_PATH.get(bwa_version)
    assert bwa, "BWA version %s is not supported" % (bwa_version)
    logger.info("In process with bwa %s" % (bwa))
    # Generate filename strings and download the files to the local filesystem
    reads_filename = reads_file
    reads_basename = strip_extensions(reads_filename, STRIP_EXTENSIONS)
    reference_filename = reference
    reference_filename = reference_filename
    logger.info("Using reference file: %s" % (reference_filename))
    print(subprocess.check_output('ls -l', shell=True))
    # generate the suffix array index file
    sai_filename = '%s.sai' % (reads_basename)
    with open(sai_filename, 'w') as sai_file:
        # Build the bwa command and call bwa
        bwa_command = "%s aln %s -t %d %s %s" \
            % (bwa, bwa_aln_params, cpu_count(),
               reference_filename, reads_filename)
        logger.info("Running bwa with %s" % (bwa_command))
        subprocess.check_call(shlex.split(bwa_command), stdout=sai_file)
    print(subprocess.check_output('ls -l', shell=True))
    # Upload the output to the DNAnexus project
    logger.info("Uploading suffix array %s" % (sai_filename))
    output = {"suffix_array_index": sai_filename}
    logger.info("Returning from process with %s" % (output))
    return output


def main(reads1, crop_length, reference,
         bwa_version, bwa_aln_params, samtools_version, debug, reads2=None):
    # Main entry-point.  Parameter defaults assumed to come from dxapp.json.
    # reads1, reference, reads2 are links to DNAnexus files or None
    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    ## This spawns only one or two subjobs for single- or paired-end,
    # respectively.  It could also download the files, chunk the reads,
    # and spawn multiple subjobs.
    ## Files are downloaded later by subjobs into their own filesystems
    # and uploaded to the project.
    ## Initialize file handlers for input files.
    paired_end = reads2 is not None
    if crop_length == 'native':
        crop_subjob = None
        unmapped_reads = [reads1, reads2]
    else:
        crop_subjob_input = {
            "reads1_file": reads1,
            "reads2_file": reads2,
            "crop_length": crop_length,
            "debug": debug
        }
        logger.info("Crop job input: %s" % (crop_subjob_input))
        crop_subjob = crop(crop_subjob_input["reads1_file"],crop_subjob_input["reads2_file"],crop_subjob_input["crop_length"],crop_subjob_input["debug"])
        unmapped_reads = [crop_subjob.get("cropped_reads1")]
        if paired_end:
            unmapped_reads.append(crop_subjob.get("cropped_reads2"))
        else:
            unmapped_reads.append(None)

    unmapped_reads = [r for r in unmapped_reads if r]

    mapping_subjobs = []
    for reads in unmapped_reads:
        mapping_subjob_input = {
            "reads_file": reads,
            "reference": reference,
            "bwa_aln_params": bwa_aln_params,
            "bwa_version": bwa_version,
            "debug": debug
        }
        logger.info("Mapping job input: %s" % (mapping_subjob_input))
        if crop_subjob:
            mapping_subjobs.append(
                process(reads,reference,bwa_aln_params,bwa_version,debug)
                )
        else:
            mapping_subjobs.append(
                process(reads,reference,bwa_aln_params,bwa_version,debug)
                )

    # Create the job that will perform the "postprocess" step.
    # depends_on=mapping_subjobs, so blocks on all mapping subjobs
    fn_input={
        "indexed_reads": [
            subjob.get("suffix_array_index")
            for subjob in mapping_subjobs],
        "unmapped_reads": unmapped_reads,
        "reference": reference,
        "bwa_version": bwa_version,
        "samtools_version": samtools_version,
        "debug": debug}
    postprocess_job = postprocess(fn_input["indexed_reads"], fn_input["unmapped_reads"], fn_input["reference"],
                fn_input["bwa_version"], fn_input["samtools_version"], fn_input["debug"])

    mapped_reads = postprocess_job.get("mapped_reads")
    mapping_statistics = postprocess_job.get("mapping_statistics")
    n_mapped_reads = postprocess_job.get("n_mapped_reads")

    output = {
        "mapped_reads": mapped_reads,
        "crop_length": crop_length,
        "mapping_statistics": mapping_statistics,
        "paired_end": paired_end,
        "n_mapped_reads": n_mapped_reads
    }
    reads_filename = reads1
    reads_basename = strip_extensions(reads_filename, STRIP_EXTENSIONS)
    map_parse = '%s.map.parse' % (reads_basename)
    with open(map_parse,"w") as fh:
        for key, val in output.items():
            if isinstance(val, list):
                fh.write(": ".join([key, ", ".join(val)])+"\n")
            else:
                fh.write(": ".join([key, str(val)])+"\n")
    logger.info("Exiting with output: %s" % (output))
    return output
'''
def main(reads1, crop_length, reference,
         bwa_version, bwa_aln_params, samtools_version, debug, reads2=None):
    # Main entry-point.  Parameter defaults assumed to come from dxapp.json.
    # reads1, reference, reads2 are links to DNAnexus files or None
    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    ## This spawns only one or two subjobs for single- or paired-end,
    # respectively.  It could also download the files, chunk the reads,
    # and spawn multiple subjobs.
    ## Files are downloaded later by subjobs into their own filesystems
    # and uploaded to the project.
    ## Initialize file handlers for input files.
    paired_end = reads2 is not None
    if paired_end:
        print("paired_end==True")
    else:
        print("paired_end==False")
    if debug:
        print("debug == True")
    else: 
        print("debug == False")
    if crop_length=="native":
        print("crop==native")
    else:
        print("crop!=native")
'''
#dxpy.run()
if __name__ == '__main__':
    main(reads1, crop_length, reference, bwa_version, bwa_aln_params, samtools_version, debug, reads2)
