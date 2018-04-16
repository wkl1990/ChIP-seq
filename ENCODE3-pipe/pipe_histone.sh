##Description: pipeline for histone ChIP-seq with unary control and unreplicated
##Author: WKL
##Environment: Shell
##Date: 2018/04/12

#set path
cd ChIP-seq
#set the args
reference="~/human_UCSC_hg19/chrAll.fasta"
chromsizes="~/male.hg19.chrom.sizes"
narrowpeakAs="~/narrowPeak.as" 
gappedpeakAs="~/gappedPeak.as" 
broadpeakAs="~/broadPeak.as"

exp_fq="data/rep1-ENCFF921SED.chr21.fq.gz"
input_fq="data/ctl1-ENCFF690VPV.chr21.fq.gz"

dirname=${exp_fq%/*}
exp_filename=${exp_fq##*/}
input_filename=${input_fq##*/}
exp_basename=`echo $exp_filename | sed -e s/.fq// -e s/.gz// -e s/.fastq// -e s/.fa// -e s/.fasta//`
input_basename=`echo $input_filename | sed -e s/.fq// -e s/.gz// -e s/.fastq// -e s/.fa// -e s/.fasta//`
peak_dirname=${exp_basename%.*}_${input_basename%.*}_peaks_macs
echo $dirname $exp_basename $input_basename

#fastqc
# ls data/*.fq.gz | while read id; do file=$(basename $id); sample=${file%.*} fastqclogfile="logfile/"${sample}".fastqc.raw.log"; fastqc --outdir fastqc/fastq_raw $file > $fastqclogfile; done

#rep1 trim and map 
python3 ENCODE3-pipe/01chip_map.py --in1=$exp_fq --ref=$reference --crop_len=True \
    --bwa="0.7.10" --bwa_params="-q 5 -l 32 -k 2" --samtools="1.2" --debug=False 1>${exp_basename}.map.log1 2>${exp_basename}.map.log2 &
#ctl1 trim and map
python3 ENCODE3-pipe/01chip_map.py --in1=$input_fq --ref=$reference --crop_len=True \
    --bwa="0.7.10" --bwa_params="-q 5 -l 32 -k 2" --samtools="1.2" --debug=False 1>${input_basename}.map.log1 2>${input_basename}.map.log2 

#rep1 filter qc (1,remove unmapped, low MAPQ reads, 2,remove duplicates, 3,compute library complexity)
rep_bam=`grep "^mapped_reads: " ${dirname}"/"${exp_basename}".map.parse" | sed 's/^mapped_reads: //'`
python3 ENCODE3-pipe/02chip_filterQC.py --in=$rep_bam --pair=False --scrub=False --picard="1.119" \
    --samtools_params="-q 30" --samtools="1.2" --debug=False 1>${exp_basename}.filterqc.log1 2>${exp_basename}.filterqc.log2 &
#ctl1 filter qc
input_bam=`grep "^mapped_reads: " ${dirname}"/"${input_basename}".map.parse" | sed 's/^mapped_reads: //'`
python3 ENCODE3-pipe/02chip_filterQC.py  --in=$input_bam --pair=False --scrub=False --picard="1.119" \
    --samtools_params="-q 30" --samtools="1.2" --debug=False 1>${input_basename}.filterqc.log1 2>${input_basename}.filterqc.log2 
unset  rep_bam input_bam

#rep1 Xcor (1,Create tagAlign file, 2,Subsample tagAlign file, 3,Calculate Cross-correlation QC scores)
rep_filtbam=`grep "^filtered_bam: " ${dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.parse | sed 's/^filtered_bam: //'`
python3 ENCODE3-pipe/03chip_xcor.py --in1=$rep_filtbam --pair=False \
    --spp="1.14" 1>${exp_basename}.xcor.log1 2>${exp_basename}.xcor.log2 &
#ctl1 Xcor
input_filtbam=`grep "^filtered_bam: " ${dirname}/${input_basename}-crop.raw.srt.filt.nodup.srt.parse | sed 's/^filtered_bam: //'`
python3 ENCODE3-pipe/03chip_xcor.py --in1=$input_filtbam --pair=False \
    --spp="1.14" 1>${input_basename}.xcor.log1 2>${input_basename}.xcor.log2 
unset rep_filtbam input_filtbam

#macs peak (1,pseudoreplicator, 2, macs for real samples and pseudoreplicator samples and pooled samples)
declare rep_tagAlign=`grep "^tagAlign_file: " ${dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.parse | sed 's/^tagAlign_file: //'`
declare input_tagAlign=`grep "^tagAlign_file: " ${dirname}/${input_basename}-crop.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.parse | sed 's/^tagAlign_file: //'`
declare xcor1=`grep "CC_scores_file: " ${dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.parse | sed 's/CC_scores_file: //'`
python3 ENCODE3-pipe/04chip_macs2.py --in1=$rep_tagAlign --ct1=$input_tagAlign --xcor1=$xcor1 --in1_PE=False \
    --chromsizes=$chromsizes --genomesize="hs" --narrowpeak=$narrowpeakAs --gappedpeak=$gappedpeakAs \
    --broadpeak=$broadpeakAs 1>${exp_basename}.macs.log1 2>${exp_basename}.macs.log2 
unset input_tagAlign

#final narrowpeak (1,find overlapped peaks, 2,compute the FRiP)
#declare rep_tagAlign=`grep "^tagAlign_file: " ${dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.parse | sed 's/^tagAlign_file: //'`
#declare xcor1=`grep "CC_scores_file: " ${dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.filt.nodup.sample.15.SE.tagAlign.gz.parse | sed 's/CC_scores_file: //'`

declare peak1=`grep "rep1pr1_narrowpeaks: " ${peak_dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.SE.tagAlign.gz.peaksfile | sed 's/rep1pr1_narrowpeaks: //'`
declare peak2=`grep "rep1pr2_narrowpeaks: " ${peak_dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.SE.tagAlign.gz.peaksfile | sed 's/rep1pr2_narrowpeaks: //'`
declare peakpool=`grep "rep1_narrowpeaks: " ${peak_dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.SE.tagAlign.gz.peaksfile | sed 's/rep1_narrowpeaks: //'`
declare signal1=`grep "rep1_fc_signal: " ${peak_dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.SE.tagAlign.gz.peaksfile | sed 's/rep1_fc_signal: //'`
ll $peak1 $peak2 $peakpool $signal1

python3 ENCODE3-pipe/05chip_overlap_peaks.py --in1=$rep_tagAlign --peak1=$peak1 --peak2=$peak2 --peakpool=$peakpool \
    --xcor1=$xcor1 --pair=False --chromsizes=$chromsizes --peaktype="narrowPeak" --peakas=$narrowpeakAs \
    --signal1=$signal1 --prefix="${peak_dirname}/final" 1>${exp_basename}.final.narrowpeak.log1 2>${exp_basename}.final.narrowpeak.log2 
unset peak1 peak2 peakpool signal1

#final broadpeak
declare peak1=`grep "rep1pr1_broadpeaks: " ${peak_dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.SE.tagAlign.gz.peaksfile | sed 's/rep1pr1_broadpeaks: //'`
declare peak2=`grep "rep1pr2_broadpeaks: " ${peak_dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.SE.tagAlign.gz.peaksfile | sed 's/rep1pr2_broadpeaks: //'`
declare peakpool=`grep "rep1_broadpeaks: " ${peak_dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.SE.tagAlign.gz.peaksfile | sed 's/rep1_broadpeaks: //'`
declare signal1=`grep "rep1_fc_signal: " ${peak_dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.SE.tagAlign.gz.peaksfile | sed 's/rep1_fc_signal: //'`
ll $peak1 $peak2 $peakpool $signal1

python3 ENCODE3-pipe/05chip_overlap_peaks.py --in1=$rep_tagAlign --peak1=$peak1 --peak2=$peak2 --peakpool=$peakpool \
    --xcor1=$xcor1 --pair=False --chromsizes=$chromsizes --peaktype="broadPeak" --peakas=$broadpeakAs \
    --signal1=$signal1 --prefix="${peak_dirname}/final" 1>${exp_basename}.final.broadpeak.log1 2>${exp_basename}.final.broadpeak.log2 
unset peak1 peak2 peakpool signal1

#final gappedpeak
declare peak1=`grep "rep1pr1_gappedpeaks: " ${peak_dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.SE.tagAlign.gz.peaksfile | sed 's/rep1pr1_gappedpeaks: //'`
declare peak2=`grep "rep1pr2_gappedpeaks: " ${peak_dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.SE.tagAlign.gz.peaksfile | sed 's/rep1pr2_gappedpeaks: //'`
declare peakpool=`grep "rep1_gappedpeaks: " ${peak_dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.SE.tagAlign.gz.peaksfile | sed 's/rep1_gappedpeaks: //'`
declare signal1=`grep "rep1_fc_signal: " ${peak_dirname}/${exp_basename}-crop.raw.srt.filt.nodup.srt.SE.tagAlign.gz.peaksfile | sed 's/rep1_fc_signal: //'`
ll $peak1 $peak2 $peakpool $signal1

python3 ENCODE3-pipe/05chip_overlap_peaks.py --in1=$rep_tagAlign --peak1=$peak1 --peak2=$peak2 --peakpool=$peakpool \
    --xcor1=$xcor1 --pair=False --chromsizes=$chromsizes --peaktype="gappedPeak" --peakas=$gappedpeakAs \
    --signal1="peaks_macs/r1.fc_signal.bw" --prefix="${peak_dirname}/final" 1>${exp_basename}.final.gappedpeak.log1 2>${exp_basename}.final.gappedpeak.log2 
unset peak1 peak2 peakpool signal1 rep_tagAlign xcor1
unset reference chromsizes narrowpeakAs gappedpeakAs broadpeakAs exp_fq input_fq dirname exp_filename input_filename exp_basename input_basename peak_dirname

