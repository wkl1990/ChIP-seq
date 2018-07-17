# ChIP-seq analysis
## Description
The codes contain the basic processes for ChIP-seq analysis, which contain QC, mapping, and peak calling. Narrow peak is usually for sequence-specific transcription factors, their cofactors, and, with some caveats, transcription start site or enhancer-associated histone marks. Broad peak contain certain chromatin marks (H3K9me3, H3K36me3, etc.) and chromatin proteins associated with transcriptional elongation or repression.
## Softwares and packages
* Quality control: [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* Mapping: [BWA](http://bio-bwa.sourceforge.net/)
* Peak calling: [MACS2](https://github.com/taoliu/MACS/), [BCP](https://cb.utdallas.edu/BCP/), [MUSIC](https://github.com/gersteinlab/MUSIC)
* Annotation and visualization: [ChIPseeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html).
## Pipeline
ENCODE3-pipe is a revision of [ENCODE3 pipeline](https://github.com/ENCODE-DCC/chip-seq-pipeline).
![ChIP-seq pipeline](none)
