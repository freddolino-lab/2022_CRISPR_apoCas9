#!/usr/bin/env bash

## ===== data organization =====
echo "$HOME"
##### CUSTOMIZE HERE (2 lines) #####
batchAccn='202203'
dataType='adaptationNonUMI' 
projectDir="$HOME"'/projects/crisprAdaptation/crisprAdaptation-'"$batchAccn"'/'
## ===== listing of sample accessions of different types =====
##### CUSTOMIZE HERE (2+ lines: sampleNumberSet) #####
sampleNumberSet=$(seq -f "%02g" 1 1 6)
batchDataDir="$HOME"'/data/crisprAdaptation-'"$batchAccn"
analysisDataDir="$batchDataDir"'/'"$dataType"'/'
scriptCodeName='crispr'"$batchAccn"'_'"$dataType"'-'
dataPrefix='crispr_'"$batchAccn"'-'


## ===== program paths =====
CONDA_PATH="$HOME"'/local/miniconda3/envs/apoCas9-py-3.7/bin/'
CUTADAPT_PATH="$CONDA_PATH"'cutadapt'
PARDRE_PATH="$HOME"'/src/ParDRe-rel2.2.5/ParDRe'
BOWTIE2_PATH='/usr/bin/bowtie2'
BOWTIE2BUILD_PATH='/usr/bin/bowtie2-build'
SAMTOOLS_PATH='/usr/bin/samtools'
BEDTOOLS_PATH="$HOME"'/src/bedtools2/bin/bedtools'
BEDTOOLS_COVERAGEBED="$HOME"'/src/bedtools2/bin/genomeCoverageBed'

GENOMEFASTA="$HOME"'/genomes/N_meningitidis/8013/bowtie2idx/GCF_000026965.1_ASM2696v1_genomic.fna'
BOWTIE2INDEX="$HOME"'/genomes/N_meningitidis/8013/bowtie2idx/GCF_000026965.1_ASM2696v1_genomic'

## ===== homemade script paths =====
padBedFilesByLength="$HOME"'/projects/2022_CRISPR_apoCas9/apoCas9-05-getPaddedBed.py'
getFastaFromBed="$HOME"'/projects/2022_CRISPR_apoCas9/apoCas9-05-getFastaFromBed.py'
getReverseComplement="$HOME"'/projects/2022_CRISPR_apoCas9/apoCas9-06-revCompFasta.py'
calcPlotPwm="$HOME"'/projects/2022_CRISPR_apoCas9/apoCas9-06-calcPwm.py'

##  ==== references          ====
refSubDir="$analysisDataDir"'00-ref/'
##  ==== raw data processing ====
rawSubDir="$analysisDataDir"'00-seq/'
filenameTablePath="$projectDir$batchAccn"'-experiments-filenameTable.csv'
rawSuffix='.fastq'

r1Suffix='_R1'
r2Suffix='_R2'

## ==== quality control ====
qcSubDir="$analysisDataDir"'01-qc/'

## ==== spacer and UMI extraction ====
## spacer
newSpacerSubDir="$analysisDataDir"'02-spacer/'
cut1SeqFqSuffix='.cut1.seq.fastq.gz'
cut2SeqFqSuffix='.cut2.seq.fastq.gz'
removeSeqFqSuffix='.remove.seq.fastq.gz'
# ## UMI
umiSubDir="$analysisDataDir"'02-umi/'
cutUmiSeqFqSuffix='.cutUmi.seq.fastq.gz'
umiSeqFqSuffix='.umi.seq.fastq.gz'

## ==== deduplication ====
dedupSubDir="$analysisDataDir"'03-dedup/'
spacerDedupSeqFqSuffix='.spacer.dedupped.seq.fastq.gz'
umiDedupSeqFqSuffix='.umi.dedupped.seq.fastq.gz'

## ==== alignment ====
alignSubDir="$analysisDataDir"'04-align/'
alignMidfix='.newSpacer.align'
bt2BuildLog="$analysisDataDir"'04.bowtie2-build.log'
alignLog="$analysisDataDir"'04.align.all.log'

## ==== padded genome positions ====
padSubDir="$analysisDataDir"'05-pad/'
paddedFaSuffix='.padded.fasta'

## ==== reverse complement of extracte sequence, PWM ====
rcSubDir="$analysisDataDir"'061-rc/'
pwmSubDir="$analysisDataDir"'062-pwm/'
rcPaddedFaSuffix='.padded.rc.fasta'

