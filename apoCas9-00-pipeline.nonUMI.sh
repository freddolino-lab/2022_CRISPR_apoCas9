#!/usr/bin/env bash
set -e

## Defining the path to scripts for this data analysis project
dataBatchName='202203'
dataTypeName='adaptationNonUMI'
### CUSTOMIZE ENDS ###

projectDir='/home/diaorch/projects/crisprAdaptation/crisprAdaptation-'"$dataBatchName"'/'
scriptCodeName='crispr'"$dataBatchName"'_'"$dataTypeName"'-'
source "$projectDir""$scriptCodeName"'00-envSetUp.sh'

# Setting up / activating conda environment 
conda activate bio_py-3.7

# Step 00:
# Creating data directory
# manually create analysis data directory so that the pipeline.log and .err can be stored for running this script
# mkdir "$analysisDataDir"

# Step 00:   
# Linking and organizing raw sequencing data files
# create symbolic links for raw data files (.fastq.gz) in working directory (dataDir/00-seq)
# mkdir "$rawSubDir"
# ln -s "$filename" "$rawSubDir""$renaming""$fileSuffix"

# Step 01:  
# Performing quality control for raw data
echo '=== Step 01: quality control ==='
# echo 'Step 01: Performing quality control on raw reads'
mkdir "$qcSubDir"
## run FastQC
mkdir "$qcSubDir"'fastqc/'
fastqc "$rawSubDir""$dataPrefix"*"$rawSuffix" -o "$qcSubDir"'fastqc/' --noextract
## run MultiQC
mkdir "$qcSubDir"'multiqc/'
multiqc "$qcSubDir"'fastqc/'* -o "$qcSubDir"'multiqc/'

# Step 02:
## trimming for spacer
echo '=== Step 02: trimming ==='
LEADS_R26='cttatgaaataaggatttcccgtcgaagtattgtagcactgcgaaatgagaaagggagctacaac' # LeadS + 5 b gap + Repeat 26
R26='attgtagcactgcgaaatgagaaagggagctacaac' # Repeat 26
PAIR1='acacagtgcaattatggcaagaaaggcgacagatattgtgtc' # "Pair 1"

mkdir "$newSpacerSubDir"
# mkdir "$umiSubDir"
for filename in "$rawSubDir""$dataPrefix"*"$r1Suffix""$rawSuffix";
do
    echo $filename
    filenameBase=${filename#$rawSubDir$dataPrefix}
    i=${filenameBase%$r1Suffix$rawSuffix}
    echo "$i"
    echo 'Sample '"$i"
    echo '=====             Spacer               ====='
    echo '=== cutting LeadS + Repeat 26 at 5 prime ==='
    "$CUTADAPT_PATH" -g "$LEADS_R26" --discard-untrimmed -o "$newSpacerSubDir""$dataPrefix""$i""$cut1SeqFqSuffix" "$filename" > "$newSpacerSubDir""$dataPrefix""$i"'.cut1.log' 2> "$newSpacerSubDir""$dataPrefix""$i"'.cut1.err'
    echo '=== cutting partial Spacer 25 at 3 prime ==='
    "$CUTADAPT_PATH" -a "$PAIR1" --minimum-length 1 -o "$newSpacerSubDir""$dataPrefix""$i""$cut2SeqFqSuffix" "$newSpacerSubDir""$dataPrefix""$i""$cut1SeqFqSuffix" > "$newSpacerSubDir""$dataPrefix""$i"'.cut2.log' 2> "$newSpacerSubDir""$dataPrefix""$i"'.cut2.err'
    echo '===     getting new Spacer sequences     ==='
    "$CUTADAPT_PATH" -a "$R26" --discard-untrimmed -o "$newSpacerSubDir""$dataPrefix""$i""$removeSeqFqSuffix" "$newSpacerSubDir""$dataPrefix""$i""$cut2SeqFqSuffix" > "$newSpacerSubDir""$dataPrefix""$i"'.remove.log' 2> "$newSpacerSubDir""$dataPrefix""$i"'.remove.err'
    "$CUTADAPT_PATH" -a "$R26" --discard-untrimmed --minimum-length 28 --maximum-length 32 -o "$newSpacerSubDir""$dataPrefix""$i""$removeSeqFqSuffix" "$newSpacerSubDir""$dataPrefix""$i""$cut1SeqFqSuffix" > "$newSpacerSubDir""$dataPrefix""$i"'.remove.log' 2> "$newSpacerSubDir""$dataPrefix""$i"'.remove.err'
done

# Step 04:
# aligning to N. meningitidis genome
echo '=== Step 04: alignment ==='
# echo '==== Building Bowtie2 index ===='
# "$BOWTIE2BUILD_PATH" "$GENOMEFASTA" "$BOWTIE2INDEX" > "$bt2BuildLog"
echo '==== Aligning ===='
readPath="$newSpacerSubDir"
readSuffix="$removeSeqFqSuffix"
logFullFilename="$alignLog"

mkdir "$alignSubDir"
printf '' > "$logFullFilename"

for filename in "$readPath"*"$readSuffix"
do
    tmp=${filename%"$readSuffix"}
    filenameBase=${tmp#"$readPath"}
    echo 'Sample name (base name): '"$filenameBase"
    echo 'Sample name (base name): '"$filenameBase" >> "$logFullFilename"
    echo '=== Aligning to custom genome ==='
    "$BOWTIE2_PATH" --very-sensitive -x "$BOWTIE2INDEX" -U "$filename" -S "$alignSubDir""$filenameBase""$alignMidfix".sam 2>> "$logFullFilename"
    echo "=== Converting SAM to sorted BAM ==="
    "$SAMTOOLS_PATH" view -Sbh "$alignSubDir""$filenameBase""$alignMidfix"'.sam' -o "$alignSubDir""$filenameBase""$alignMidfix"'.bam'
    "$SAMTOOLS_PATH" sort -o "$alignSubDir""$filenameBase""$alignMidfix"'.sorted.bam' "$alignSubDir""$filenameBase""$alignMidfix"'.bam' 
    rm "$alignSubDir""$filenameBase""$alignMidfix"'.sam'
    rm "$alignSubDir""$filenameBase""$alignMidfix"'.bam'
    "$SAMTOOLS_PATH" index "$alignSubDir""$filenameBase""$alignMidfix"'.sorted.bam'
    echo "=== Filter sorted bam file for uniquely mapped reads ==="
    echo "===== forward ====="
    "$SAMTOOLS_PATH" view -h -F 20 "$alignSubDir""$filenameBase""$alignMidfix".'sorted.bam' | grep -v 'XS:i' > "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.forward.sam'
    if grep -v -q '^@' "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.forward.sam';
    then
        "$SAMTOOLS_PATH" view -Sbh "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.forward.sam' -o "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.forward.bam'
    else
        echo 'ERROR: no uniquely aligned forward reads for sample '"$filenameBase"
        exit 1
    fi
    rm "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.forward.sam'
    "$SAMTOOLS_PATH" index "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.forward.bam'
    echo "===== reverse ====="
    "$SAMTOOLS_PATH" view -h -f 16 -F 4 "$alignSubDir""$filenameBase""$alignMidfix"'.sorted.bam' | grep -v 'XS:i' > "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.reverse.sam'
    if grep -v -q '^@' "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.reverse.sam';
    then
        "$SAMTOOLS_PATH" view -Sbh "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.reverse.sam' -o "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.reverse.bam'
    else
        echo 'ERROR: no uniquely aligned reverse reads for sample '"$filenameBase"
        exit 1
    fi
    rm "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.reverse.sam'
    "$SAMTOOLS_PATH" index "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.reverse.bam'
done

# Step 05:
# getting padded sequences
mkdir "$padSubDir"
echo '=== Step 05: padding ==='

for filename in "$alignSubDir"*"$alignMidfix"'.uniq.forward.bam'
do
    echo "$filename"
    tmp=${filename%%"$alignMidfix"*}
    filenameBase=${tmp##"$alignSubDir"}
    echo 'Sample name (base name): '"$filenameBase"
    echo '=== Padding forward alignments ==='
    "$BEDTOOLS_COVERAGEBED" -bg -5 -ibam "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.forward.bam' >"$padSubDir""$filenameBase"'.forward.bedgraph'
    python2 "$padBedFilesByLength" "$padSubDir""$filenameBase"'.forward.bedgraph' 25 50 forward >"$padSubDir""$filenameBase"'.forward.padded.bed'
    # reverse
    echo "=== Padding reverse alignments ==="
    "$BEDTOOLS_COVERAGEBED" -bg -5 -ibam "$alignSubDir""$filenameBase""$alignMidfix"'.uniq.reverse.bam' >"$padSubDir""$filenameBase"'.reverse.bedgraph'
    python2 "$padBedFilesByLength" "$padSubDir""$filenameBase"'.reverse.bedgraph' 25 50 reverse >"$padSubDir""$filenameBase"'.reverse.padded.bed'
    # merge
    echo "=== Merging BED files of padded locations ==="
    cat "$padSubDir""$filenameBase"'.forward.padded.bed' "$padSubDir""$filenameBase"'.reverse.padded.bed' >"$padSubDir""$filenameBase"'.padded.bed'
    echo "=== Retrieving genomic sequences by BED into FASTA ==="
    python "$getFastaFromBed" "$padSubDir""$filenameBase"'.padded.bed' "$GENOMEFASTA" >"$padSubDir""$filenameBase""$paddedFaSuffix"
done

# Step 06
# visualizaing PAMs
# Part 1 and 2
# reverse complementing extracted sequences, calculate PWM
mkdir "$rcSubDir"
mkdir "$pwmSubDir"
echo '=== Step 06: visualization ==='
for filename in "$padSubDir"*"$paddedFaSuffix"
do
    echo "$filename"
    tmp=${filename%%"$paddedFaSuffix"}
    filenameBase=${tmp#"$padSubDir"}
    echo 'Sample name (base name): '"$filenameBase"
    # Part 1: reverse complement 
    python3 "$getReverseComplement" "$filename" "$rcSubDir""$filenameBase""$rcPaddedFaSuffix"
    # Part 2: PWM calculation and plotting
    python3 "$calcPlotPwm" "$rcSubDir""$filenameBase""$rcPaddedFaSuffix" "$pwmSubDir"
done

# Final:
# print success message
echo 'Script run was successful.'
