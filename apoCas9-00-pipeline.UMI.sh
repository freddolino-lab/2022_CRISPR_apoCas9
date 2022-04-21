#!/usr/bin/env bash
set -e

## Defining the path to scripts for this data analysis project
dataBatchName='202111'
dataTypeName='adaptationUMI'

projectDir='/home/diaorch/projects/crisprAdaptation/crisprAdaptation-'"$dataBatchName"'/'
scriptCodeName='crispr'"$dataBatchName"'_'"$dataTypeName"'-'
source "$projectDir""$scriptCodeName"'00-envSetUp.sh'

# Setting up / activating conda environment 
conda activate bio_py-3.7

# Step 00 to Step 01 are the same as non-UMI sample analysis

# Step 02:
## trimming for spacer
echo '=== Step 02: trimming ==='
LEADS_R26='cttatgaaataaggatttcccgtcgaagtattgtagcactgcgaaatgagaaagggagctacaac' # LeadS + 5 b gap + Repeat 26
R26='attgtagcactgcgaaatgagaaagggagctacaac' # Repeat 26
pSp25='gtaaaggtaatgcgccgcgcatg' # partial Spacer 25
P7_PREFIX='actacgcacgcgacga'
UMI_LENGTH=8

mkdir "$newSpacerSubDir"
mkdir "$umiSubDir"
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
    echo '=== cutting partial Spacer 25 at 3 prime ===' # also filter out empty reads
    "$CUTADAPT_PATH" -a "$pSp25" --minimum-length 1 -o "$newSpacerSubDir""$dataPrefix""$i""$cut2SeqFqSuffix" "$newSpacerSubDir""$dataPrefix""$i""$cut1SeqFqSuffix" > "$newSpacerSubDir""$dataPrefix""$i"'.cut2.log' 2> "$newSpacerSubDir""$dataPrefix""$i"'.cut2.err'
    echo '===     getting new Spacer sequences     ===' # by removing R26 at 3 prime for adapted reads
    "$CUTADAPT_PATH" -a "$R26" --discard-untrimmed -o "$newSpacerSubDir""$dataPrefix""$i""$removeSeqFqSuffix" "$newSpacerSubDir""$dataPrefix""$i""$cut2SeqFqSuffix" > "$newSpacerSubDir""$dataPrefix""$i"'.remove.log' 2> "$newSpacerSubDir""$dataPrefix""$i"'.remove.err'
    echo '=====              UMI                 ====='
    r2Filename="${filename/R1/R2}"
    echo '===        cutting p7 at 5prime          ==='
    "$CUTADAPT_PATH" -g ^"$P7_PREFIX" --discard-untrimmed -o "$umiSubDir""$dataPrefix""$i""$cutUmiSeqFqSuffix" "$r2Filename" > "$umiSubDir""$dataPrefix""$i"'.cutUmi.log' 2> "$umiSubDir""$dataPrefix""$i"'.cutUmi.err'
    echo '===        getting UMI sequences         ==='
    "$CUTADAPT_PATH" --length  "$UMI_LENGTH" --minimum-length "$UMI_LENGTH" -o "$umiSubDir""$dataPrefix""$i""$umiSeqFqSuffix" "$umiSubDir""$dataPrefix""$i""$cutUmiSeqFqSuffix" > "$umiSubDir""$dataPrefix""$i"'.umi.log' 2> "$umiSubDir""$dataPrefix""$i"'.umi.err'
done

# Step 03:
# deduplicating spacer-umi combination
echo '=== Step 03: deduplication ==='
mkdir "$dedupSubDir"
NTHREAD=4
NZIP=9
for filename in "$rawSubDir""$dataPrefix"*"$r1Suffix""$rawSuffix";
do
    echo $filename
    filenameBase=${filename#$rawSubDir$dataPrefix}
    i=${filenameBase%$r1Suffix$rawSuffix}
    echo "$i"
    echo 'Sample '"$i"
    spacerSeqFilename="$newSpacerSubDir""$dataPrefix""$i""$removeSeqFqSuffix"
    echo "$spacerSeqFilename"
    umiSeqFilename="$umiSubDir""$dataPrefix""$i""$umiSeqFqSuffix"
    echo "$umiSeqFilename"
    "$PARDRE_PATH" -i "$spacerSeqFilename" -p "$umiSeqFilename" -o "$dedupSubDir""$dataPrefix""$i""$spacerDedupSeqFqSuffix" -r "$dedupSubDir""$dataPrefix""$i""$umiDedupSeqFqSuffix" -t "$NTHREAD" -z "$NZIP" > "$dedupSubDir""$dataPrefix""$i"'.dedup.log' 2> "$dedupSubDir""$dataPrefix""$i"'.dedup.err'
done

# Step 04 to Step 06 are the same as common non-UMI sample analysis

# Final:
# print success message
echo 'Script run was successful.'
