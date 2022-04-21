#!/usr/bin/env bash
set -e

## Defining the path to scripts for this data analysis project
dataBatchName='202203'
dataTypeName='adaptationNonUMI'

projectDir='/home/diaorch/projects/crisprAdaptation/crisprAdaptation-'"$dataBatchName"'/'
scriptCodeName='crispr'"$dataBatchName"'_'"$dataTypeName"'-'
source "$projectDir""$scriptCodeName"'00-envSetUp.sh'

# Setting up / activating conda environment 
conda activate bio_py-3.7

# Step 00 to Step 01 are the same as common non-UMI sample analysis

# Step 02:
## trimming for spacer
echo '=== Step 02: trimming ==='
# CHANGES IN SEQUENCING READ STRUCTURE DUE TO CHANGES IN SEQ PLATFORM AND PROVIDER
# sequencing batch 202203
# Fig 5D Cje sample: only part of Repeat 26 falls in the sequencing read
R26='attgtagcactgcgaaatgagaa' 

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
    echo '===     getting new Spacer sequences     ==='
    "$CUTADAPT_PATH" -a "$R26" --discard-untrimmed --minimum-length 28 --maximum-length 32 -o "$newSpacerSubDir""$dataPrefix""$i""$removeSeqFqSuffix" "$newSpacerSubDir""$dataPrefix""$i""$cut1SeqFqSuffix" > "$newSpacerSubDir""$dataPrefix""$i"'.remove.log' 2> "$newSpacerSubDir""$dataPrefix""$i"'.remove.err'
done

# Step 04 to Step 06 are the same as common non-UMI sample analysis

# Final:
# print success message
echo 'Script run was successful.'
