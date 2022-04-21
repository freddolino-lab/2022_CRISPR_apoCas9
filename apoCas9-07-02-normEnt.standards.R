# rm(list = ls())
batchAccn <- '202106'
dataPath <- paste(Sys.getenv('HOME'), '/data/crisprAdaptation-', batchAccn, '/', sep = '')

# PAM entropy values for one batch of sequencing data
pamEntropyTbl <- read.table(file = paste(dataPath, 'analysis/20211006-crispr_', batchAccn, '-pamEntropy.weighted.csv', sep = ''), head = TRUE, sep = ',', colClasses = c('character', 'numeric'))
# print(head(pamEntropyTbl202106))
sampleInfo <- gsub(x = pamEntropyTbl202106[, 'sampleName'], pattern = 'crispr_202106-', replacement = '')
sampleAccn <- sapply(X = sampleInfo, FUN = function(x) strsplit(x = x, split = '-', fixed = TRUE)[[1]][1])
names(sampleAccn) <- NULL
# print(sampleAccn)

# # PAM entropy values, combined all batches
# pamEntropyTbl <- pamEntropyTbl202103
pamEntropyVec <- pamEntropyTbl202106[, 'entropy']
names(pamEntropyVec) <- sampleAccn
print(pamEntropyVec)

# read in sample annotations
annotTbl <- read.table(file = paste(dataPath, '00-ref/202106-experiments-sampleLabel.csv', sep = ''), 
                       sep = ',', header = TRUE, colClasses = 'character', stringsAsFactors = FALSE)
# print(head(annotTbl))

# calculate the standards ("bounds")
strainStandard <- c('L' = '299', 'H' = '308') # low/high entropy strains based on genotype
# print(annotTbl[annotTbl[, 'strain'] %in% strainStandard, ])
batchSampleStdH <- annotTbl[annotTbl[, 'strain'] == strainStandard['H'], 'sample']
batchSampleStdL <- annotTbl[annotTbl[, 'strain'] == strainStandard['L'] & annotTbl[, 'IPTG'] == '1' & annotTbl[, 'aTc'] == '2', 'sample']
print(list(batchSampleStdH, batchSampleStdL))
# print(pamEntropyVec)
# high entropy
if (length(batchSampleStdH) > 0){
  batchEntStdH <- mean(pamEntropyVec[batchSampleStdH])
}
print(batchEntStdH)
# low entropy
if (length(batchSampleStdL) > 0){
  batchEntStdL <- mean(pamEntropyVec[batchSampleStdL])
}
print(batchEntStdL)

boundsDf <- data.frame('batch' = batchAccn,
                       'high' = batchEntStdH,
                       'low' = batchEntStdL, 
stringsAsFactors = FALSE)
print(str(boundsDf))
write.table(boundsDf, 
            file = 'analysis/pamEntropy.standards.csv', 
            sep = ',', col.names = TRUE, row.names = FALSE)
