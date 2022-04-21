batchAccn <- '202106'
dataPath <- paste(Sys.getenv('HOME'), '/data/crisprAdaptation-', batchAccn, '/', sep = '')

pamEntropyTbl <- read.table(file = paste(dataPath, 'analysis/', batchAccn, '-pamEntropy.weighted.csv', sep = ''), head = TRUE, sep = ',', colClasses = c('character', 'numeric'))
pamEntropyTbl[, 'batch'] <- sapply(X = pamEntropyTbl[, 'sampleName'], FUN = function(x) gsub(x = strsplit(x = x, split = '-', fixed = TRUE)[[1]][1], pattern = 'crispr_', replacement = '', fixed = TRUE))
pamEntropyTbl[, 'sample'] <- sapply(X = pamEntropyTbl[, 'sampleName'], FUN = function(x) strsplit(x = x, split = '-', fixed = TRUE)[[1]][2])

pamEntropyBoundsTbl <- read.table(file = paste(dataPath, 'analysis/', batchAccn, '-pamEntropy.standards.csv', sep = ''), head = TRUE, sep = ',', colClasses = c('character', rep('numeric', 2)))
print(pamEntropyBoundsTbl)

b <- batchAccn
l <- pamEntropyBoundsTbl[pamEntropyBoundsTbl[, 'batch'] == b, 'low']
u <- pamEntropyBoundsTbl[pamEntropyBoundsTbl[, 'batch'] == b, 'high']
pamScaled <- (pamEntropyTbl[, 'entropy'] - l) / (u - l) 
print(summary(pamScaled))
pamEntropyTbl[, 'normed_entropy'] <- pamScaled

resTbl <- pamEntropyTbl[, c('batch', 'sample', 'entropy', 'normed_entropy')]
colnames(resTbl)[3] <- 'pam_entropy'
print(resTbl)
write.table(resTbl, 
            file = 'analysis/pamEntropy.norm.csv', 
            sep = ',', col.names = TRUE, row.names = FALSE, quote = TRUE)

