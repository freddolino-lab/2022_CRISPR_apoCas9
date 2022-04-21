
library('reshape2')
library('ggplot2')
library('RColorBrewer')
library('gridExtra')
library('grid')
# rm(list = ls())
`%strcat%` <- function(a, b) paste(a, b, sep = '')

## modification block starts
batchAccn <- '202106'
## modification block ends

dataPrefix <- 'crispr_' %strcat% batchAccn %strcat% '-'
homePath <- Sys.getenv("HOME") %strcat% '/'
dataPath <- homePath %strcat% 'data/crisprAdaptation-' %strcat% batchAccn %strcat% '/'
filepath <- dataPath %strcat% c('adaptationUMI/', 'adaptationNonUMI/') %strcat% '062-pwm/'
filenameList <- dir(path = filepath, pattern = '.csv', full.names = TRUE)
lenSeq <- 75
# print(filenameList)
outputPath <- dataPath %strcat% 'analysis/'

collapseDominatingBase <- function(pwmFormatted){
  domBases <- vector(mode = 'character', length = ncol(pwmFormatted))
  domPercentages <- vector(mode = 'numeric', length = ncol(pwmFormatted))
  alphabet <- rownames(pwmFormatted)
  for (i in 1:ncol(pwmFormatted)){
    domRow <- which.max(pwmFormatted[, i])
    domBases[i] <- alphabet[domRow]
    domPercentages[i] <- max(pwmFormatted[, i])
  }
  return(list(domBases, domPercentages))
}

getSampleName <- function(filename){
  sn <- strsplit(x = basename(filename), split = '.', fixed = TRUE)[[1]][1]
  return(sn)
}

# get dominating base information 
dataDf <- data.frame(matrix(data = NA, nrow = 0, ncol = 4))
colnames(dataDf) <- c('sample', 'position', 'predominantBase', 'predominantPercentage')
for (filename in filenameList){
  print(filename)
  pwm <- read.csv(file = filename, sep = ',', header = FALSE, stringsAsFactors = FALSE)
  rownames(pwm) <- c('A', 'C', 'G', 'T')
  predominantBase <- collapseDominatingBase(pwm)
  sampleName <- getSampleName(filename)
  singleSampleDf <- cbind.data.frame(c(sampleName), 
                                     c(1:ncol(pwm)), 
                                     predominantBase[[1]], 
                                     predominantBase[[2]], 
                                     stringsAsFactors = FALSE)
  colnames(singleSampleDf) <- c('sample', 'position', 'predominantBase', 'predominantPercentage')
  dataDf <- rbind(dataDf, singleSampleDf)
}
# print(head(dataDf))

dataDf <- dataDf[order(dataDf[, 'sample'], dataDf[, 'position']), ]
dataDf$sample <- gsub(x = dataDf$sample, 
                      pattern = dataPrefix, replacement = '', fixed = TRUE)

# calculating (as a function) bit entropy for one position
b <- 2 # base for log in entropy, base 2 corresponding to the unit "bits"
calculateWeightVectorEntropy <- function(w, # weight vector
                                         b = 2 # log base in entropy
                                         ){
  s <- -1 * sum(w * log(x = w, base = b), na.rm = TRUE)
  return(s)
}

# calculating bit entropy for PAM positions
pamPositionLeft <- 55
pamPositionRight <- 58
e <- c()
eNames <- c()
for (filename in filenameList){
  print(filename)
  pwm <- read.csv(file = filename, sep = ',', header = FALSE, stringsAsFactors = FALSE)
  rownames(pwm) <- c('A', 'C', 'G', 'T')
  sampleName <- getSampleName(filename)
  eNames <- c(eNames, sampleName)
  pamPositionMtx <- pwm[, pamPositionLeft:pamPositionRight]
  pamPositionEntropySet <- apply(X = pamPositionMtx, MARGIN = 2, 
                                 FUN = calculateWeightVectorEntropy)
  e <- c(e, sum(pamPositionEntropySet))
}
pamEntropy <- data.frame(sampleName = eNames, entropy = e, stringsAsFactors = FALSE)
pamEntropy <- pamEntropy[order(pamEntropy[, 'sampleName']), ]
print(pamEntropy)
write.table(x = pamEntropy, 
            file = 'analysis/pamEntropy.raw.csv', 
            col.names = TRUE, row.names = FALSE, quote = TRUE, sep = ',')
