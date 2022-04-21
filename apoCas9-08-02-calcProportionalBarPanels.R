rm(list = ls())
`%strcat%` <- function(x, y) paste(x, y, sep = '')
concat <- function(...){
  paste(list(...), collapse = '')
}

# library('ggplot2')
library('Biostrings', verbose = F)
source('apoCas9-08-00-funcLib.R', sep = '')
dataPath <- '~/data/apoCas9/plots/'

#### input processing modified from heat map plotting ####
# # define plots of interest and position ranges for PAM
basePamRange <- c(55:58)
figSampleTbl <- read.table(file = dataPath %strcat% 'ref/figure_sample_list.csv', 
                           sep = ',', colClasses = 'character', header = T)
# print(unique(figSampleTbl[, 'figure']))
proportionFigSet <- c('Fig.S10A', 'Fig.S10B')
# print(sum(figSampleTbl[, 'figure'] == proportionFigSet))

GENOMESIZE <- c('NC_017501.1' = 2277550)
kmerRefFilename <- paste(dataPath, '00-ref/Nm_genomic_4mers.tsv', sep = '')
# print(kmerRefFilename)

baseAlphabet <- c('A', 'C', 'G', 'T')
allPamPrm <- apply(X = gtools::permutations(n = 4, r = 4, v = baseAlphabet, 
                                            set = TRUE, repeats.allowed = TRUE), 
                   MARGIN = 1, FUN = function(x)paste0(x, collapse=''))
# print(allPamPrm)
# print(length(allPamPrm))

countPamFromFasta <- function(fn, chrom, 
                              pamStart = min(basePamRange), 
                              pamLen = length(basePamRange )){
  faSet <- readDNAStringSet(fn)
  faChromSet <- faSet[grepl(pattern = chrom, x = names(faSet))]
  # print(faChromSet)
  # # get PAM
  pam <- as.character(subseq(x = faChromSet, start = pamStart, width = pamLen))
  # # get insertion counts for the position
  w <- as.numeric(gsub(pattern = '.*\\[(.*)\\].*', replacement = '\\1', x = names(x = faChromSet), perl = TRUE))
  # get PAM sum counts
  w_df <- aggregate(x = w, by = list(pam), FUN = sum) # weighted; insertion events
  p_df <- as.data.frame(table(pam)) # counts not weighted; insertion positions
  df <- merge(x = w_df, y = p_df, by.x = 'Group.1', by.y = 'pam', all = TRUE, sort = TRUE)
  colnames(df) <- c('pam', 'event', 'site')
  dfOrdered <- df[match(x = allPamPrm, table = df), ]
  return(df)
}

for (figName in proportionFigSet){
  print(figName)
  subFigSampleTbl <- figSampleTbl[figSampleTbl[, 'figure'] == figName, ]
  print(subFigSampleTbl)
  
  sampleFullnameSet <- vector('character', nrow(subFigSampleTbl))
  for (i in 1:nrow(subFigSampleTbl)){
    sampleFullname <- extractSampleFullName(
      batchPathSet = formatPwmDataPath(batch = subFigSampleTbl[i, 'computational']), 
      sampleNumbering = subFigSampleTbl[i, 'sample_numbering'], 
      inputMode = 'rc')
    sampleFullnameSet[i] <- sampleFullname
  }
  print(sampleFullnameSet)
  
  chrom <- 'NC_017501.1'
  for (s in sampleFullnameSet){
    print(s)
    
    pamCount <- countPamFromFasta(fn = s, chrom = chrom)
    missingPam <- allPamPrm[! allPamPrm %in% pamCount[, 'pam']]
    if (length(missingPam) > 0){
      pamCountNorm <- rbind(pamCount, 
                            data.frame('pam' = missingPam, 
                                       'event' = 0, 'site' = 0))
    } else {
      pamCountNorm <- pamCount
    }
    pamCountNorm[, 'eventProportion'] <- pamCountNorm[, 'event'] / sum(pamCountNorm[, 'event'])
    pamCountNorm[, 'siteProportion'] <- pamCountNorm[, 'site'] / sum(pamCountNorm[, 'site'])
  
    pamCountRes <- data.frame(
      sample = getSampleId(s), 
      chrom = chrom, 
      pamCountNorm
    )
    pamCountResFormat <- pamCountRes
    pamCountResFormat[, c('eventProportion', 'siteProportion')] <- 
      round(x = pamCountResFormat[, c('eventProportion', 'siteProportion')], 
            digits = 4)
    
    pamCountResFormat <- pamCountResFormat[order(pamCountResFormat[, 'eventProportion'], decreasing = TRUE), ]
    print(head(pamCountResFormat, 10))
    colnames(pamCountResFormat) <- c('Sample', 'Uptake_Source', 'PAM', 
                               'Event_Count', 'Position_Count', 
                               'Event_Proportion', 'Postition_Proportion'
    )
    print(head(pamCountResFormat))
    print(tail(pamCountResFormat))
    
    write.csv(x = pamCountResFormat, 
              file = concat('proportion/pam_proportion_tbl/', getSampleId(s),'.pam_proportion.csv'),  
              row.names = FALSE
    )
  }
}
