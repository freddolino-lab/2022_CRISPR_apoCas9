# Project/code: CRISPR adaptation
# Start date: Sun Dec  1 19:48:16 2019
# Objective: visualize the utilization of different 4mer PAMs
# --------------
# Author: diaorch
# Modification date:  Sun Dec  1 19:48:16 2019
# Modification (if modified from another procjet/script):
#   0. original (not modified from other sources)
#   1. 20200706: copied and mortified from ~/projects/crisprAdaptation/crisprAdaptation-201902-tracrCas9_mismatch0/crispr201902_mism0-06-pamProportion.viz.R for plotting PAM proportion bar charts for manuscript
#   2. 20200712:  copied and mortified from ~/projects/crisprAdaptation/crispr-manuscript/202007/crisprWrite202007-01-plotProportionalBarPanels.3G.R to clean up and adapt to another set of samples
#   3. 20211003:  copied and mortified from ~/projects/crisprAdaptation/crispr-manuscript/202007/crisprWrite202007-01-plotProportionalBarPanels.R to replot and write summarized version of proportion tables
# --------------

rm(list = ls())

library('gtools')
library('Biostrings')
library('ggplot2')
library('egg')
library('grid')

#### input processing modified from heat map plotting ####
# # define plots of interest and position ranges for PAM
basePamRange <- c(55:58)
figSampleTbl <- read.table(file = dataPath %strcat% 'ref/figure_sample_list.csv', 
                           sep = ',', colClasses = 'character', header = T)
# print(unique(figSampleTbl[, 'figure']))
proportionFigSet <- c('Fig.S10A', 'Fig.S10B')
# print(sum(figSampleTbl[, 'figure'] == proportionFigSet))

inputPath <- dataPath %strcat% 'pam_proportion_tbl'
proportionTblFilenameSet <- list.files(path = inputPath, pattern = '.pam_proportion.csv', full.name = T)
getProportionSampleId <- function(filename){
  sn <- strsplit(x = basename(filename), split = '.', fixed = TRUE)[[1]][1] # sample name
  batchInfo <- strsplit(x = sn, split = '-', fixed = TRUE)[[1]][3]
  sampleNumbering <- strsplit(x = strsplit(x = sn, split = '-', fixed = TRUE)[[1]][4], split = '_', fixed = TRUE)[[1]][1]
  return(paste(batchInfo, sampleNumbering, sep = '-'))
}
names(proportionTblFilenameSet) <- sapply(
  X = proportionTblFilenameSet, FUN = getProportionSampleId)
print(proportionTblFilenameSet)

baseAlphabet <- c('A', 'C', 'G', 'T')
allPamPrm <- apply(X = gtools::permutations(n = 4, r = 4, v = baseAlphabet, 
                                            set = TRUE, repeats.allowed = TRUE), 
                   MARGIN = 1, FUN = function(x)paste0(x, collapse=''))
# print(allPamPrm)
# print(length(allPamPrm))

exactGATT <- 'GATT'
nearGATT <- c('GACT', 'GCTT', 'GTTT', 'GATA')  
pamGroup <- rep('non-optimal', length(allPamPrm))
names(pamGroup) <- allPamPrm
pamGroup[exactGATT] <- 'optimal'
pamGroup[nearGATT] <- 'sub-optimal'
# print(pamGroup)
# print(table(pamGroup))

GENOMESIZE <- c('NC_017501.1' = 2277550)
genomicBackground <- read.table(
  file = 'ref/genomic_kmer_content-NC_017501.1.csv'
  header = TRUE, sep = ',', quote = '"')
colnames(genomicBackground) <- c('kmer', 'kmer_count', 
                                 'rc_kmer', 'rc_kmer_count', 
                                 'total_count')
# print(head(genomicBackground))
kmerOrder <- genomicBackground$kmer 
# print(length(kmerOrder))

rulers <- data.frame(
  'y' = c(
          sum(genomicBackground$total_count[genomicBackground$kmer == 'GATT']) / (2 * GENOMESIZE), 
          sum(genomicBackground$total_count[genomicBackground$kmer %in% c('GATT', nearGATT)]) / (2 * GENOMESIZE)), 
  'linetype' = c('a', 'b'), # c('dotted', 'dashed', 'solid'), 
  'name' = c('Optimal', 'Optimal+sub-optimal'), stringsAsFactors = FALSE)

sampleTbl <- figSampleTbl[figSampleTbl[, 'figure'] %in% figSet, ]

for (figName in proportionFigSet){
  subsetSampleTbl <- figSampleTbl[figSampleTbl[, 'figure'] == figName, ]
  sampleOfInterest <- subsetSampleTbl[, 'sample_name_internal']
  print(sampleOfInterest)

  # input PAM occurrences
  allPamCountTbl <- NULL
  groupRes <- NULL
  for (i in 1:length(sampleOfInterest)){
    s <- sampleOfInterest[i]
    print(s)
    pathFilenameSet <- list.files(
      path = inputPath, pattern = '.pam_proportion.csv', full.name = T)
    filename <- pathFilenameSet[grepl(pattern = s, 
                                      x = pathFilenameSet, perl = T)]
    if (length(filename) > 1){
      stop('More than one PAM proportion table were found: ' %strcat% s)
    }
    pamCountTbl <- read.table(filename, header = T, sep = ',', quote = '"')
    pamCount <- pamCountTbl[, 'Event_Count']
    names(pamCount) <- pamCountTbl[, 'PAM']
    orderedPamCount <- pamCount[kmerOrder]
    names(orderedPamCount) <- NULL
    
    pamCountDf <- data.frame(
      sampleInfo = s, 
      pam = kmerOrder, 
      pamCount = orderedPamCount, 
      pamCountFraction = 
        orderedPamCount / sum(orderedPamCount, na.rm = TRUE)
    )
    if (is.null(allPamCountTbl)){
      allPamCountTbl <- cbind('sampleInfo' = s, groupTbl)
    } else {
      allPamCountTbl <- rbind(allPamCountTbl, cbind('sampleInfo' = s, groupTbl))
    }
    # print(head(allPamCountTbl))
    # print(dim(allPamCountTbl))

    # # mark PAM group
    fieldSet <- c('pam', 'pamCount', 'pamCountFraction')
    exact <- cbind(group = 'optimal', pamCountDf[pamCountDf$pam == 'GATT', fieldSet], stringsAsFactors = FALSE)
    # print(exact)
    near <- cbind(group = 'sub-optimal', pamCountDf[pamCountDf$pam %in% nearGATT, fieldSet], stringsAsFactors = FALSE)
    near <- near[match(x = nearGATT, table = near[, 'pam']), ]
    # print(near)
    far <- cbind(group = 'remainder', pamCountDf[!(pamCountDf$pam %in% c(nearGATT, 'GATT')), fieldSet], stringsAsFactors = FALSE)
    far <- far[order(far[, 'pamCount'], decreasing = T), ]
    # print(far)
    # print(dim(far))
    # print(head(far))
    
    groupTbl <- rbind(exact, near, far)
    # print(head(groupTbl))
    usage <- c(
               sum(exact[, 'pamCountFraction'], na.rm = TRUE), 
               sum(near[, 'pamCountFraction'], na.rm = TRUE), 
               sum(far[, 'pamCountFraction'], na.rm = TRUE))
    # print(usage)
    res <- data.frame(
      'sample' = s, 
      'group' = c('exact', 'near', 'far'), 
      'pamCountFraction' = usage, 
      stringsAsFactors = FALSE
    )
    print(dim(res))
    print(res)
    print(sum(res[, 'pamCountFraction']))
    if(is.null(groupRes)){
      groupRes <- res
    } else {
      groupRes <- rbind(groupRes, res)
    }
  }
  
  #### plot bar chart
  plotDataStacked <- groupRes
  
  plotDataStacked[, 'group'] <- factor(x = plotDataStacked[, 'group'], levels = rev(c('exact', 'near', 'far')))
  plotDataStacked[, 'sample'] <- factor(x = plotDataStacked[, 'sample'], levels = sampleOfInterest)
  print(str(plotDataStacked))
  pStacked <- ggplot(data = plotDataStacked) +
    geom_bar(aes(x = sample, y = pamCountFraction,
                 fill = group), stat = 'identity', width = 0.75) +
    geom_hline(data = rulers, aes(yintercept = y, linetype = linetype), 
               size = 1, color = 'black') + 
    scale_x_discrete(name = 'Sample ID', limits = rev(sampleOfInterest)) + 
    scale_y_continuous(name = '% of PAM insertion events', 
                       limits = c(0, 1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1), minor_breaks = seq(0.05, 0.95, 0.05),
                       labels = paste(as.character(c(0, 0.25, 0.5, 0.75, 1) * 100), '%', sep = ''), 
                       expand = c(0, 0)) +
    scale_fill_manual(
      name = '4-mer groups: ',
      breaks = c('exact', 'near', 'far'), 
      labels = c('exact' = 'Optimal\n(GATT)', 
                 'near' = paste('Sub-optimal\n(', paste(nearGATT, collapse = ','), ')', sep = ''),
                 'far' = 'Non-optimal\n(Remaining 251 4-mers)'),
      values = c('#E69F00', '#CC79A7', '#009E73')) +

    scale_linetype_manual(
      name = '% of 4-mer on genome:', 
      values = c('a' = 'dotted', 'b' = 'dashed'), 
      breaks = c('a', 'b'), labels = rulers[rulers[, 'linetype'] %in% c('a', 'b'), 'name']) + 
    coord_flip() +
    labs(title = figName) +
    theme_bw() +
    theme(
      text = element_text(color = 'black', size = 20 * 0.5 * .pt), 
      line = element_line(color = 'black', size = 1),
      plot.background = element_blank(),
      plot.margin = margin(2, 2, 2, 2, 'cm'),
      panel.background = element_rect(color = 'black', size = 1),
      panel.grid.major = element_line(color = 'grey'),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_line(color = 'grey', linetype = 'dashed'),
      legend.position = 'bottom',
      legend.box = 'vertical', 
      legend.box.just = 'left', 
      legend.title = element_text(margin = margin(0.7, 0.1, 0.1, 0.1, 'cm')),
      axis.title.x = element_text(margin = margin(0.2, 0.2, 0.2, 0.2, 'cm')),
      axis.title.y = element_text(margin = margin(0.2, 0.2, 0.2, 0.2, 'cm')), 
      axis.text.x = element_text(color = 'black'),
      axis.text.y = element_text(color = 'black')
    ) +
    guides(
      fill = guide_legend(label.position = 'bottom', title.position = 'top', title.hjust = 0, reverse = FALSE), 
      linetype = guide_legend(label.position = 'right', title.position = 'top', title.hjust = 0, reverse = FALSE, direction = 'vertical', keywidth = unit(4, 'cm'), override.aes = list(size = 2)))
  pStackedFixed <- set_panel_size(pStacked, 
                                  height = unit(13, 'cm'), width = unit(5, 'cm')
  )
  
  ggsave(pStacked, 
         filename = concat('proportion/', figName,'.proportion.pdf'),
         width = unit(13, 'cm'), height = unit(13, 'cm'), device = cairo_pdf
  )
  
  #### save data
  samplePlotId <- 1:length(sampleOfInterest)
  names(samplePlotId) <- sampleOfInterest
  # print(samplePlotId)
  
  outTbl <- cbind(
    data.frame('fig' = figName, 
               'figRow' = samplePlotId[allPamCountTbl[, 'sampleInfo']], 
               'batch_sample' = allPamCountTbl[, 'sampleInfo']
               ), 
    allPamCountTbl[, c('pam', 'pamCount', 'pamCountFraction')], 
    data.frame('group' = pamGroup[allPamCountTbl[, 'pam']]))
  print(head(outTbl))
  colnames(outTbl) <- c('Figure', 'Row_in_figure', 'Batch_sample', 
                        'PAM', 'PAM_count', 'PAM_proportion', 'PAM_group')
  print(head(outTbl))
  write.table(outTbl, 
              file = concat('proportion/', figName,'.proportion.csv'),
              sep = ',', row.names = FALSE, col.names = TRUE, quote = TRUE
  )
  
}
