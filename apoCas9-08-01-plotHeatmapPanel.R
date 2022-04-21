rm(list = ls())
`%strcat%` <- function(x, y) paste(x, y, sep = '')
concat <- function(...){
  paste(list(...), collapse = '')
}

library('reshape2')
library('ggplot2')
library('grid')
library('egg')
library('cowplot')
source('apoCas9-08-00-funcLib.R', sep = '')
dataPath <- '~/data/apoCas9/plots/'
# # combine PWM to plot data with dominating bases
collapseDominatingBase <- function(pwmFormatted){
  domBases <- vector(mode = 'character', length = ncol(pwmFormatted))
  domPepwmentages <- vector(mode = 'numeric', length = ncol(pwmFormatted))
  alphabet <- rownames(pwmFormatted)
  for (i in 1:ncol(pwmFormatted)){
    domRow <- which.max(pwmFormatted[, i])
    domBases[i] <- alphabet[domRow]
    domPepwmentages[i] <- max(pwmFormatted[, i])
  }
  return(list(domBases, domPepwmentages))
}

readPwmFromFilenameSet <- function(filenameSet){
  dataDf <- data.frame(matrix(data = NA, nrow = 0, ncol = 4))
  colnames(dataDf) <- c('sample', 'position', 'dominatingBase', 'dominatingPercentage')
  for (filename in filenameSet){
    pwm <- read.table(file = filename, sep = ',', header = FALSE)
    rownames(pwm) <- c('A', 'C', 'G', 'T')
    dominating <- collapseDominatingBase(pwm)
    sampleId <- getSampleId(filename)
    singleSampleDf <- cbind.data.frame(c(sampleId),
                                       c(1:ncol(pwm)),
                                       dominating[[1]],
                                       dominating[[2]])
    colnames(singleSampleDf) <- c('sample', 'position', 'dominatingBase', 'dominatingPercentage')
    dataDf <- rbind(dataDf, singleSampleDf)
  }
  return(dataDf)
}

readAllNucPwmFromFilenameSet <- function(filenameSet){
  filenameSet <- sampleFullnameSet  
  dataDf <- data.frame(matrix(data = NA, nrow = 0, ncol = 4))
  for (filename in filenameSet){
    pwm <- read.table(file = filename, sep = ',', header = FALSE)
    pamPwm <- pwm[, basePlotRange]
    sampleId <- getSampleId(filename)
    singleSampleDf <- cbind.data.frame(c(sampleId), c('A', 'C', 'G', 'T'), pamPwm)
    colnames(singleSampleDf) <- c('sample', 'nucleotide', paste('pos', as.character(basePlotRange), sep = ''))
    dataDf <- rbind(dataDf, singleSampleDf)
  }
  return(dataDf)
}
formatHeatmapWithLegends <- function(hmFull, yesSpecialEntLegend = F, specialEntLegend = NULL){
  hmFixed <- set_panel_size(hmFull, width = unit(2 * length(basePlotRange), 'cm'), height = unit(2 * 0.75 * length(sampleFullnameSet), 'cm'))
  
  if (yesSpecialEntLegend){
    legendEnt <- specialEntLegend
  } else {
    legendEntPlotData <- data.frame(xmin = c(1:4) - 0.5, xmax = c(1:4) + 0.5, ymin = 0.5, ymax = 1.5, a = c(0, 0.33, 0.67, 1))
    legendEnt <- ggplot(data = legendEntPlotData) + 
      geom_rect(aes(alpha = a, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'black', color = 'black') + 
      scale_x_continuous(breaks = 1:nrow(legendEntPlotData), label = rev(as.character(legendEntPlotData[, 'a']))) +
      annotate(geom = 'text', x = c(1:4), y = 2, label = c('\u2A7E1', '0.67', '0.33', '\u2A7D0'), size = 20 * 0.5 * 0.8, color = 'black', hjust = 0.5) + 
      annotate(geom = 'text', x = 1.5, y = 3, label = 'Most random', size = 20 * 0.5 * 0.8, color = 'black', hjust = 1) + 
      annotate(geom = 'text', x = 3.5, y = 3, label = 'Most informative', size = 20 * 0.5 * 0.8, color = 'black', hjust = 0) + 
      scale_y_reverse(position = 'right') +
      scale_alpha(range = c(0, 1), guide = FALSE) + 
      labs(title = 'Normalized PAM Entropy') + 
      coord_cartesian(clip = 'off') + 
      theme_void() + 
      theme(
        text = element_text(size = 20 * 0.5 * .pt * 0.8),
        plot.title = element_text(hjust = 0),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, unit = 'cm'))
  }
  
  legendEntFixed <- set_panel_size(legendEnt, width = unit(legendKeyWidth * 2 * 4, 'cm'), height = unit(legendKeyWidth * 1 * 2, 'cm'))
  
  legendBasePlotData <- data.frame(xmin = rep(c(1:4) - 0.5, 4), xmax = rep(c(1:4) + 0.5, 4), ymin = rep(c(1:4) - 0.5, each = 4), ymax = rep(c(1:4) + 0.5, each = 4), a = rep(c(seq(0.25, 1, 0.25)), 4), b = factor(rep(c('A', 'C', 'G', 'T'), each = 4), levels = c('A', 'C', 'G', 'T')))
  basePercentageLabel <- paste(as.character(seq(0.25, 1, 0.25) * 100), '%', sep = '')
  legendBase <- ggplot(data = legendBasePlotData) +
    geom_rect(aes(alpha = a, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = b), color = 'grey52') + 
    scale_y_reverse(breaks = c(1:4), labels = c('A', 'C', 'G', 'T'), position = 'right') +
    scale_x_continuous(breaks = 1:4, label = basePercentageLabel) + 
    scale_alpha(range = c(0, 1), limits = c(0.25, 1), guide = FALSE) + 
    scale_fill_manual(values = hmPalette, guide = FALSE) + 
    labs(title = 'Dominant Base') + 
    theme_void() + 
    theme(
      text = element_text(size = 20 * 0.5 * .pt * 0.75),
      plot.title = element_text(hjust = 0),
      plot.margin = margin(0.25, 0.25, 0.25, 0.25, unit = 'cm'),
      axis.text.x = element_text(color = 'black', hjust = 0.5), 
      axis.text.y = element_text(color = 'black'))
  legendBaseFixed <- set_panel_size(legendBase, width = unit(legendKeyWidth * 4 * 2, 'cm'), height = unit(legendKeyWidth * 4 * 0.75, 'cm'))
  legendCol <- plot_grid(legendEntFixed, legendBaseFixed, nrow = 2, align = 'v')
  p <- plot_grid(hmFixed, legendCol, ncol = 2, rel_widths = c(1.5, 1))
  return(p)
}

saveHeatmap <- function(save_plot, save_format){
  save_suffix <- save_format
  if (save_format == 'pdf'){
    save_device <- cairo_pdf
  } else if (save_format == 'svg') {
    save_device <- cairo_ps
  } else if (save_format == 'png') {
    save_device <- png
    fullWidth <- fullWidth + 20
  }
  ggsave(
    plot = save_plot, 
    filename = concat('heatmap/manuscript.', figName, '.heatmap.csv'), 
    device = save_device, 
    height = fullHeight, width = fullWidth, unit = 'cm', bg = 'white')
}

# # pre-set heat map color palette
hmPalette <- c('A' = '#6ACC64', 'C' = '#4878D0', 
               'G' = '#F0D373', 'T' = '#D65F5F')

# # define plots of interest
basePlotRange <- c(51:60)

# # read sample list for each figure-panel
figSampleTbl <- read.table(file = 'ref/figure_sample_list.csv', sep = ',', colClasses = 'character', header = TRUE)

# # list all figures to be plotted
figureList <- unique(figSampleTbl[, 'figure'])
heatmapFigLabel <- getSubFigLabels(plotType = 'heatmap')

# # read normalized PAM entropy
normPamEntTbl <- read.table(file = dataPath %strcat% 'pamEntropy.norm.csv', 
                            sep = ',', colClasses = c(rep('character', 2), rep('numeric', 2)), quote = '"', header = TRUE)
# print(head(normPamEntTbl))
# print(table(normPamEntTbl[, 'batch']))
normPamEntTbl[, 'sampleInfo'] <- paste(normPamEntTbl[, 'batch'], normPamEntTbl[, 'sample'], sep = '-')

# # plot all heat maps
allFigDataTbl <- NULL
for (figName in figureList){
  print(figName)
  figTitle <- 'Fig.' %strcat% heatmapFigLabel[figName]
  subFigSampleTbl <- figSampleTbl[figSampleTbl[, 'figure'] == figName, ]
  print(subFigSampleTbl)
  
  sampleFullnameSet <- vector('character', nrow(subFigSampleTbl))
  for (i in 1:nrow(subFigSampleTbl)){
    # format sample names to the full path to PWM (csv) files
    # using customized functions from a sourced script
    sampleFullname <- getSampleFullName(getPwmDataPath(subFigSampleTbl[i, 'sample_name_internal']))
    sampleFullnameSet[i] <- sampleFullname
  }
  pwmFigRanged <- readPwmFromFilenameSet(filenameSet = sampleFullnameSet)
  
  pwmNucRanged <- readAllNucPwmFromFilenameSet(filenameSet = sampleFullnameSet)
  write.table(pwmNucRanged,
              file = concat('heatmap/manuscript.', figName, '.heatmap.pwm.csv'),
              sep = ',', quote = TRUE, row.names = FALSE)
 
  
  # # process data for plotting heat map
  hmPlotData <- pwmFigRanged
  hmPlotData[, 'sample'] <- factor(hmPlotData[, 'sample'], levels = unique(hmPlotData[, 'sample']))
  hmPlotData[, 'position'] <- factor(hmPlotData[, 'position'], levels = basePlotRange)
  print(str(hmPlotData))
  print(summary(hmPlotData[, 'dominatingPercentage']))
  
  # # annotate normalized PAM entropy
  subNormPamEntTbl <- normPamEntTbl[normPamEntTbl[, 'sampleInfo'] %in% hmPlotData[, 'sample'], ]
  print(subNormPamEntTbl)
  
  # # make data frame for PAM entropy and adaptation efficiencies data plotting as tiles
  annotationTbl <- subNormPamEntTbl
  annotationTbl[, 'sampleInfo'] <- factor(annotationTbl[, 'sampleInfo'], levels = levels(hmPlotData[, 'sample']))
  annotationTbl[, 'rect_ymax'] <- (nrow(annotationTbl) - as.numeric(annotationTbl[,'sampleInfo']) + 1) + 0.5
  annotationTbl[, 'rect_ymin'] <- (nrow(annotationTbl) - as.numeric(annotationTbl[,'sampleInfo']) + 1) - 0.5
  annotationTbl[, 'rect_ent_xmax'] <- length(basePlotRange) + 1 + 0.5
  annotationTbl[, 'rect_ent_xmin'] <- length(basePlotRange) + 1 - 0.5
  annotationTbl[, 'ent_label'] <- format(x = round(annotationTbl[, 'normed_entropy'], digits = 2), digits = 2)
  normEntScaledAlpha <- sapply(X = annotationTbl[, 'normed_entropy'], FUN = function(x) 1 - max(min(ifelse(is.na(x), 1, x), 1), 0))
  annotationTbl[, 'ent_scaled'] <- normEntScaledAlpha / 1 * 0.75 + 0.25 # adjust to use the same alpha scale as base percentages
  
  # # plot
  hm <- ggplot(data = hmPlotData, aes(x = position, y = sample)) +
    geom_tile(aes(fill = dominatingBase, alpha = dominatingPercentage), color = NA) +
    geom_rect(data = annotationTbl, 
              aes(ymin = rect_ymin, ymax = rect_ymax, xmax = rect_ent_xmax, xmin = rect_ent_xmin, 
                  alpha = ent_scaled), 
              fill = 'black', color = 'black', size = 1, inherit.aes = FALSE) + 
    geom_text(data = annotationTbl, aes(y = sampleInfo, label = ent_label), 
              x = length(basePlotRange) + 1 + 0.5 + 0.05, color = 'black', size = 20 * 0.5, hjust = 0) + 
    scale_x_discrete(labels = as.character(c(1:length(basePlotRange))), expand = c(0, 0)) + 
    scale_y_discrete(limits = rev(levels(hmPlotData[, 'sample'])), expand = c(0, 0)) +
    scale_fill_manual(name = 'Dominant Base',
                      values = hmPalette,
                      breaks = c('A', 'C', 'G', 'T'),
                      labels = c('A', 'C', 'G', 'T'), guide = FALSE) +
    scale_alpha_continuous(name = 'Normalized PAM entropy', 
                           breaks = c(0.25, 0.5, 0.75, 1),
                           labels = c('0.00', '0.33', '0.67', '1.00'),
                           range = c(0, 1), limits = c(0.25, 1), na.value = 0.25, 
                           guide = FALSE
                           ) +
    coord_fixed(ratio = 0.75, # ratio = y / x 
                clip = 'off') + 
    theme(
      text = element_text(size = 20 * 0.5 * .pt), 
      line = element_line(size = 2), 
      plot.margin = margin(0.25, 0.25, 0.25, 0.25, unit = 'cm'),
      panel.background = element_rect(fill = NA, color = 'black', size = 1, linetype = 'solid'),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(color = 'black', hjust = 0.5, margin = margin(t = 0.2, unit = 'cm')),
      axis.text.y = element_text(color = 'black', margin = margin(t = 0, r = 0.1, b = 0, l = 0, unit = 'cm')),
      legend.position = 'right', 
      legend.text = element_text(margin = margin(t = 0.1, r = 0, b = 0.1, l = 0, unit = 'cm'), vjust = 0.5), 
      legend.key = element_rect(fill = NA),
      legend.key.height = unit(1, 'cm'), 
      legend.key.width = unit(1, 'cm')
    ) +
    labs(title = figTitle)
  fullHeight <- (2 * 0.75 * max(length(sampleFullnameSet), 6)+ 1.5 + 2.5) 
  hmFormat <- formatHeatmapWithLegends(hmFull = hm)
  saveHeatmap(save_plot = hmFormat, save_format = 'pdf')
  # saveHeatmap(save_plot = hmFormat, save_format = 'svg')
  # saveHeatmap(save_plot = hmFormat, save_format = 'png')

  # # # set-up supplementary table
  # # organize suppl table contents
  print(subNormPamEntTbl)
  
  outTbl <- NULL
  for (i in 1:length(levels(hmPlotData[, 'sample']))){
    s <- levels(hmPlotData[, 'sample'])[i]
    print(s)
    samplePwm <- pwmFigRanged[pwmFigRanged[, 'sample'] == s, ]
    samplePwm <- samplePwm[order(samplePwm[, 'position'], decreasing = FALSE), ]
    relPos <- samplePwm[, 'position'] - min(samplePwm[, 'position']) + 1
    cname <- c()
    cvalue <- c()
    for (rpos in relPos){
      cname <- c(cname, paste('pos', as.character(rpos), c('base', 'percentage'), sep = '_'))
      cvalue <- c(cvalue, as.character(samplePwm[rpos, c('dominatingBase')]), samplePwm[rpos, c('dominatingPercentage')])
    }
    names(cvalue) <- cname
    sampleOutTbl <- cbind.data.frame('fig' = figTitle, 'fig_sample_row' = i, 
                          t(as.data.frame(cvalue)))
    rownames(sampleOutTbl) <- NULL
    sampleOutTbl['normed_pam_ent'] <- subNormPamEntTbl[subNormPamEntTbl[, 'sampleInfo'] == s, 'normed_entropy']
    if (is.null(outTbl)){
      outTbl <- sampleOutTbl
    } else {
      outTbl <- rbind(outTbl, sampleOutTbl)
    }
  }
  print(outTbl)
  if (is.null(allFigDataTbl)){
    allFigDataTbl <- outTbl
    outColnames <- colnames(outTbl)
    outColnames[outColnames == 'fig'] <- 'Figure'
    outColnames[outColnames == 'fig_sample_row'] <- 'Row.in.Figure'
    outColnames <- gsub(x = outColnames, pattern = 'pos_', replacement = 'Col.', fixed = TRUE) # Column
    outColnames <- gsub(x = outColnames, pattern = 'base', replacement = 'Nuc.Base.', fixed = TRUE) # Nucleotide Base
    outColnames <- gsub(x = outColnames, pattern = 'percentage', replacement = 'Rel.Freq.', fixed = TRUE) # Relative Frequency
    outColnames[outColnames == 'normed_pam_ent'] <- 'Norm.PAM.Ent.' # normalized PAM Entropy
  } else {
    allFigDataTbl <- rbind(allFigDataTbl, outTbl)
  }
  write.table(as.data.frame(outTbl),
              file = concat('heatmap/manuscript', '.', figName, '.heatmap.csv'),
              sep = ',', quote = TRUE, row.names = FALSE,
              col.names = outColnames)
}

write.table(allFigDataTbl,
            file = 'heatmap/manuscript.heatmap_data.csv',
            sep = ',', quote = TRUE, row.names = FALSE,
            col.names = outColnames)

