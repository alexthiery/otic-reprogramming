data_dir = './data/'
reformatted_dir = './data/reformatted/'
dir.create(reformatted_dir, recursive = TRUE)

library(tidyverse)
library(ReadqPCR)
library(NormqPCR)
library(rstatix)
library(RColorBrewer)
library(scales)

# read in all csv files and re-format
qPCR_data <- lapply(list.files(data_dir, pattern = '*.csv', full.names = TRUE), function(x){
  as.data.frame(t(read.csv(x, row.names = 1))) %>%
    rownames_to_column() %>%
    # remove rows with NA
    drop_na() %>%
    pivot_longer(cols = which(colnames(.) != 'rowname')) %>%
    setNames(c('Sample', 'Detector', 'Cq'))
})
  names(qPCR_data) <- sub('*.csv', '', list.files(data_dir, pattern = '*.csv'))

# write re-formatted data to re-import with read.qPCR
lapply(names(qPCR_data), function(x) {write_delim(qPCR_data[[x]], paste0(reformatted_dir, x, '.txt'))} )

# read in new data with read.qPCR package
qPCR_data <- lapply(list.files(reformatted_dir, pattern = '*.txt', full.names = TRUE), function(x) read.qPCR(x))
names(qPCR_data) <- sub('*.txt', '', list.files(reformatted_dir, pattern = '*.txt'))

# calculate delta ct - if GAPDH then use both RPLP1 and GAPDH, otherwise just RPLP1
deltact <- lapply(qPCR_data, function(x) if('GAPDH' %in% featureNames(x)){deltaCq(qPCRBatch = x, hkgs=c('RPLP1', 'GAPDH'), calc="geom")} else {deltaCq(qPCRBatch = x, hkgs='RPLP1')})


# automatically generate comparisons matrix and add column/rownames
contM <- list()
for(x in names(qPCR_data)){
  df <- cbind(ifelse(grepl('Control', sampleNames(qPCR_data[[x]])), 1, 0), ifelse(grepl('Control', sampleNames(qPCR_data[[x]])), 0, 1))
  colnames(df) <-  c('Control', unique(sub('_.*', '', sampleNames(qPCR_data[[x]])))[unique(sub('_.*', '', sampleNames(qPCR_data[[x]]))) != 'Control'])
  rownames(df) <- sampleNames(qPCR_data[[x]])
  contM[[x]] <- df
}


# calculate ddct
ddct <- lapply(names(qPCR_data), function(x) {
  if('GAPDH' %in% featureNames(qPCR_data[[x]])){
    deltaDeltaCq(qPCRBatch = qPCR_data[[x]], hkgs=c('RPLP1', 'GAPDH'), contrastM = contM[[x]], case=colnames(contM[[x]])[colnames(contM[[x]]) != 'Control'], control='Control', hkgCalc="geom")
  } else {
    deltaDeltaCq(qPCRBatch = qPCR_data[[x]], hkgs='RPLP1', contrastM = contM[[x]], case=colnames(contM[[x]])[colnames(contM[[x]]) != 'Control'], control='Control')
    }
  }
)

names(ddct) <- names(qPCR_data)


# specify genes and associated colours of interest
library(viridis)

colors = brewer.pal(n = 5, name = "Pastel1")

colors = viridis(5)
names(colors) = c('LMX1A', 'FOXG1', 'PAX2', 'SOHO1', 'ZBTB16')

for(x in names(ddct)){
  # convert numeric factor values to numeric, and character factor values to character
  ddct[[x]][,2:8] <- lapply(ddct[[x]][,2:8], function(y) as.numeric(as.character(y)))
  ddct[[x]][,1] <- as.character(ddct[[x]][,1])
  
  # remove genes for plotting which are not of interest
  ddct[[x]] <- ddct[[x]][ddct[[x]][, 'ID'] %in% names(colors),]
}



# t-tests on delta ct
significance <- lapply(deltact, function(x){
  as.data.frame(t(exprs(x))) %>%
    rownames_to_column() %>%
    mutate(rowname = sub('_.*', '', rowname)) %>%
    gather(Gene, dct, 2:length(.)) %>%
    group_by(Gene) %>%
    # remove genes which are not to be tested
    filter(Gene %in% names(colors)) %>%
    t_test(dct ~ rowname, p.adjust.method = 'none') %>%
    mutate(padj = p.adjust(p, method = 'BH')) %>%
    mutate(p.adj.signif = ifelse(padj < 0.01, '**', ifelse(padj < 0.05, '*', 'ns')))
})

plot_data <- list()
# add significance to ddct dataframe
for(x in names(ddct)){
  plot_data[[x]] <- ddct[[x]][,c(1, 6:8)]
  plot_data[[x]] <- cbind(plot_data[[x]], sample = x)
  plot_data[[x]][, 'sig'] <- significance[[x]] %>%
    filter(Gene %in% ddct[[x]][, 'ID']) %>%
    arrange(ddct[[x]][, 'ID']) %>%
    pull(p.adj.signif)
}

# combine all experiments into a single data frame for plotting
plot_data <- reduce(plot_data, rbind)

# change experiment order for facet wrap order
plot_data$sample <- sub('_MO', ' aON', plot_data$sample)
plot_data$sample <- factor(plot_data$sample, levels = c('SOX8 aON', 'PAX2 aON', 'LMX1A aON', 'ZBTB16 aON'))


png('./qPCR_barplot.png', width = 18, height = 12, res = 200, units = 'cm')
ggplot(plot_data, aes(y = `2^-ddCt`, x = ID, fill = ID, label = ifelse(sig != 'ns', sig, ''))) +
  # hack - conditionally plot max min error bars and then hide symmetrical portion behind bar - doesnt affect results/plot, purely done for aesthetics
  geom_errorbar(aes(ymin = ifelse(`2^-ddCt` < 1, `2^-ddCt.min`, 1.1), ymax = ifelse(`2^-ddCt` > 1, `2^-ddCt.max`, 0.9)), width = 0.5) +
  geom_bar(stat='identity') +
  geom_text(aes(y = ifelse(`2^-ddCt` < 1, `2^-ddCt.min`, `2^-ddCt.max`)), nudge_y = ifelse(plot_data$`2^-ddCt` < 1, -0.2, 0.2) , size = 7) +
  # original error bars plotted on top of the bars
  #geom_errorbar(aes(ymin = `2^-ddCt.min`, ymax = `2^-ddCt.max`), width = 0.5) +
  scale_fill_manual(values = colors[plot_data$ID]) +
  scale_y_continuous(trans = log2_trans()) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank()) +
  facet_wrap(~sample, ncol = 4, strip.position = "bottom") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=10),
        strip.placement = "outside") +
  ylab(bquote(paste('fold change (', 2^-~Delta*Delta*Ct, ' )')))
graphics.off()


