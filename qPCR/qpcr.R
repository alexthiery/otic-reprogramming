


library(tidyverse)
# install.packages('pcr')
# library(pcr)

# # calculate ddct using PCR package
# temp <- as.data.frame(t(read.csv('~/Desktop/qPCR/data/SOX8MO.csv', row.names = 1))) %>%
#   # remove rows with NA
#   drop_na()
# 
# qPCR_data <- lapply(list.files('~/Desktop/qPCR/data/', full.names = T), function(x) {
#   as.data.frame(t(read.csv(x, row.names = 1))) %>%
#     # remove rows with NA
#     drop_na()
#   })
# 
# names(qPCR_data) <- sub('_.*', '', list.files('~/Desktop/qPCR/data/'))
# 
# group_var <- lapply(qPCR_data, function(x) sub('_.*', '', rownames(x)))
# 
# res <- lapply(names(qPCR_data), function(x){
#   pcr_analyze(qPCR_data[[x]],
#               group_var = group_var[[x]],
#               reference_gene = 'RPLP1',
#               reference_group = 'Control')
#   
# })

# BiocManager::install("ReadqPCR")
# BiocManager::install("NormqPCR")

library(ReadqPCR)
library(NormqPCR)
library(rstatix)
library(RColorBrewer)

# 
# vignette("ReadqPCR")
# vignette("NormqPCR")


qPCR_data <- as.data.frame(t(read.csv('~/Desktop/qPCR/data/SOX8_MO.csv', row.names = 1))) %>%
  rownames_to_column() %>%
  # remove rows with NA
  drop_na() %>%
  pivot_longer(cols = which(colnames(.) != 'rowname'))



colnames(qPCR_data) = c('Sample', 'Detector', 'Cq')

write_delim(qPCR_data, '~/Desktop/qPCR/temp.txt')

qPCR_data <- read.qPCR('~/Desktop/qPCR/temp.txt')

deltact <- deltaCt(qPCRBatch = qPCR_data, hkgs=c('RPLP1', 'GAPDH'), calc="geom")



# automatically generate comparisons matrix
contM <- cbind(ifelse(grepl('Control', sampleNames(qPCR_data)), 1, 0), ifelse(grepl('Control', sampleNames(qPCR_data)), 0, 1))
# add column names regardless of test sample name
colnames(contM) <- c('Control', unique(sub('_.*', '', sampleNames(qPCR_data)))[unique(sub('_.*', '', sampleNames(qPCR_data))) != 'Control'])
rownames(contM) <- sampleNames(qPCR_data)

ddct <- deltaDeltaCq(qPCRBatch = qPCR_data, hkgs=c('GAPDH', 'RPLP1'),
                     contrastM=contM, case=colnames(contM)[colnames(contM) != 'Control'],
                     control="Control", hkgCalc = 'geom')

# convert dataframe to numeric
ddct[,2:8] <- lapply(ddct[,2:8], function(x) as.numeric(as.character(x)))
ddct[,1] <- as.character(ddct[,1])

# if value is below 1 then -1/value
# ddct[,6:8] <- mapply(function(x){ifelse(x < 1, -1/x, x)}, ddct[,6:8])

# specify genes and associated colours of interest
colors = brewer.pal(n = 5, name = "Dark2")
names(colors) = c('LMX1A', 'FOXG1', 'PAX2', 'SOHO1', 'ZBTB16')


# remove genes for plotting which are not of interest
ddct <- ddct[ddct$ID %in% names(colors),]



# t-tests on delta ct
significance <- as.data.frame(t(exprs(deltact))) %>%
  rownames_to_column() %>%
  mutate(rowname = sub('_.*', '', rowname)) %>%
  gather(Gene, dct, 2:length(.)) %>%
  group_by(Gene) %>%
  pairwise_t_test(dct ~ rowname, p.adjust.method = 'bonferroni')


# add significance to ddct dataframe
ddct$sig <- significance %>%
  filter(Gene %in% ddct$ID) %>%
  arrange(ddct$ID) %>%
  pull(p.adj.signif)



ggplot(ddct, aes(y = `2^-ddCt`, x = ID, fill = ID, label = ifelse(sig != 'ns', sig, ''))) +
  geom_bar(stat='identity') +
  geom_text(aes(y = `2^-ddCt.max` + 0.1), size = 10) +
  geom_errorbar(aes(ymin = `2^-ddCt.min`, ymax = `2^-ddCt.max`), width = 0.5) +
  scale_fill_manual(values = colors[ddct$ID]) +
  theme_classic()












colors['PAX2']


# manual
temp <- t(as.data.frame(t(read.csv('~/Desktop/qPCR/data/SOX8_MO.csv', row.names = 1))) %>% drop_na())



# function for calculating geometric mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


# subtract housekeeping gene (delta ct)
# t(apply(temp, 1, function(x) x-temp['RPLP1',]))





# subtract geometric mean of two housekeeping genes (delta ct)
temp <- t(apply(temp, 1, function(x) x-apply(temp[c('GAPDH'),], 2, gm_mean)))


# anova on delta ct
as.data.frame(t(temp)) %>%
  rownames_to_column() %>%
  mutate(rowname = sub('_.*', '', rowname)) %>%
  gather(Gene, dct, 2:length(.)) %>%
  group_by(Gene) %>%
  pairwise_t_test(dct ~ rowname, p.adjust.method = 'bonferroni')


# calcultate average control delta ct
control_av <- as.data.frame(t(temp)) %>%
  rownames_to_column() %>%
  mutate(rowname = sub('_.*', '', rowname)) %>%
  group_by(rowname) %>%
  dplyr::summarize_all(mean) %>%
  filter(rowname == 'Control') %>%
  column_to_rownames(var = 'rowname') %>%
  as_vector()
  
# delta ct - average control delta ct == delta delta ct
ddct = apply(t(temp), 1, function(x) x - control_av)

# calculate average 2^-ddct
as.data.frame(t(2^-ddct)) %>%
  rownames_to_column() %>%
  mutate(rowname = sub('_.*', '', rowname)) %>%
  group_by(rowname) %>%
  dplyr::summarize_all(gm_mean)


# BiocManager::install("ddCt")
# library(ddCt)
# 
# # manually calculate ddct
# 
# dct <- as.data.frame(t(read.csv('~/Desktop/qPCR/data/SOX8MO.csv', row.names = 1))) %>%
#   # remove rows with NA
#   drop_na() %>%
#   rownames_to_column() %>%
#   # select GAPDH as reference
#   pivot_longer(cols = which(colnames(.) != c('rowname', 'GAPDH'))) %>%
#   # calculate delta CT
#   mutate(delta_ct = value - GAPDH) %>%
#   mutate(treatment = sub('_.*', '', rowname))
# 
# 
# av_dct <- dct %>%
#   filter(treatment == 'Control') %>%
#   group_by(name) %>%
#   dplyr::summarize(Mean = mean(delta_ct, na.rm=TRUE))
# 
# 
# ddct <- dct %>%
#   group_by(name) %>%
#   mutate(ddct = delta_ct - as.numeric(av_dct[av_dct$name %in% name, 2])) %>%
#   mutate(relative_expression = 2^-ddct)
# 
# 
# ddct %>%
#   group_by(name, treatment) %>%
#   dplyr::summarize(Mean = mean(relative_expression, na.rm=TRUE))


