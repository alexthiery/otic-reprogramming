# BiocManager::install("ReadqPCR")
# BiocManager::install("NormqPCR")
# install.packages('rstatix')

data_dir = './data/'
reformatted_dir = './data/reformatted/'
dir.create(reformatted_dir, recursive = TRUE)

library(tidyverse)
library(ReadqPCR)
library(NormqPCR)
library(rstatix)
library(RColorBrewer)

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


# calculate delta ct
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



