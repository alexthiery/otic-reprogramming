#!/usr/bin/env Rscript
library(ggplot2)
output_path = './output/'
dir.create(output_path)

# import ATAC peaks intersected with +K27Ac -K27me3
peaks <- read.delim(list.files('./', pattern=".txt", full.names = TRUE), sep = "\t")

# extract and simplify annotations for categorisation
peaks <- as.factor(sub(' .*', "", peaks$Annotation))

# order frequency
freq_data <- as.data.frame(prop.table(table(peaks))[order(prop.table(table(peaks)))])
colnames(freq_data) = c('peaks', 'Frequency')


# plot frequency plot of peak annotations
png(paste0(output_path, "peak_annotation_distribution.png"), height = 10, width = 10, units = "cm", res = 400)
ggplot(freq_data, aes(x = peaks, y = Frequency)) +
  geom_bar(stat='identity', fill='steelblue') +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5))
graphics.off()

