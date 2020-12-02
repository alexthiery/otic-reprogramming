#!/usr/bin/env Rscript

library(gprofiler2)
library(dplyr)
library(ggplot2)
library(extrafont)

output_path = './output/'
dir.create(output_path, recursive = T)

putative_enhancers <- read.delim(list.files(pattern = '*txt', full.names = TRUE))

# run functional enrichment analysis using GO:biological process and KEGG terms
fea_res <- gost(putative_enhancers$Entrez.ID, organism = 'ggallus', sources = c('GO:BP', 'KEGG'))

# generate URL for full results
# gost(putative_enhancers$Entrez.ID, organism = 'ggallus', sources = c('GO:BP', 'KEGG'), as_short_link = TRUE)

go_terms <- c("GO:0007399", "KEGG:04310", "GO:0048839", "GO:0050767", "GO:0043408", "KEGG:04330")

# select enriched terms of interest and generate bar plot
plot_data <- fea_res$result %>%
  filter(term_id %in% go_terms) %>%
  select(c(p_value, term_name, term_id)) %>%
  mutate(-log10(p_value)) %>%
  arrange(desc(`-log10(p_value)`)) %>%
  mutate(term_name = paste0(term_name, ' (', term_id, ")")) %>%
  mutate(term_name = factor(term_name, levels = term_name))

png(paste0(output_path, "functional_enrichment.png"), height = 10, width = 15, family = 'Arial', units = 'cm', res = 400)
ggplot(plot_data, aes(x = term_name, y = -log10(p_value), label = term_name)) +
  geom_bar(stat='identity', width=0.5, fill='steelblue') +
  coord_flip() +
  geom_text(aes(y = 0), hjust = 'left', vjust = -2, size = 3.5) +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20)) +
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
graphics.off()
