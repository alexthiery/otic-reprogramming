#!/usr/bin/env Rscript

output_path = './output'
dir.create(output_path)

library(gprofiler2)
library(dplyr)
library(ggplot2)

putative_enhancers <- read.delim(list.files('./', pattern = '.txt', full.names = TRUE))

# run functional enrichment analysis using GO:biological process and KEGG terms
fea_res <- gost(putative_enhancers$Entrez.ID, organism = 'ggallus', sources = c('GO:BP', 'KEGG'))

# generate URL for full results
# gost(putative_enhancers$Entrez.ID, organism = 'ggallus', sources = c('GO:BP', 'KEGG'), as_short_link = TRUE)

# select enriched terms of interest and generate bar plot
plot_data <- output$result %>%
  filter(term_id %in% c("GO:0022008", "GO:0007399", "GO:0030182", "GO:0007411", "GO:0060070", "KEGG:04330", "KEGG:04310")) %>%
  select(c(p_value, term_name, term_id)) %>%
  mutate(-log10(p_value)) %>%
  arrange(desc(`-log10(p_value)`)) %>%
  mutate(term_name = paste0(term_name, ' (', term_id, ")")) %>%
  mutate(term_name = factor(term_name, levels = term_name))

png(paste0(output_path, "functional_enrichment.png"), height = 20, width = 25, units = "cm", res = 400)
ggplot(plot_data, aes(x = term_name, y = -log10(p_value), label = term_name)) +
  geom_bar(stat='identity', width=0.5, fill='steelblue') +
  coord_flip() +
  geom_text(aes(y = 0), hjust = 'left', vjust = -2, size = 3.5) +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10))
graphics.off()