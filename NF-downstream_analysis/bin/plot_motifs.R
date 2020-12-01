#!/usr/bin/env Rscript

library(ggseqlogo)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(extrafont)
font_import(prompt = FALSE)

output_path = "./output/"
dir.create(output_path, recursive = T)


######## read in data
# read in logo data
motif_logos = list()
for(i in 1:20){
  motif_logos[[paste(i)]] <- t(read.delim(paste0('./ATAC_motif_enrichment/knownResults/known', i, '.motif'))[1:4])
  rownames(motif_logos[[paste(i)]]) = c('A', 'C', 'G', 'T')
}

# read in motif info
motif_meta = read.delim(paste0('./ATAC_motif_enrichment/knownResults.txt'))[1:20,c(1,3)]

# strip name
motif_meta[,1] <- sub("\\(.*", "", motif_meta[,1])


####### prepare grobs
# gene names
gene = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
  text(x = 0.5, y = 0.5, "Gene", cex = 15, col = "black", font=2)))

motif_names <- lapply(motif_meta[,1], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                 text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})


# motifs
motif = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                 text(x = 0.5, y = 0.5, "Motif", cex = 15, col = "black", font=2)))

motif_logos = lapply(motif_logos, function(x) {ggseqlogo(x, method = 'prob') + theme_void() + theme(plot.margin = unit(c(2,0,2,0), "cm"))})


# pvalues
pval = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                      text(x = 0.5, y = 0.5, "p-value", cex = 15, col = "black", font=2)))

motif_pval <- lapply(motif_meta[,2], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                             text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

######## plot grobs
png(paste0(output_path, 'top20_motifs.png'), width = 180, height = 300, family = 'Arial', units = 'cm', res = 400)
grid.arrange(grobs=c(gene, motif_names, motif, motif_logos, pval, motif_pval), ncol=3, widths = c(1, 3, 1), as.table=FALSE)
graphics.off()




########  plot selected motifs
motifs_of_interest <- c('Sox3', 'Sox2', 'Sox10', 'TEAD3', 'Six2', 'Six1', 'Sox9', 'AP-2alpha')
motifs_of_interest <- which(motif_meta$Motif.Name %in% motifs_of_interest)

png(paste0(output_path, 'selected_motifs.png'), width = 150, height = 150, family = 'Arial', units = 'cm', res = 400)
grid.arrange(grobs=c(gene, motif_names[motifs_of_interest], motif, motif_logos[motifs_of_interest], pval, motif_pval[motifs_of_interest]), ncol=3, widths = c(1, 4, 1), as.table=FALSE)
graphics.off()





