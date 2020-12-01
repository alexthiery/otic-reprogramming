# function to plot tsne from antler object

tsne_plot <- function(curr_m, gms, basename, cols=clust.colors[curr_m$cellClusters[['hclust']]$cell_ids], seed=1, pca=FALSE, perplexity=12, eta=200, plot_folder=output_path, main = NULL, ...){
  set.seed(seed)
  data.logscaled =  t(scale(t(log(curr_m$getReadcounts(data_status='Normalized')[unlist(gms),]+1)), center=TRUE, scale=TRUE))
  tsne_xy = Rtsne::Rtsne(t(data.logscaled), pca=pca, perplexity=perplexity, eta=eta, max_iter=2000, verbose=T)$Y
  tsne_xy = cbind(as.data.frame(tsne_xy), cols = as.factor(m2$cellClusters$Mansel$cell_ids))
  png(paste0(plot_folder, basename, '_TSNE.png'), ...)
  ggplot(tsne_xy, aes(x=V1, y=V2, color=cols)) +
    geom_point() +
    scale_colour_manual(values = clust.colors) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "none")
  dev.off()
}