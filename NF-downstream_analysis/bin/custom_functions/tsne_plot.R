# function to plot tsne from antler object

tsne_plot <- function(curr_m, gms, basename, cols=clust.colors[curr_m$cellClusters[['hclust']]$cell_ids], seed=1, pca=FALSE, perplexity=12, eta=200, plot_folder=output_path, main = NULL){
  set.seed(seed)
  data.logscaled =  t(scale(t(log(curr_m$getReadcounts(data_status='Normalized')[unlist(gms),]+1)), center=TRUE, scale=TRUE))
  tsne_xy = Rtsne::Rtsne(t(data.logscaled), pca=pca, perplexity=perplexity, eta=eta, max_iter=2000, verbose=T)$Y
  pdf(paste0(plot_folder, '/', basename, '_TSNE.pdf'))
  plot(tsne_xy, col=cols, pch=16, main = main)
  dev.off()
}