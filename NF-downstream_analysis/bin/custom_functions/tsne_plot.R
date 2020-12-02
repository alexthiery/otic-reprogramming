# function to plot tsne from antler object
tsne_plot <- function(curr_m, gms, colour_by, colours, seed=1, pca=FALSE, perplexity=12, eta=200, main = NULL){
  set.seed(seed)
  data.logscaled =  t(scale(t(log(curr_m$getReadcounts(data_status='Normalized')[unlist(gms),]+1)), center=TRUE, scale=TRUE))
  tsne_xy = Rtsne::Rtsne(t(data.logscaled), pca=pca, perplexity=perplexity, eta=eta, max_iter=2000, verbose=T)$Y
  tsne_xy = cbind(as.data.frame(tsne_xy), cols = as.factor(m2$cellClusters$Mansel$cell_ids))

  if(length(unique(colour_by)) > length(unique(colours)) & length(unique(colours)) == 2){
    cat('\nPlotting with continuous colour scale\n\n')
    
    plot_data = cbind(as.data.frame(tsne_xy), colour_by = as.numeric(colour_by))
    plot = ggplot(plot_data, aes(x=V1, y=V2, color=colour_by)) +
      geom_point() +
      scale_colour_gradient(low = colours[1], high = colours[2])
  } else {
    cat('\nPlotting with discrete colour scale\n\n')
    
    plot_data = cbind(as.data.frame(tsne_xy), colour_by = as.factor(colour_by))
    plot = ggplot(plot_data, aes(x=V1, y=V2, color=colour_by)) +
      geom_point() +
      scale_colour_manual(values = colours)
  }
  return(plot)
}
