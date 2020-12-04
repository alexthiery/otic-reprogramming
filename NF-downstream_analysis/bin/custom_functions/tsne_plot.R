# function to plot tsne from antler object
tsne_plot <- function(curr_m, gms, colour_by, colours, seed=1, pca=FALSE, perplexity=12, eta=200, main = NULL){
  set.seed(seed)
  data.logscaled =  t(scale(t(log(curr_m$getReadcounts(data_status='Normalized')[unlist(gms),]+1)), center=TRUE, scale=TRUE))
  tsne_xy = Rtsne::Rtsne(t(data.logscaled), pca=pca, perplexity=perplexity, eta=eta, max_iter=2000, verbose=T)$Y

  if(length(unique(colour_by)) > length(unique(colours)) & length(unique(colours)) == 2){
    cat('\nPlotting with continuous colour scale\n\n')
    
    plot_data = data.frame('t-SNE 1' = tsne_xy[,1], 't-SNE 2' = tsne_xy[,2], 'colour_by' = as.numeric(colour_by), check.names = FALSE)
    plot = ggplot(plot_data, aes(x=`t-SNE 1`, y=`t-SNE 2`, color=colour_by)) +
      geom_point() +
      scale_colour_gradient(low = colours[1], high = colours[2]) +
      theme_classic() +
      theme(legend.position = "none", axis.ticks=element_blank(), axis.text = element_blank())
  } else {
    cat('\nPlotting with discrete colour scale\n\n')
    
    plot_data = data.frame('t-SNE 1' = tsne_xy[,1], 't-SNE 2' = tsne_xy[,2], 'colour_by' = factor(colour_by), check.names = FALSE)
    plot = ggplot(plot_data, aes(x=`t-SNE 1`, y=`t-SNE 2`, color=colour_by)) +
      geom_point() +
      scale_colour_manual(values = colours) +
      theme_classic() +
      theme(legend.position = "none", axis.ticks=element_blank(), axis.text = element_blank())
  }
  return(plot)
}
