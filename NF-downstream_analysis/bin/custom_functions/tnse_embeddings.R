# function which returns tsne embeddings from antler object

tsne_embeddings <- function(curr_m, gms, seed=1, pca=F, perplexity=12, eta=200){
  set.seed(seed)
  data.logscaled =  t(scale(t(log(curr_m$getReadcounts(data_status='Normalized')[unlist(gms),]+1)), center=TRUE, scale=TRUE))
  tsne_xy = Rtsne::Rtsne(t(data.logscaled), pca=pca, perplexity=perplexity, eta=eta, max_iter=2000, verbose=T)$Y
  rownames(tsne_xy) <- colnames(data.logscaled)
  return(tsne_xy)
}