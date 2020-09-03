plot_tsne_coexpression <- function(curr_m, gms, gene1, gene2, basename = paste0(gene1, "_", gene2), seed=1, pca=F, perplexity=12, eta=200, plot_folder=output_path, main = NULL){
  genes <- c(gene1, gene2)
  a <- as.integer(100*log10(1+curr_m$getReadcounts(data_status='Normalized')[genes[1],]) / max(log10(1+curr_m$getReadcounts(data_status='Normalized')[genes[1],])))
  b <- as.integer(100*log10(1+curr_m$getReadcounts(data_status='Normalized')[genes[2],]) / max(log10(1+curr_m$getReadcounts(data_status='Normalized')[genes[2],])))
  
  dat <- data.frame(a, b, row.names = colnames(curr_m$getReadcounts(data_status='Normalized')))
  col.mat <- expand.grid(a=seq(0,100,by=1), b=seq(0,100,by=1))
  col.mat <- within(col.mat, mix <- rgb(green = a, red = b, blue = b, maxColorValue = 100))
  cell.cols <- unlist(apply(dat, 1, function(x){filter(col.mat, a == x[[1]] & b == x[[2]])[[3]]}))
  
  col.mat[,1:2] <- col.mat[,1:2]/100
  key.plot <- ggplot(col.mat, aes(x = col.mat[,1], y = col.mat[,2])) +
    xlab(genes[1]) +
    ylab(genes[2]) +
    geom_tile(aes(fill = mix)) +
    scale_fill_identity() +
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  
  set.seed(seed)
  data.logscaled =  t(scale(t(log(curr_m$getReadcounts(data_status='Normalized')[unlist(gms),]+1)), center=TRUE, scale=TRUE))
  tsne_xy = Rtsne::Rtsne(t(data.logscaled), pca=pca, perplexity=perplexity, eta=eta, max_iter=2000, verbose=T)$Y
  
  tsne_xy <- as.data.frame(tsne_xy, names(cell.cols))
  tsne.plot <- ggplot(tsne_xy, aes(x = tsne_xy[,1], y = tsne_xy[,2], color = rownames(tsne_xy))) +
    geom_point() +
    xlab("tSNE_1")+
    ylab("tSNE_2")+
    scale_color_manual(values=cell.cols)+
    scale_fill_manual(values=cell.cols)+
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  
  pdf(paste0(plot_folder, '/', basename, '_co-expression_TSNE.pdf'), height = 5, width = 7.5)
  gridExtra::grid.arrange(tsne.plot, key.plot, layout_matrix = rbind(c(1,1,2),
                                                                     c(1,1,NA)))
  graphics.off()
}