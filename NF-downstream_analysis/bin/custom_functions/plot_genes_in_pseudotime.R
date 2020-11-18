pseudotime_multiplot <- function(data, otic_branch = 1, epi_branch = 2, basename = 'relative_monocle_exp', gene_list,
         out_path = output_path, width = 10, height = 7, separate_plots = FALSE){
    
    otic = my_plot_genes_in_pseudotime(data[gene_list, which(pData(data)$State != epi_branch)], color_by = "timepoint", relative_expr=FALSE)
    epi = my_plot_genes_in_pseudotime(data[gene_list, which(pData(data)$State != otic_branch)], color_by = "timepoint", relative_expr=FALSE)
    
    smooth_curves = rbind.data.frame(
        cbind(otic, "branch"="Otic"),
        cbind(epi, "branch"="Epibranchial")
    )
    
    if(separate_plots == TRUE){
        for(branch in c('Otic', 'Epibranchial')){
            if(branch == 'Otic'){
                dat = smooth_curves[smooth_curves$branch == "Otic",]
            }else{
                dat = smooth_curves[smooth_curves$branch == "Epibranchial",]
            }
            pdf(paste0(out_path, basename, '_', branch, '.pdf'), width=width, height=height)
            print(plot(ggplot(aes(Pseudotime, expression, color=feature_label), data = dat) +
                           geom_line(aes(x = Pseudotime, y = expectation),
                                     data = dat, size=1) +
                           scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                         labels = scales::trans_format("log10", scales::math_format(10^.x)))))
            graphics.off()
        }
    }else{
        pdf(paste0(out_path, basename, '.pdf'), width=width, height=height)
        print(plot(ggplot(aes(Pseudotime, expression, color=feature_label), data = smooth_curves) +
                       geom_line(aes(x = Pseudotime, y = expectation),
                                 data = smooth_curves, size=1) +
                       scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
                       facet_grid(~ branch)))
        graphics.off()
    }
}



my_plot_genes_in_pseudotime <- function (cds_subset, min_expr = NULL, cell_size = 0.75, nrow = NULL,
                                         ncol = 1, panel_order = NULL, color_by = "State", trend_formula = "~ sm.ns(Pseudotime, df=3)",
                                         label_by_short_name = TRUE, relative_expr = TRUE, vertical_jitter = NULL,
                                         horizontal_jitter = NULL)
{
    f_id <- NA
    Cell <- NA
    if (cds_subset@expressionFamily@vfamily %in% c("negbinomial",
                                                   "negbinomial.size")) {
        integer_expression <- TRUE
    }
    else {
        integer_expression <- FALSE
        relative_expr <- TRUE
    }
    if (integer_expression) {
        cds_exprs <- exprs(cds_subset)
        if (relative_expr) {
            if (is.null(sizeFactors(cds_subset))) {
                stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
            }
            cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
        }
        cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    }
    else {
        cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
    }
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    if (integer_expression) {
        cds_exprs$adjusted_expression <- cds_exprs$expression
    }
    else {
        cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    }
    if (label_by_short_name == TRUE) {
        if (is.null(cds_exprs$gene_short_name) == FALSE) {
            cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
            cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
        }
        else {
            cds_exprs$feature_label <- cds_exprs$f_id
        }
    }
    else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    cds_exprs$f_id <- as.character(cds_exprs$f_id)
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)
    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
    model_expectation <- genSmoothCurves(cds_subset, cores = 1,
                                         trend_formula = trend_formula, relative_expr = T, new_data = new_data)
    colnames(model_expectation) <- colnames(cds_subset)
    expectation <- plyr::ddply(cds_exprs, plyr::.(f_id, Cell), function(x) data.frame(expectation = model_expectation[x$f_id,
                                                                                                                      x$Cell]))
    cds_exprs <- merge(cds_exprs, expectation)
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
    if (is.null(panel_order) == FALSE) {
        cds_exprs$feature_label <- factor(cds_exprs$feature_label,
                                          levels = panel_order)
    }
    
    return(cds_exprs[, c('feature_label', 'timepoint', 'expression', "expectation", 'Pseudotime')])
    
    
    cds_exprs$timepoint <- factor(cds_exprs$timepoint, levels=sort(unique(cds_exprs$timepoint)))
    
    q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
    if (is.null(color_by) == FALSE) {
        q <- q + geom_point(aes_string(color = color_by), size = I(cell_size),
                            position = position_jitter(horizontal_jitter, vertical_jitter))
    }
    else {
        q <- q + geom_point(size = I(cell_size), position = position_jitter(horizontal_jitter,
                                                                            vertical_jitter))
    }
    q <- q + geom_line(aes(x = Pseudotime, y = expectation),
                       data = cds_exprs)
    q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow,
                                          ncol = ncol, scales = "free_y")
    if (min_expr < 1) {
        q <- q + expand_limits(y = c(min_expr, 1))
    }
    if (relative_expr) {
        q <- q + ylab("Relative Expression")
    }
    else {
        q <- q + ylab("Absolute Expression")
    }
    q <- q + xlab("Pseudo-time")
    # q <- q + monocle:::monocle_theme_opts()
    q <- q + scale_colour_manual(values=c("#BBBDC1", "#6B98E9", "#05080D"))# breaks=c(8.5, 11, 15))
    q <- q + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    q
}


plot_genes_in_pseudotime <- function (cds_subset, min_expr = NULL, cell_size = 0.75, nrow = NULL,
    ncol = 1, panel_order = NULL, color_by = "State", trend_formula = "~ sm.ns(Pseudotime, df=3)",
    label_by_short_name = TRUE, relative_expr = TRUE, vertical_jitter = NULL,
    horizontal_jitter = NULL)
{
    f_id <- NA
    Cell <- NA
    if (cds_subset@expressionFamily@vfamily %in% c("negbinomial",
        "negbinomial.size")) {
        integer_expression <- TRUE
    }
    else {
        integer_expression <- FALSE
        relative_expr <- TRUE
    }
    if (integer_expression) {
        cds_exprs <- exprs(cds_subset)
        if (relative_expr) {
            if (is.null(sizeFactors(cds_subset))) {
                stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
            }
            cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
        }
        cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    }
    else {
        cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
    }
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    if (integer_expression) {
        cds_exprs$adjusted_expression <- cds_exprs$expression
    }
    else {
        cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    }
    if (label_by_short_name == TRUE) {
        if (is.null(cds_exprs$gene_short_name) == FALSE) {
            cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
            cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
        }
        else {
            cds_exprs$feature_label <- cds_exprs$f_id
        }
    }
    else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    cds_exprs$f_id <- as.character(cds_exprs$f_id)
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)
    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
    model_expectation <- genSmoothCurves(cds_subset, cores = 1,
        trend_formula = trend_formula, relative_expr = T, new_data = new_data)
    colnames(model_expectation) <- colnames(cds_subset)
    expectation <- plyr::ddply(cds_exprs, plyr::.(f_id, Cell), function(x) data.frame(expectation = model_expectation[x$f_id,
        x$Cell]))
    cds_exprs <- merge(cds_exprs, expectation)
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
    if (is.null(panel_order) == FALSE) {
        cds_exprs$feature_label <- factor(cds_exprs$feature_label,
            levels = panel_order)
    }

    return(cds_exprs[, c('feature_label', 'timepoint', 'expression', "expectation", 'Pseudotime')])


    cds_exprs$timepoint <- factor(cds_exprs$timepoint, levels=sort(unique(cds_exprs$timepoint)))

    q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
    if (is.null(color_by) == FALSE) {
        q <- q + geom_point(aes_string(color = color_by), size = I(cell_size),
            position = position_jitter(horizontal_jitter, vertical_jitter))
    }
    else {
        q <- q + geom_point(size = I(cell_size), position = position_jitter(horizontal_jitter,
            vertical_jitter))
    }
    q <- q + geom_line(aes(x = Pseudotime, y = expectation),
        data = cds_exprs)
    q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow,
        ncol = ncol, scales = "free_y")
    if (min_expr < 1) {
        q <- q + expand_limits(y = c(min_expr, 1))
    }
    if (relative_expr) {
        q <- q + ylab("Relative Expression")
    }
    else {
        q <- q + ylab("Absolute Expression")
    }
    q <- q + xlab("Pseudo-time")
    # q <- q + monocle:::monocle_theme_opts()
    q <- q + scale_colour_manual(values=c("#BBBDC1", "#6B98E9", "#05080D"))# breaks=c(8.5, 11, 15))
    q <- q + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

    q
}



genes_to_TFs <- function(Antler_obj, gene_list, GO_terms = c('GO:0003700', 'GO:0043565')){
    
    # subset ensembl IDs based on genes in gene modules
    gene_ids <- pData(featureData(Antler_obj$expressionSet))[pData(featureData(Antler_obj$expressionSet))$current_gene_names %in% gene_list,]
    
    ensembl = useMart("ensembl",dataset="ggallus_gene_ensembl")
    
    ensembl_list <- getBM(attributes=c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"), 
                          filters = 'ensembl_gene_id', 
                          values = gene_ids$ensembl_gene_id, 
                          mart = ensembl)
    
    # subset genes based on transcription factor GO terms
    # GO:0043565 = sequence specific DNA binding and GO:0003700 = DNA-binding transcription factor activity
    gene_subset <- list()
    for(go in GO_terms){
        gene_subset[[go]] <- ensembl_list$ensembl_gene_id[ensembl_list$go_id %in% go]
    }
    gene_subset <- unique(unlist(gene_subset))
    
    # get current gene names from ensembl IDs
    gene_sub <- pData(featureData(Antler_obj$expressionSet))$current_gene_names[pData(featureData(Antler_obj$expressionSet))$ensembl_gene_id %in% gene_subset]
    
    # filter gene modules by transcription factors
    return(gene_sub)
}


modules_to_TFs <- function(Antler_obj, gene_modules = Antler_obj$topCorr_DR$genemodules.selected, ...){
    gene_sub <- genes_to_TFs(Antler_obj, gene_list = unlist(gene_modules), ...)
    return(lapply(gene_modules, function(x) x[x %in% gene_sub]))
}

