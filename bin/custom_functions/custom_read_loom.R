# read in loom data, with ensembl ID as rownames instead of gene name

custom_read_loom <- function (file, engine = "hdf5r") 
{
  if (engine == "h5") {
    cat("reading loom file via h5...\n")
    f <- h5::h5file(file, mode = "r")
    cells <- f["col_attrs/CellID"][]
    genes <- f["row_attrs/Accession"][]
    dl <- c(spliced = "/layers/spliced", unspliced = "/layers/unspliced", 
            ambiguous = "/layers/ambiguous")
    if ("/layers/spanning" %in% h5::list.datasets(f)) {
      dl <- c(dl, c(spanning = "/layers/spanning"))
    }
    dlist <- lapply(dl, function(path) {
      m <- as(f[path][], "dgCMatrix")
      rownames(m) <- genes
      colnames(m) <- cells
      return(m)
    })
    h5::h5close(f)
    return(dlist)
  }
  else if (engine == "hdf5r") {
    cat("reading loom file via hdf5r...\n")
    f <- hdf5r::H5File$new(file, mode = "r")
    cells <- f[["col_attrs/CellID"]][]
    genes <- f[["row_attrs/Accession"]][]
    dl <- c(spliced = "layers/spliced", unspliced = "layers/unspliced", 
            ambiguous = "layers/ambiguous")
    if ("layers/spanning" %in% hdf5r::list.datasets(f)) {
      dl <- c(dl, c(spanning = "layers/spanning"))
    }
    dlist <- lapply(dl, function(path) {
      m <- as(t(f[[path]][, ]), "dgCMatrix")
      rownames(m) <- genes
      colnames(m) <- cells
      return(m)
    })
    f$close_all()
    return(dlist)
  }
  else {
    warning("Unknown engine. Use hdf5r or h5 to import loom file.")
    return(list())
  }
}
