compareCluster_disorder <- function(geneClusters, fun="enricher_disorder", ...) {
  fun_name <- fun
  fun <- eval(parse(text=fun))

  clProf <- plyr::llply(geneClusters,
                  .fun=function(i) {
                    x=suppressMessages(fun(i, ...))
                    as.data.frame(x)
                  }
  )
  clusters.levels = names(geneClusters)
  clProf.df <- plyr::ldply(clProf, rbind)
  
  if (nrow(clProf.df) == 0) {
    stop("No enrichment found in any of gene cluster, please check your input...")
  }
  
  clProf.df <- plyr::rename(clProf.df, c(.id="Cluster"))
  clProf.df$Cluster = factor(clProf.df$Cluster, levels=clusters.levels)
  
  
  ##colnames(clProf.df)[1] <- "Cluster"
  new("compareClusterResult",
      compareClusterResult = clProf.df,
      geneClusters = geneClusters,
      fun = fun_name,
      .call = match.call(expand.dots=TRUE)
  )
}
