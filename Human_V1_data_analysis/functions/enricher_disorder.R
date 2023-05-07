enricher_disorder <- function(gene, universe, TERM2GENE) {
  USER_DATA <- DOSE:::build_Anno(TERM2GENE, NA)
  PATHID2EXTID <- get("PATHID2EXTID", envir = USER_DATA)
  
  PATHID <- names(PATHID2EXTID)
  
  qPATHID2ExtID <- lapply(PATHID2EXTID, intersect, gene)
  uPATHID2EXTID <- lapply(PATHID2EXTID, intersect, universe)
  k <- sapply(qPATHID2ExtID, length)
  M <- sapply(uPATHID2EXTID, length)
  N <- rep(length(universe), length(M))
  n <- rep(length(gene), length(M))
  args.df <- data.frame(numWdrawn=k-1, ## White balls drawn
                        numW=M,        ## White balls
                        numB=N-M,      ## Black balls
                        numDrawn=n)    ## balls drawn
  pvalues <- apply(args.df, 1, function(n)
    phyper(n[1], n[2], n[3], n[4], lower.tail=FALSE)
  )
  GeneRatio <- apply(data.frame(a=k, b=n), 1, function(x)
    paste(x[1], "/", x[2], sep="", collapse="")
  )
  
  BgRatio <- apply(data.frame(a=M, b=N), 1, function(x)
    paste(x[1], "/", x[2], sep="", collapse="")
  )
  GeneRatio_num <- apply(data.frame(a=k, b=n), 1, function(x)
    x[1]/ x[2]
  )
  OddsRatio <- apply(data.frame(a=k, b=n, c=M, d=N), 1, function(x)
    x[1]/x[2]/((x[3]-x[1])/(x[4]-x[2]))
  )
  
  Over <- data.frame(ID = as.character(PATHID),
                     GeneRatio = GeneRatio,
                     BgRatio = BgRatio,
                     pvalue = pvalues,
                     GeneRatio_num = GeneRatio_num,
                     OddsRatio = OddsRatio,
                     stringsAsFactors = FALSE)
  
  p.adj <- p.adjust(Over$pvalue, method="BH")
  geneID <- sapply(qPATHID2ExtID, function(i) paste(i, collapse="/"))
  Over <- data.frame(Over,
                     p.adjust = p.adj,
                     geneID = geneID,
                     Count = k,
                     stringsAsFactors = FALSE)
  
  Over$Description <- Over$ID
  nc <- ncol(Over)
  Over <- Over[, c(1,nc, 2:(nc-1))]
  Over <- Over[order(pvalues),]
  Over$ID <- as.character(Over$ID)
  Over$Description <- as.character(Over$Description)
  row.names(Over) <- as.character(Over$ID)
  return(Over)
}