library(scales)
library(DEP)
library(ComplexHeatmap)
get_annotation_brewer_pal <- function(dep, indicate) {
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.character(indicate))
  
  # Check indicate columns
  col_data <- colData(dep) %>%
    as.data.frame()
  columns <- colnames(col_data)
  if(all(!indicate %in% columns)) {
    stop("'",
         paste0(indicate, collapse = "' and/or '"),
         "' column(s) is/are not present in ",
         deparse(substitute(dep)),
         ".\nValid columns are: '",
         paste(columns, collapse = "', '"),
         "'.",
         call. = FALSE)
  }
  if(any(!indicate %in% columns)) {
    indicate <- indicate[indicate %in% columns]
    warning("Only used the following indicate column(s): '",
            paste0(indicate, collapse = "', '"),
            "'")
  }
  
  # Get annotation
  anno <- select(col_data, indicate)
  
  # Annotation color
  names <- colnames(anno)
  anno_col <- vector(mode="list", length=length(names))
  names(anno_col) <- names
  for(i in names) {
    var = anno[[i]] %>% unique() %>% sort()
    if(length(var) == 1)
      cols <- c("black")
    if(length(var) == 2)
      cols <- c("orangered", "cornflowerblue")
    if(length(var) >2 & length(var) < 8)
      cols <- RColorBrewer::brewer.pal(length(var), "Set2")
    if(length(var) > 8)
      cols <- RColorBrewer::brewer.pal(length(var), "Set3")
    names(cols) <- var
    anno_col[[i]] <-  cols
  }
  
  # HeatmapAnnotation object
  HeatmapAnnotation(df = anno,
                    col = anno_col,
                    show_annotation_name = TRUE)
}
