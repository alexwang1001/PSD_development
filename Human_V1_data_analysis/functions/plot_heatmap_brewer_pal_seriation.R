library(DEP)
library(seriation)
plot_heatmap_brewer_pal_seriation <- function(dep, type = c("contrast", "centered"),
                                           kmeans = FALSE, k = 6,
                                           col_limit = 6, indicate = NULL,
                                           clustering_distance = c("euclidean", "maximum", "manhattan", "canberra",
                                                                   "binary", "minkowski", "pearson", "spearman", "kendall", "gower"),
                                           row_font_size = 6, col_font_size = 10, plot = TRUE, ...) {
  
  # Show error if inputs are not the required classes
  if(is.integer(k)) k <- as.numeric(k)
  if(is.integer(col_limit)) col_limit <- as.numeric(col_limit)
  if(is.integer(row_font_size)) row_font_size <- as.numeric(row_font_size)
  if(is.integer(col_font_size)) col_font_size <- as.numeric(col_font_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(type),
                          is.logical(kmeans),
                          is.numeric(k),
                          length(k) == 1,
                          is.numeric(col_limit),
                          length(col_limit) == 1,
                          is.numeric(row_font_size),
                          length(row_font_size) == 1,
                          is.numeric(col_font_size),
                          length(col_font_size) == 1,
                          is.logical(plot),
                          length(plot) == 1)
  
  # Show error if inputs do not contain required columns
  type <- match.arg(type)
  clustering_distance <- match.arg(clustering_distance)
  
  # Extract row and col data
  row_data <- rowData(dep, use.names = FALSE)
  col_data <- colData(dep) %>%
    as.data.frame()
  
  # Show error if inputs do not contain required columns
  if(any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop(paste0("'label', 'condition' and/or 'replicate' columns are not present in '",
                deparse(substitute(dep)), "'"),
         call. = FALSE)
  }
  if(length(grep("_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun test_diff() to obtain the required columns."),
         call. = FALSE)
  }
  if(!"significant" %in% colnames(row_data)) {
    stop(paste0("'significant' column is not present in '",
                deparse(substitute(dep)),
                "'.\nRun add_rejections() to obtain the required column."),
         call. = FALSE)
  }
  
  # Heatmap annotation
  if(!is.null(indicate) & type == "contrast") {
    warning("Heatmap annotation only applicable for type = 'centered'",
            call. = FALSE)
  }
  if(!is.null(indicate) & type == "centered") {
    ha1 <- get_annotation_brewer_pal(dep, indicate)
  } else {
    ha1 <- NULL
  }
  
  # Filter for significant proteins only
  filtered <- dep[row_data$significant, ]
  
  # Check for missing values
  if(any(is.na(assay(filtered)))) {
    warning("Missing values in '", deparse(substitute(dep)), "'. ",
            "Using clustering_distance = 'gower'",
            call. = FALSE)
    clustering_distance <- "gower"
    obs_NA <- TRUE
  } else {
    obs_NA <- FALSE
  }
  
  # Get centered intensity values ('centered')
  if(type == "centered") {
    rowData(filtered)$mean <- rowMeans(assay(filtered), na.rm = TRUE)
    df <- assay(filtered) - rowData(filtered, use.names = FALSE)$mean
  }
  # Get contrast fold changes ('contrast')
  if(type == "contrast") {
    df <- rowData(filtered, use.names = FALSE) %>%
      data.frame() %>%
      column_to_rownames(var = "name") %>%
      select(ends_with("_diff"))
    colnames(df) <-
      gsub("_diff", "", colnames(df)) %>%
      gsub("_vs_", " vs ", .)
    df <- as.matrix(df)
  }
  
  # Facultative kmeans clustering
  if(kmeans & obs_NA) {
    warning("Cannot perform kmeans clustering with missing values",
            call. = FALSE)
    kmeans <- FALSE
  }
  if(kmeans & !obs_NA) {
    set.seed(1)
    df_kmeans <- kmeans(df, k)
    if(type == "centered") {
      # Order the k-means clusters according to the maximum fold change
      # in all samples averaged over the proteins in the cluster
      order <- data.frame(df) %>%
        cbind(., cluster = df_kmeans$cluster) %>%
        mutate(row = apply(.[, seq_len(ncol(.) - 1)], 1, function(x) max(x))) %>%
        group_by(cluster) %>%
        summarize(index = sum(row)/n()) %>%
        arrange(desc(index)) %>%
        pull(cluster) %>%
        match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
    if(type == "contrast") {
      # Order the k-means clusters according to their average fold change
      order <- data.frame(df) %>%
        cbind(df, cluster = df_kmeans$cluster) %>%
        gather(condition, diff, -cluster) %>%
        group_by(cluster) %>%
        summarize(row = mean(diff)) %>%
        arrange(desc(row)) %>%
        pull(cluster) %>%
        match(seq_len(k), .)
      df_kmeans$cluster <- order[df_kmeans$cluster]
    }
  }
  
  if(ncol(df) == 1) {
    col_clust = FALSE
  } else {
    col_clust = TRUE
  }
  if(nrow(df) == 1) {
    row_clust = FALSE
  } else {
    row_clust = TRUE
  }
  if(clustering_distance == "gower") {
    clustering_distance <- function(x) {
      dist <- cluster::daisy(x, metric = "gower")
      dist[is.na(dist)] <- max(dist, na.rm = TRUE)
      return(dist)
    }
  }
  
  # Legend info
  legend <- ifelse(type == "contrast",
                   "log2 Fold change",
                   "log2 Centered intensity")
  
  # Heatmap
  o1 = seriate(dist(df), method = "OLO")
  o2 = seriate(dist(t(df)), method = "OLO")
  ht1 = Heatmap(df,
                col = circlize::colorRamp2(
                  seq(-col_limit, col_limit, (col_limit/5)),
                  rev(RColorBrewer::brewer.pal(11, "RdBu"))),
                split = if(kmeans) {df_kmeans$cluster} else {NULL},
                cluster_rows = as.dendrogram(o1[[1]]),
                cluster_columns = as.dendrogram(o2[[1]]),
                row_names_side = "left",
                column_names_side = "top",
                clustering_distance_rows = clustering_distance,
                clustering_distance_columns = clustering_distance,
                heatmap_legend_param = list(color_bar = "continuous",
                                            legend_direction = "horizontal",
                                            legend_width = unit(5, "cm"),
                                            title_position = "lefttop"),
                name = legend,
                row_names_gp = gpar(fontsize = row_font_size),
                column_names_gp = gpar(fontsize = col_font_size),
                top_annotation = ha1,
                ...)
  if(plot) {
    # Plot
    draw(ht1, heatmap_legend_side = "top")
  } else {
    # Return data.frame
    colnames(df) <- gsub(" ", "_", colnames(df))
    df <- df[, unlist(column_order(ht1))]
    if(kmeans) {
      df <- cbind(df, k = df_kmeans$cluster)
    }
    return <- df[unlist(row_order(ht1)),]
    data.frame(protein = row.names(return), return) %>%
      mutate(order = row_number())
  }
}
