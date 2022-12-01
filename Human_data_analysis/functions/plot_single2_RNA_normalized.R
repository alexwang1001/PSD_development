plot_single2_RNA_normalized <- function(se, genes, scale = FALSE, plot = TRUE, ylim = NULL) {
  # Show error if inputs are not the required classes
  library(SummarizedExperiment)
  library(tidyverse)
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(genes),
                          is.logical(plot),
                          length(plot) == 1)
  
  row_data <- rowData(se, use.names = FALSE)
  
  # Show error if an unvalid gene name is given
  if(all(!genes %in% row_data$SYMBOL)) {
    if(length(genes) == 1) {
      rows <- grep(substr(genes, 1, nchar(genes) - 1),row_data$SYMBOL)
      possibilities <- row_data$SYMBOL[rows]
    } else {
      rows <- lapply(genes, function(x)
        grep(substr(x, 1, nchar(x) - 1),row_data$SYMBOL))
      possibilities <- row_data$SYMBOL[unlist(rows)]
    }
    
    if(length(possibilities) > 0) {
      possibilities_msg <- paste0(
        "Do you mean: '",
        paste0(possibilities, collapse = "', '"),
        "'")
    } else {
      possibilities_msg <- NULL
    }
    stop("please run `plot_single()` with a valid gene names in the 'genes' argument\n",
         possibilities_msg,
         call. = FALSE)
  }
  if(any(!genes %in% row_data$SYMBOL)) {
    genes <- genes[genes %in% row_data$SYMBOL]
    warning("Only used the following gene(s): '",
            paste0(genes, collapse = "', '"),
            "'")
  }
  
  # Single gene
  subset <- se[genes]
  
  # Plot centered log-intensity values
  # Obtain protein-centered fold change values if scale = TRUE
  if (scale) {
    means <- rowMeans(assay(subset), na.rm = TRUE)
    ranges <- rowRanges(assay(subset), na.rm = TRUE)
    ranges_value <- ranges[,2]-ranges[,1]
    df <- data.frame((assay(subset) - means)/ranges_value) %>%
      rownames_to_column() %>%
      gather(sample, val, -rowname) %>%
      left_join(., data.frame(colData(subset)), by = "sample")
    colnames(df)[1] <- "Gene"
    if (length(genes) > 1) {
      p <- ggplot(data = df, mapping = aes(x = log2_age_days, y = val, color = Gene, fill = Gene)) +
        geom_smooth(se = F, show.legend = F) +
        geom_point(size = 0.8, alpha = 0.8, show.legend = T) +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_line(colour = "gray80", linetype = "dashed"),
              panel.border = element_rect(fill = NA),
              panel.background = element_rect(fill = "gray98", color = "black")
        ) +
        theme(axis.text=element_text(size=12, color = "black"),
              axis.text.x = element_text(face = c("plain", "plain", "bold", "plain", "plain", "plain"),hjust = c(0.8,0.2,0.5,0.5,0.5,0.5)),
              axis.title = element_text(size=14) ) +
        geom_vline(xintercept = 8.0552824355) +
        labs(x = "Post-conceptional age (log-transformed)",
             y = "Scaled expression") +
        scale_x_continuous(breaks = c(6.8073549221,7.1996723448,8.0552824355,9.30149619498,10.7532167492,11.6375305515,12.885315061), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20"))
    } else {
      p <- ggplot(data = df, mapping = aes(x = log2_age_days, y = val), color = "black", fill = "black") +
        geom_smooth(color = "blue", se = F, show.legend = F) +
        geom_point(size = 0.8, alpha = 1, show.legend = F) +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_line(colour = "gray80", linetype = "dashed"),
              panel.border = element_rect(fill = NA),
              panel.background = element_rect(fill = "gray98", color = "black")
        ) +
        theme(axis.text=element_text(size=12, color = "black"),
              axis.text.x = element_text(face = c("plain", "plain", "bold", "plain", "plain", "plain"),hjust = c(0.8,0.2,0.5,0.5,0.5,0.5)),
              axis.title = element_text(size=14) ) +
        geom_vline(xintercept = 8.0552824355) +
        labs(x = "Post-conceptional age (log-transformed)",
             y = "Scaled expression") +
        scale_x_continuous(breaks = c(6.8073549221,7.1996723448,8.0552824355,9.30149619498,10.7532167492,11.6375305515,12.885315061), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20"))
    }
    
    if(plot) {
      return(p)
    } else {
      df <- df %>%
        select(Gene, sample, age, log2_age_days, val)
      colnames(df) <- c("Gene", "sample", "age",
                        "log2_age_days", "level")
    }
    return(df)
  } else {
    df <- data.frame(assay(subset)) %>%
      rownames_to_column() %>%
      gather(sample, val, -rowname) %>%
      left_join(., data.frame(colData(subset)), by = "sample")
    colnames(df)[1] <- "Gene"
    if (is.null(ylim)) {
      if (length(genes) > 1) {
        p <- ggplot(data = df, mapping = aes(x = log2_age_days, y = val, color = Gene, fill = Gene)) +
          geom_smooth(se = F, show.legend = F) +
          geom_point(size = 0.8, alpha = 1, show.legend = T) +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.major.x = element_line(colour = "gray80", linetype = "dashed"),
                panel.border = element_rect(fill = NA),
                panel.background = element_rect(fill = "gray98", color = "black")
          ) +
          theme(axis.text=element_text(size=12, color = "black"),
                axis.text.x = element_text(face = c("plain", "plain", "bold", "plain", "plain", "plain"),hjust = c(0.8,0.2,0.5,0.5,0.5,0.5)),
                axis.title = element_text(size=14) ) +
          geom_vline(xintercept = 8.0552824355) +
          labs(x = "Post-conceptional age (log-transformed)",
               y = "Scaled expression") +
          scale_x_continuous(breaks = c(6.8073549221,7.1996723448,8.0552824355,9.30149619498,10.7532167492,11.6375305515,12.885315061), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20"))
      } else {
        p <- ggplot(data = df, mapping = aes(x = log2_age_days, y = val), color = "black", fill = "black") +
          geom_smooth(color = "blue", se = F, show.legend = F) +
          geom_point(size = 0.8, alpha = 0.8, show.legend = F) +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.major.x = element_line(colour = "gray80", linetype = "dashed"),
                panel.border = element_rect(fill = NA),
                panel.background = element_rect(fill = "gray98", color = "black")
          ) +
          theme(axis.text=element_text(size=12, color = "black"),
                axis.text.x = element_text(face = c("plain", "plain", "bold", "plain", "plain", "plain"),hjust = c(0.8,0.2,0.5,0.5,0.5,0.5)),
                axis.title = element_text(size=14) ) +
          geom_vline(xintercept = 8.0552824355) +
          labs(x = "Post-conceptional age (log-transformed)",
               y = "Scaled expression") +
          scale_x_continuous(breaks = c(6.8073549221,7.1996723448,8.0552824355,9.30149619498,10.7532167492,11.6375305515,12.885315061), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20"))
      }
    } else {
      if (length(genes) > 1) {
        p <- ggplot(data = df, mapping = aes(x = log2_age_days, y = val, color = Gene, fill = Gene)) +
          geom_smooth(se = F, show.legend = F) +
          geom_point(size = 0.8, alpha = 1, show.legend = T) +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.major.x = element_line(colour = "gray80", linetype = "dashed"),
                panel.border = element_rect(fill = NA),
                panel.background = element_rect(fill = "gray98", color = "black")
          ) +
          theme(axis.text=element_text(size=12, color = "black"),
                axis.text.x = element_text(face = c("plain", "plain", "bold", "plain", "plain", "plain"),hjust = c(0.8,0.2,0.5,0.5,0.5,0.5)),
                axis.title = element_text(size=14) ) +
          geom_vline(xintercept = 8.0552824355) +
          labs(x = "Post-conceptional age (log-transformed)",
               y = "Scaled expression") +
          scale_x_continuous(breaks = c(6.8073549221,7.1996723448,8.0552824355,9.30149619498,10.7532167492,11.6375305515,12.885315061), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20")) +
          scale_y_continuous(limits = ylim)
      } else {
        p <- ggplot(data = df, mapping = aes(x = log2_age_days, y = val), color = "black", fill = "black") +
          geom_smooth(color = "blue", se = F, show.legend = F) +
          geom_point(size = 0.8, alpha = 0.8, show.legend = F) +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.major.x = element_line(colour = "gray80", linetype = "dashed"),
                panel.border = element_rect(fill = NA),
                panel.background = element_rect(fill = "gray98", color = "black")
          ) +
          theme(axis.text=element_text(size=12, color = "black"),
                axis.text.x = element_text(face = c("plain", "plain", "bold", "plain", "plain", "plain"),hjust = c(0.8,0.2,0.5,0.5,0.5,0.5)),
                axis.title = element_text(size=14) ) +
          geom_vline(xintercept = 8.0552824355) +
          labs(x = "Post-conceptional age (log-transformed)",
               y = "Scaled expression") +
          scale_x_continuous(breaks = c(6.8073549221,7.1996723448,8.0552824355,9.30149619498,10.7532167492,11.6375305515,12.885315061), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20")) +
          scale_y_continuous(limits = ylim)
      }
    }
    
    if(plot) {
      return(p)
    } else {
      df <- df %>%
        select(Gene, sample, age, log2_age_days, val)
      colnames(df) <- c("Gene", "sample", "age",
                        "log2_age_days", "level")
    }
    return(df)
  }
}
