plot_single2 <- function(dep, proteins, scale = FALSE, plot = TRUE, ylim = NULL) {
  # Show error if inputs are not the required classes
  library(DEP)
  library("tidyverse")
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(proteins),
                          is.logical(plot),
                          length(plot) == 1)
  
  row_data <- rowData(dep, use.names = FALSE)
  
  # Show error if an unvalid protein name is given
  if(all(!proteins %in% row_data$name)) {
    if(length(proteins) == 1) {
      rows <- grep(substr(proteins, 1, nchar(proteins) - 1),row_data$name)
      possibilities <- row_data$name[rows]
    } else {
      rows <- lapply(proteins, function(x)
        grep(substr(x, 1, nchar(x) - 1),row_data$name))
      possibilities <- row_data$name[unlist(rows)]
    }
    
    if(length(possibilities) > 0) {
      possibilities_msg <- paste0(
        "Do you mean: '",
        paste0(possibilities, collapse = "', '"),
        "'")
    } else {
      possibilities_msg <- NULL
    }
    stop("please run `plot_single()` with a valid protein names in the 'proteins' argument\n",
         possibilities_msg,
         call. = FALSE)
  }
  if(any(!proteins %in% row_data$name)) {
    proteins <- proteins[proteins %in% row_data$name]
    warning("Only used the following protein(s): '",
            paste0(proteins, collapse = "', '"),
            "'")
  }
  
  # Single protein
  subset <- dep[proteins]
  
  # Plot centered log-intensity values
  # Obtain protein-centered fold change values if scale = TRUE
  if (scale) {
    means <- rowMeans(assay(subset), na.rm = TRUE)
    df <- data.frame(assay(subset) - means) %>%
      rownames_to_column() %>%
      gather(ID, val, -rowname) %>%
      left_join(., data.frame(colData(subset)), by = "ID")
    colnames(df)[1] <- "Protein"
    if (length(proteins) > 1) {
      p <- ggplot(data = df, mapping = aes(x = Log2AgeDays, y = val, color = Protein, fill = Protein)) +
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
             y = "Scaled abundance") +
        scale_x_continuous(breaks = c(6.8073549221,7.1996723448,8.0552824355,9.30149619498,10.7532167492,11.6375305515,12.885315061), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20"))
    } else {
      p <- ggplot(data = df, mapping = aes(x = Log2AgeDays, y = val), color = "black", fill = "black") +
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
             y = "Scaled abundance") +
        scale_x_continuous(breaks = c(6.8073549221,7.1996723448,8.0552824355,9.30149619498,10.7532167492,11.6375305515,12.885315061), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20"))
    }
    
    if(plot) {
      return(p)
    } else {
      df <- df %>%
        select(Protein, ID, condition, Log2AgeDays, val)
      colnames(df) <- c("protein", "ID", "condition",
                        "Log2AgeDays", "level")
    }
    return(df)
  } else {
    df <- data.frame(assay(subset)) %>%
      rownames_to_column() %>%
      gather(ID, val, -rowname) %>%
      left_join(., data.frame(colData(subset)), by = "ID")
    colnames(df)[1] <- "Protein"
    if (is.null(ylim)) {
      if (length(proteins) > 1) {
        p <- ggplot(data = df, mapping = aes(x = Log2AgeDays, y = val, color = Protein, fill = Protein)) +
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
               y = expression("Abundance ("~log[2]~")")) +
          scale_x_continuous(breaks = c(6.8073549221,7.1996723448,8.0552824355,9.30149619498,10.7532167492,11.6375305515,12.885315061), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20"))
      } else {
        p <- ggplot(data = df, mapping = aes(x = Log2AgeDays, y = val), color = "black", fill = "black") +
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
               y = expression("Abundance ("~log[2]~")")) +
          scale_x_continuous(breaks = c(6.8073549221,7.1996723448,8.0552824355,9.30149619498,10.7532167492,11.6375305515,12.885315061), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20"))
      }
    } else {
      if (length(proteins) > 1) {
        p <- ggplot(data = df, mapping = aes(x = Log2AgeDays, y = val, color = Protein, fill = Protein)) +
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
               y = expression("Abundance ("~log[2]~")")) +
          scale_x_continuous(breaks = c(6.8073549221,7.1996723448,8.0552824355,9.30149619498,10.7532167492,11.6375305515,12.885315061), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20")) +
          scale_y_continuous(limits = ylim)
      } else {
        p <- ggplot(data = df, mapping = aes(x = Log2AgeDays, y = val), color = "black", fill = "black") +
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
               y = expression("Abundance ("~log[2]~")")) +
          scale_x_continuous(breaks = c(6.8073549221,7.1996723448,8.0552824355,9.30149619498,10.7532167492,11.6375305515,12.885315061), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20")) +
          scale_y_continuous(limits = ylim)
      }
    }
    
    if(plot) {
      return(p)
    } else {
      df <- df %>%
        select(Protein, ID, condition, Log2AgeDays, val)
      colnames(df) <- c("protein", "ID", "condition",
                        "Log2AgeDays", "level")
    }
    return(df)
  }
}

