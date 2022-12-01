
#human plot
plot_human <- function(dep, proteins) {
  row_data <- rowData(dep, use.names = FALSE)
  subset <- dep[proteins]
  df <- data.frame(assay(subset)) %>%
    rownames_to_column() %>%
    gather(ID, val, -rowname) %>%
    left_join(., data.frame(colData(subset)), by = "ID")
  colnames(df)[1] <- "Protein"
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
      theme(axis.text=element_text(size=10, color = "black"),
            axis.text.x = element_text(face = c("plain", "plain", "bold", "plain", "plain", "plain"),hjust = c(0.8,0.2,0.5,0.5,0.5,0.5)),
            axis.title = element_text(size=14) ) +
      geom_vline(xintercept = 8.164907) +
      labs(x = "Post-conceptional age (log-transformed)",
           y = expression("Abundance ("~log[2]~")")) +
      scale_x_continuous(breaks = c(6.9772799235,7.330917,8.164907,9.333155350,10.7648715907,11.6438561898,12.8879821331), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20"))
  } else {
    p <- ggplot(data = df, mapping = aes(x = Log2AgeDays, y = val), color = "black", fill = "black") +
      geom_smooth(color = "red", se = F, show.legend = F) +
      geom_point(size = 0.8, alpha = 0.8, show.legend = F) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(colour = "gray80", linetype = "dashed"),
            panel.border = element_rect(fill = NA),
            panel.background = element_rect(fill = "gray98", color = "black")
      ) +
      theme(axis.text=element_text(size=10, color = "black"),
            axis.text.x = element_text(face = c("plain", "plain", "bold", "plain", "plain", "plain"),hjust = c(0.8,0.2,0.5,0.5,0.5,0.5)),
            axis.title = element_text(size=14) ) +
      geom_vline(xintercept = 8.164907) +
      labs(x = "Post-conceptional age (log-transformed)",
           y = expression("Abundance ("~log[2]~")")) +
      scale_x_continuous(breaks = c(6.9772799235,7.330917,8.164907,9.333155350,10.7648715907,11.6438561898,12.8879821331), labels = c("GW18", "GW23", "Birth", "Year01", "Year04", "Year08", "Year20"))
  }
  return(p)
}

#macaque plot
plot_macaque <- function(dep, proteins) {
  row_data <- rowData(dep, use.names = FALSE)
  subset <- dep[proteins]
  df <- data.frame(assay(subset)) %>%
    rownames_to_column() %>%
    gather(ID, val, -rowname) %>%
    left_join(., data.frame(colData(subset)), by = "ID")
  colnames(df)[1] <- "Protein"
  if (length(proteins) > 1) {
    p <- ggplot(data = df, mapping = aes(x = Log2AgeDays, y = val, color = Protein, fill = Protein)) +
      geom_smooth(se = F, show.legend = F, span = 1.1) +
      geom_point(size = 0.8, alpha = 1, show.legend = T) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(colour = "gray80", linetype = "dashed"),
            panel.border = element_rect(fill = NA),
            panel.background = element_rect(fill = "gray98", color = "black")
      ) +
      theme(axis.text=element_text(size=10, color = "black"),
            axis.text.x = element_text(face = c("plain", "plain", "bold", "plain", "plain", "plain"),hjust = c(0.5,0.5,0.5,0.5,0.5,0.5)),
            axis.title = element_text(size=14) ) +
      geom_vline(xintercept = 7.375039431) +
      labs(x = "Post-conceptional age (log-transformed)",
           y = expression("Abundance ("~log[2]~")")) +
      scale_x_continuous(breaks = c(6.22881869,6.781359714,7.375039431,9.052568051, 10.30035256, 11.89784546), labels = c("E75", "E110", "Birth", "Year01", "Year03", "Year10"))
  } else {
    p <- ggplot(data = df, mapping = aes(x = Log2AgeDays, y = val), color = "black", fill = "black") +
      geom_smooth(color = "darkgreen", se = F, show.legend = F, span = 1.1) +
      geom_point(size = 0.8, alpha = 0.8, show.legend = F) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(colour = "gray80", linetype = "dashed"),
            panel.border = element_rect(fill = NA),
            panel.background = element_rect(fill = "gray98", color = "black")
      ) +
      theme(axis.text=element_text(size=10, color = "black"),
            axis.text.x = element_text(face = c("plain", "plain", "bold", "plain", "plain", "plain"),hjust = c(0.5,0.5,0.5,0.5,0.5,0.5)),
            axis.title = element_text(size=14) ) +
      geom_vline(xintercept = 7.375039431) +
      labs(x = "Post-conceptional age (log-transformed)",
           y = expression("Abundance ("~log[2]~")")) +
      scale_x_continuous(breaks = c(6.22881869,6.781359714,7.375039431,9.052568051, 10.30035256, 11.89784546), labels = c("E75", "E110", "Birth", "Year01", "Year03", "Year10"))
  }
  return(p)
}

#mouse plot
plot_mouse <- function(dep, proteins) {
  row_data <- rowData(dep, use.names = FALSE)
  subset <- dep[proteins]
  df <- data.frame(assay(subset)) %>%
    rownames_to_column() %>%
    gather(ID, val, -rowname) %>%
    left_join(., data.frame(colData(subset)), by = "ID")
  colnames(df)[1] <- "Protein"
  if (length(proteins) > 1) {
    p <- ggplot(data = df, mapping = aes(x = Log2AgeDays, y = val, color = Protein, fill = Protein)) +
      geom_smooth(se = F, show.legend = F, span = 1.1) +
      geom_point(size = 0.8, alpha = 1, show.legend = T) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(colour = "gray80", linetype = "dashed"),
            panel.border = element_rect(fill = NA),
            panel.background = element_rect(fill = "gray98", color = "black")
      ) +
      theme(axis.text=element_text(size=10, color = "black"),
            axis.title = element_text(size=14) ) +
      labs(x = "Post-conceptional age (log-transformed)",
           y = expression("Abundance ("~log[2]~")")) +
      scale_x_continuous(breaks = c(4.24792751344,4.80735492206,5.20945336563,5.78135971352), labels = c("P0", "P9", "P18", "P36"))
  } else {
    p <- ggplot(data = df, mapping = aes(x = Log2AgeDays, y = val), color = "black", fill = "black") +
      geom_smooth(color = "blue", se = F, show.legend = F, span = 1.1) +
      geom_point(size = 0.8, alpha = 0.8, show.legend = F) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(colour = "gray80", linetype = "dashed"),
            panel.border = element_rect(fill = NA),
            panel.background = element_rect(fill = "gray98", color = "black")
      ) +
      theme(axis.text=element_text(size=10, color = "black"),
            axis.title = element_text(size=14) ) +
      labs(x = "Post-conceptional age (log-transformed)",
           y = expression("Abundance ("~log[2]~")")) +
      scale_x_continuous(breaks = c(4.24792751344,4.80735492206,5.20945336563,5.78135971352), labels = c("P0", "P9", "P18", "P36"))
  }
  return(p)
}

#get dataset for downloading
get_dataset <- function(dep, proteins) {
  row_data <- rowData(dep, use.names = FALSE)
  subset <- dep[proteins]
  df <- data.frame(assay(subset)) %>%
    rownames_to_column() %>%
    gather(ID, val, -rowname) %>%
    left_join(., data.frame(colData(subset)), by = "ID")
  colnames(df)[1] <- "Protein"
  df <- df %>%
    select(Protein, ID, condition, Log2AgeDays, val)
  colnames(df) <- c("Protein", "ID", "Condition",
                    "Log2AgeDays", "Abundance")
  df <- df[order(df$Protein, df$ID),]
  return(df)
}