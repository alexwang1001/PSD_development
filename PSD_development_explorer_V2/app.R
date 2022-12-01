library("shiny")
library("tidyverse")
library("SummarizedExperiment")
options(repos = BiocManager::repositories())

source("helpers.R")
human_dep <- readRDS("data/human_dep.rds")
macaque_dep <- readRDS("data/macaque_dep.rds")
mouse_dep <- readRDS("data/mouse_dep.rds")

human_input <- read.csv("data/human_input.csv")
macaque_input <- read.csv("data/macaque_input.csv")
mouse_input <- read.csv("data/mouse_input.csv")

human_choices <- human_input$Input
macaque_choices <- macaque_input$Input
mouse_choices <- mouse_input$Input

# User interface ----
ui <- fluidPage(
  titlePanel("Postsynaptic density (PSD) development explorer"),
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
              ".shiny-output-error:before { visibility: hidden; }"),
  sidebarLayout(
    sidebarPanel(
      selectizeInput('human_protein', label = "Human protein symbol", choices = NULL, multiple = TRUE),
      selectizeInput('macaque_protein', label = "Macaque protein symbol", choices = NULL, multiple = TRUE),
      selectizeInput('mouse_protein', label = "Mouse protein symbol", choices = NULL, multiple = TRUE),
      # Button
      downloadButton("download_human_selected", "Download_human_selected"),
      br(),
      br(),
      downloadButton("download_macaque_selected", "Download_macaque_selected"),
      br(),
      br(),
      downloadButton("download_mouse_selected", "Download_mouse_selected"),
      br(),
      br(),
      downloadButton("download_all", "Download_whole_datasets"),
      br(),
      br(),
      tags$a(href="https://www.biorxiv.org/content/10.1101/2022.10.24.513541v1", 
             "Link to our publication"),
      width = 5
    ),
    mainPanel(
      textOutput("human_text"),
      plotOutput("human_plot", width = 420, height = 240),
      br(),
      textOutput("macaque_text"),
      plotOutput("macaque_plot", width = 420, height = 240),
      br(),
      textOutput("mouse_text"),
      plotOutput("mouse_plot", width = 420, height = 240),
      width = 7
      )
  )
)


# Server logic ----
server <- function(input, output, session) {
  updateSelectizeInput(session, 'human_protein', choices = human_choices, selected = "DLG4-Homo sapiens", server = TRUE)
  updateSelectizeInput(session, 'macaque_protein', choices = macaque_choices, selected = "DLG4-Macaca mulatta", server = TRUE)
  updateSelectizeInput(session, 'mouse_protein', choices = mouse_choices, selected = "Dlg4-Mus musculus", server = TRUE)
  # Reactive gene names for selected genes ----
  human_protein_table <- reactive({
    human_input[human_input$Input %in% input$human_protein,]
  })
  macaque_protein_table <- reactive({
    macaque_input[macaque_input$Input %in% input$macaque_protein,]
  })
  mouse_protein_table <- reactive({
    mouse_input[mouse_input$Input %in% input$mouse_protein,]
  })
  # Output texts and plots for selected genes ----
  output$human_text <- renderText({
    if (length(human_protein_table()$SYMBOL) != 0) {
      paste0("Human PSD - ", paste(human_protein_table()$SYMBOL, collapse = ", "))
    } else {
      "No human PSD protein selected"
    }
  })
  output$human_plot <- renderPlot({
     if (length(human_protein_table()$SYMBOL) != 0) {
      plot_human(human_dep, human_protein_table()$SYMBOL) + guides(color = guide_legend(override.aes = list(size = 3)))
    } else {
      plot.new()
    }
  })
  output$macaque_text <- renderText({
    if (length(macaque_protein_table()$SYMBOL) != 0) {
      paste0("Macaque PSD - ", paste(macaque_protein_table()$SYMBOL, collapse = ", "))
    } else {
      "No macaque PSD protein selected"
    }
  })
  output$macaque_plot <- renderPlot({
    if (length(macaque_protein_table()$SYMBOL) != 0) {
      plot_macaque(macaque_dep, macaque_protein_table()$SYMBOL) + guides(color = guide_legend(override.aes = list(size = 3)))
    } else {
      plot.new()
    }
  })
  output$mouse_text <- renderText({
    if (length(mouse_protein_table()$SYMBOL) != 0) {
      paste0("Mouse PSD - ", paste(mouse_protein_table()$SYMBOL, collapse = ", "))
    } else {
      "No mouse PSD protein selected"
    }
  })
  output$mouse_plot <- renderPlot({
    if (length(mouse_protein_table()$SYMBOL) != 0) {
      plot_mouse(mouse_dep, mouse_protein_table()$SYMBOL) + guides(color = guide_legend(override.aes = list(size = 3)))
    } else {
      plot.new()
    }
  })
  # Downloadable csv of selected dataset ----
  output$download_human_selected <- downloadHandler(
    filename = function() {
      "Human_dataset_selected.csv"
    },
    content = function(file) {
      if (length(human_protein_table()$SYMBOL) != 0) {
        human_dataset <- get_dataset(human_dep, human_protein_table()$SYMBOL)
        write.csv(human_dataset, file, row.names = FALSE)
      } else {
        NULL
      }
    }
  )
  output$download_macaque_selected <- downloadHandler(
    filename = function() {
      "Macaque_dataset_selected.csv"
    },
    content = function(file) {
      if (length(macaque_protein_table()$SYMBOL) != 0) {
        macaque_dataset <- get_dataset(macaque_dep, macaque_protein_table()$SYMBOL)
        write.csv(macaque_dataset, file, row.names = FALSE)
      } else {
        NULL
      }
    }
  )
  output$download_mouse_selected <- downloadHandler(
    filename = function() {
      "Mouse_dataset_selected.csv"
    },
    content = function(file) {
      if (length(mouse_protein_table()$SYMBOL) != 0) {
        mouse_dataset <- get_dataset(mouse_dep, mouse_protein_table()$SYMBOL)
        write.csv(mouse_dataset, file, row.names = FALSE)
      } else {
        NULL
      }
    }
  )
  
  output$download_all <- downloadHandler(
    filename <- function() {
      "Whole_datasets.zip"
    },
    
    content <- function(file) {
      file.copy("data/Whole_datasets.zip", file)
    },
    contentType = "application/zip"
  )
}


# Run app ----
shinyApp(ui, server)
