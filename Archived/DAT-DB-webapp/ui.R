library(shiny)
library(shinythemes)  
library(DBI)
library(DT)
library(shinyBS)  # For modals and tooltips
library(shinyWidgets)  # For actionBttn
library(data.table)
library(shinyjs)  # Include the shinyjs package
library(RColorBrewer)
library(plotly)
library(duckdb)
library(shinyBS)
library(Biostrings)
library(dplyr)
library(tidyr)


ui <- fluidPage(
  shinyjs::useShinyjs(),  
  theme = shinytheme("sandstone"),  
  
  # Custom CSS for styling headers and overall layout
  tags$head(
    tags$style(HTML("
    /* Set global font */
    body {
      font-family: 'Arial', sans-serif; /* Replace with your preferred font */
    }
    
    /* Customize header font size and style */
    .custom-header {
      font-size: 20px;
      font-weight: bold;
      background-color: #f4f4f9;
      border-left: 6px solid #ff9800;
      padding: 10px;
      margin-bottom: 15px;
      box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
      color: #333;
      border-radius: 5px;
    }
    
    .button-container {
      margin-bottom: 20px;
    }
    
    .equal-height {
      height: 500px;
    }
    
    img {
      width: 100%;
      height: 100%;
      object-fit: contain;
    }
    
    .panel-container {
      border: 1px solid #ddd;
      padding: 20px;
      border-radius: 8px;
      margin-top: 20px;
      background-color: #f9f9f9;
    }

    /* Set the width for Select Data panel */
    .left-panel-top {
      width: 300px;  /* Adjust as necessary to match the width of input fields */
      max-width: 100%;  /* Ensure the width is responsive */
      margin-left: auto;  /* Center the container horizontally */
      margin-right: auto;
    }

  ")),
    tags$title("DAT-DB"),  # Shorter title for the browser tab
  ),
  
  # Main title
  titlePanel(HTML("<strong>D</strong>isease-specific <strong>A</strong>ntibody <strong>T</strong>ryptic peptides <strong>D</strong>atabases for <strong>B</strong>ottom-up proteomics (<strong>DAT-DB</strong>)")),
  navbarPage("",
             
             # Overview Tab
             tabPanel("Overview",
                      fluidPage(
                        tags$div(
                          p(strong("Welcome to DAT-DB app!")),
                          p("The identification of human antibodies using bottom-up proteomics relies on database searches that match experimental peptide fragments to theoretical values derived from protein sequences in databases. 
                            Standard protein databases, such as UniProt and NCBI-RefSeq, contain a limited number of antibody sequences compared to the human body’s capacity to generate billions and currently lack disease-specific antibody sequences. 
                            This limitation could lead to the misidentification of antibodies in samples. 
                            To overcome this challenge, DAT-DB, a database of disease-specific antibodies, provides researchers with antibody tryptic peptides derived from next-generation sequencing of antibody repertoires."),
                          p(strong("Author: Xuan-Tung Trinh (txt@bmb.sdu.dk)")),
                          
                          # Add buttons for toggling images and plot, left-aligned with spacing
                          fluidRow(
                            column(12, 
                                   tags$div(
                                     style = "display: flex; justify-content: flex-start; align-items: flex-start; gap: 20px;",  # Left-align and add space between buttons
                                     
                                     # First button
                                     actionButton("toggle_image1", "What is DAT-DB?", style = "width: 200px; height: 50px; border-radius: 50px; opacity: 0.75"),  # Set button size
                                     
                                     # Second button
                                     actionButton("toggle_image2", "How was DAT-DB created?", style = "width: 200px; height: 50px; border-radius: 50px; opacity: 0.75"),  # Set button size
                                     
                                     # Third button
                                     actionButton("toggle_plot", "What does DAT-DB have?", style = "width: 200px; height: 50px; border-radius: 50px; opacity: 0.75"),  # Set button size
                                     
                                     # Fourth button
                                     actionButton("toggle_image3", "What does DAT-DB give?", style = "width: 200px; height: 50px; border-radius: 50px; opacity: 0.75")  # New button
                                   )
                            )
                          ),
                          
                          # Image 1 with title (below the button)
                          fluidRow(
                            column(12,
                                   tags$div(
                                     id = "image1_container",  # Add an ID for toggling
                                     style = "text-align: center; padding: 15px; height: auto; display: none;",  # Initially hidden
                                     tags$h3("What is DAT-DB?", style = "font-weight: bold;"),  # Title for Image1
                                     img(src = "Image1.svg", style = "max-width: 700px; width: 100%; height: auto;")  # Display Image1
                                   )
                            )
                          ),
                          
                          # Image 2 with title (below the button)
                          fluidRow(
                            column(12,
                                   tags$div(
                                     id = "image2_container",  # Add an ID for toggling
                                     style = "text-align: center; padding: 15px; height: auto; display: none;",  # Initially hidden
                                     
                                     # Title with info button next to it
                                     div(
                                       style = "display: flex; justify-content: center; align-items: center;",  # Align title and info button
                                       
                                       # Title for Image2
                                       tags$h3("How was DAT-DB created?", style = "font-weight: bold; margin-right: 10px;"),  
                                       
                                       # Info button to trigger modal
                                       actionBttn(
                                         inputId = "info_image2_button",
                                         label = NULL, 
                                         icon = icon("info-circle"), 
                                         style = "material-circle", 
                                         color = "primary", 
                                         size = "xs"
                                       )
                                     ),
                                     
                                     # Display Image2
                                     img(src = "Image2.svg", style = "max-width: 600px; width: 100%; height: auto;")
                                   )
                            )
                          ),
                          
                          # Create the modal that will pop up when the info button is clicked
                          bsModal(
                            id = "infoModal", 
                            title = "More Information", 
                            trigger = "info_image2_button",  # The button that triggers the modal
                            size = "medium",
                            
                            # Modal content
                            p("More information about the workflow is available on our GitHub page:"),
                            tags$a(href = "https://github.com/trinhxt/SDU_Immunoinformatics", 
                                   target = "_blank",  # Opens the link in a new tab
                                   "Visit our GitHub for more details!")
                          ),
                          
                          # Scatter Plot with title (below the button), center-aligned
                          fluidRow(
                            column(12,
                                   tags$div(
                                     id = "plot_container",  # Add an ID for toggling
                                     style = "padding: 15px; height: auto; display: none;",  # Initially hidden
                                     
                                     # Title centered
                                     tags$h3("What does DAT-DB have?", style = "text-align: center; font-weight: bold;"),  # Title for Scatter Plot
                                     
                                     # Scatter plot container with Flexbox centering
                                     tags$div(
                                       style = "display: flex; justify-content: center;",  # Flexbox to center the plotly plot
                                       plotlyOutput("scatterPlot", height = "600px", width = "900px")  # Set size for scatter plot
                                     )
                                   )
                            )
                          ),
                          
                          # Image 3 with title (below the new button)
                          fluidRow(
                            column(12,
                                   tags$div(
                                     id = "image3_container",  # Add an ID for toggling
                                     style = "text-align: center; padding: 15px; height: auto; display: none;",  # Initially hidden
                                     tags$h3("What does DAT-DB give?", style = "font-weight: bold;"),  # Title for Image3
                                     img(src = "Image3.svg", style = "max-width: 900px; width: 100%; height: auto;")  # Display Image3
                                   )
                            )
                          )
                        )
                      )
             ),
             
             # Download Peptides
             tabPanel("Download Peptides",
                      fluidPage(
                        
                        # Create the grid layout
                        fluidRow(
                          
                          # Left-side: 25% width for Filter panel
                          column(3,
                                 div(class = "panel-container left-panel-top",
                                     # Filter Area
                                     div(
                                       tags$div(class = "custom-header", "Select Data"),
                                       
                                       # Dropdown to select Disease (with None option)
                                       selectInput("disease", "Select Disease:", choices = c("None")),
                                       
                                       # Dropdown to select BSource (with All option)
                                       selectInput("bsource", "B-cells source:", choices = c("All")),
                                       
                                       # Dropdown to select BType (with All option)
                                       selectInput("btype", "B-cells type:", choices = c("All")),
                                       
                                       # Dropdown to select Isotype (with All option)
                                       selectInput("isotype", "Antibody isotype:", choices = c("All")),
                                       
                                       # Unique peptides (Number of Antibodies Range)
                                       div(
                                         class = "input-with-button",
                                         tags$label("Unique peptides?"),
                                         actionBttn("info_antibody", label = NULL, icon = icon("info-circle"), style = "pill", color = "primary", size = "xs")
                                       ),
                                       selectInput("antibody_range", label = NULL, choices = list("Yes", "No", "Either"), selected = "Yes"),
                                       
                                       # Number of antibodies slider (conditionally displayed)
                                       div(
                                         # Hide the specific "Number of antibodies" section by default
                                         shinyjs::hidden(
                                           div(
                                             id = "n_antibody_range",  # Add a unique ID to the section
                                             class = "input-with-button",
                                             tags$label("Number of antibodies:"),
                                             actionBttn("info_N_antibody", label = NULL, icon = icon("info-circle"), style = "pill", color = "primary", size = "xs"),
                                             sliderInput("n_antibody_range", label = NULL, min = 1, max = 10, value = c(1, 10))
                                           )
                                         )
                                       ),
                                       
                                       # Present in CDR3 region
                                       div(
                                         class = "input-with-button",
                                         tags$label("Peptides in CDR3 region?"),
                                         actionBttn("info_cdr3", label = NULL, icon = icon("info-circle"), style = "pill", color = "primary", size = "xs")
                                       ),
                                       selectInput("cdr3_value", label = NULL, choices = list("Yes", "No", "Either"), selected = "Yes"),
                                       
                                       # Number of patients slider
                                       div(
                                         class = "input-with-button",
                                         tags$label("Number of patients:"),
                                         actionBttn("info_n_patient", label = NULL, icon = icon("info-circle"), style = "pill", color = "primary", size = "xs")
                                       ),
                                       sliderInput("n_patient_range", label = NULL, min = 1, max = 10, value = c(1, 10))
                                     )
                                 )
                          ),
                          
                          # Right-side: 75% width for filtered data and download buttons
                          column(9,
                                 div(class = "panel-container right-panel",
                                     tags$div(class = "custom-header", "Download Data"),
                                     fluidRow(
                                       # piechart1 on the left
                                       column(6,  
                                              div(class = "panel-container",
                                                  plotlyOutput("piechart1", width = "100%")  # Add piechart1 here
                                              )
                                       ),
                                       # piechart2 on the right
                                       column(6,
                                              div(class = "panel-container",
                                                  shinyjs::hidden(plotlyOutput("piechart2"))  # piechart2 remains hidden initially
                                              )
                                       )
                                     ),
                                     shinyjs::hidden(uiOutput("filteredSummary")),  # Filtered summary remains outside the pie chart row
                                     
                                     # Download buttons next to summary
                                     fluidRow(
                                       div(class = "button-container",
                                           column(2, shinyjs::hidden(downloadButton("downloadDataCSV", "Download CSV"))),
                                           column(2, shinyjs::hidden(downloadButton("downloadDataFASTA", "Download FASTA"))),
                                           
                                           # Button for downloadDataFASTAplus and info button next to it
                                           column(3,
                                                  div(
                                                    class = "input-with-button",
                                                    shinyjs::hidden(downloadButton("downloadDataFASTAplus", "Download FASTA+")),
                                                    actionBttn("info_fasta_plus", label = NULL, icon = icon("info-circle"), style = "pill", color = "primary", size = "xs")
                                                  )
                                           )
                                       )
                                     )
                                 )
                          )
                          
                        )
                      )
             )
  )
)
