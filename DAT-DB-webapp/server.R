server <- function(input, output, session) {
  
  # Read the metadata2.csv file from the www folder
  db_list <- read.csv(file.path("www", "metadata2.csv"), stringsAsFactors = FALSE)
  
  # List the database files in the "pepDB" folder (including extension)
  db_files <- list.files(path = "pepDB", pattern = "\\.duckdb$", full.names = FALSE)
  
  # Filter the db_list to include only rows where File_db matches the files in pepDB
  filtered_db_list <- db_list[db_list$File_db %in% db_files, ]
  
  # Extract unique diseases based on the filtered db_list
  disease_choices <- c("None", sort(unique(filtered_db_list$Disease)))
  
  # Initially hide the images and plot
  shinyjs::show("image1_container")
  shinyjs::hide("image2_container")
  shinyjs::hide("image3_container")
  shinyjs::hide("plot_container")
  
  # When Image 1 button is clicked, show Image 1 and hide the others
  observeEvent(input$toggle_image1, {
    shinyjs::show("image1_container")
    shinyjs::hide("image2_container")
    shinyjs::hide("image3_container")
    shinyjs::hide("plot_container")
  })
  
  # When Image 2 button is clicked, show Image 2 and hide the others
  observeEvent(input$toggle_image2, {
    shinyjs::hide("image1_container")
    shinyjs::show("image2_container")
    shinyjs::hide("image3_container")
    shinyjs::hide("plot_container")
  })
  
  # When Image 3 button is clicked, show Image 3 and hide the others
  observeEvent(input$toggle_image3, {
    shinyjs::hide("image1_container")
    shinyjs::hide("image2_container")
    shinyjs::show("image3_container")
    shinyjs::hide("plot_container")
  })
  
  # When Plot button is clicked, show Plot and hide the others
  observeEvent(input$toggle_plot, {
    shinyjs::hide("image1_container")
    shinyjs::hide("image2_container")
    shinyjs::hide("image3_container")
    shinyjs::show("plot_container")
  })
  
  # Reactive output for the scatter plot
  output$scatterPlot <- renderPlotly({
    load("www/scatterPlot.RData")
    scatterPlot
  })
  
  # Non-reactive variable to store the database connection (avoid reactive context issues)
  con <- NULL
  
  # Ensure that the connection is closed when the session ends
  session$onSessionEnded(function() {
    if (!is.null(con)) {
      dbDisconnect(con)
    }
  })
  
  
  # Populate the disease selectInput based on the sorted metadata2.csv
  observe({
    updateSelectInput(session, "disease", choices = disease_choices, selected = "None")
    
    # Initially hide download buttons and outputs
    shinyjs::hide("downloadDataCSV")
    shinyjs::hide("downloadDataFASTA")
    shinyjs::hide("downloadDataFASTAdecoy")
    shinyjs::hide("filteredSummary")
    shinyjs::hide("piechart2")
    shinyjs::hide("piechart1")
    shinyjs::hide("info_fasta_decoy")
  })
  
  
  # Load the selected database and store the connection
  observeEvent(input$disease, {
    if (input$disease == "None") {
      shinyjs::hide("downloadDataCSV")
      shinyjs::hide("downloadDataFASTA")
      shinyjs::hide("downloadDataFASTAdecoy")
      shinyjs::hide("filteredSummary")
      shinyjs::hide("piechart2")
      shinyjs::hide("piechart1")
      shinyjs::hide("info_fasta_decoy")
      
      # Clear the database connection and reset filteredData
      if (!is.null(con)) {
        dbDisconnect(con)
        con <<- NULL
      }
      
      # Clear any reactive data (filteredData, plot, etc.)
      filteredData <- reactiveVal(NULL)  # Reset the reactive data
      gc()
      
      return(NULL)
    }
    
    req(input$disease)  
    
    gc()
    
    
    # Wrap everything in a progress bar
    withProgress(message = 'Processing...', value = 0, {
      incProgress(0.2, detail = "Loading disease data...")
      
      # Open a new database connection
      selected_file <- db_list$File_db[db_list$Disease == input$disease]
      
      if (length(selected_file) != 1 || selected_file == "") {
        sendSweetAlert(session, title = "Error", text = "Invalid database file for the selected disease.", type = "error")
        return(NULL)
      }
      
      # Close previous connection if one exists
      if (!is.null(con)) {
        dbDisconnect(con)
        con <<- NULL  # Reset the connection variable
        gc()  
      }
      
      # Connect to the selected database
      con <<- tryCatch({
        dbConnect(duckdb::duckdb(), dbdir = file.path("pepDB", selected_file))
      }, error = function(e) {
        sendSweetAlert(session, title = "Error", text = "Error opening database file.", type = "error")
        NULL
      })
      
      
      
      incProgress(0.4, detail = "Fetching data for filters...") # Progess bar
      # Fetch unique values for BSource, BType, and Isotype from the selected database
      bsource_values <- dbGetQuery(con, "SELECT DISTINCT BSource FROM DATDB")
      btype_values <- dbGetQuery(con, "SELECT DISTINCT BType FROM DATDB")
      isotype_values <- dbGetQuery(con, "SELECT DISTINCT Isotype FROM DATDB")
      
      # Add "All" option to BSource, BType, and Isotype dropdowns
      bsource_choices <- c("All", bsource_values$BSource)
      btype_choices <- c("All", btype_values$BType)
      isotype_choices <- c("All", isotype_values$Isotype)
      
      # Update BSource, BType, and Isotype dropdowns dynamically
      updateSelectInput(session, "bsource", choices = bsource_choices)
      updateSelectInput(session, "btype", choices = btype_choices)
      updateSelectInput(session, "isotype", choices = isotype_choices)
      
      # Immediately update and show piechart1 and 2 and filteredSummary after disease selection
      shinyjs::show("filteredSummary")
      shinyjs::show("piechart2")
      shinyjs::show("piechart1")
      shinyjs::show("downloadDataCSV") # Show download buttons when a disease is selected
      shinyjs::show("downloadDataFASTA")  # Show download buttons when a disease is selected
      shinyjs::show("downloadDataFASTAdecoy")
      shinyjs::show("info_fasta_decoy")
      gc()
      
      
      
      incProgress(0.6, detail = "Filtering data...") # Progress bar
      # Reactive expression to filter data based on user input and selected database
      filteredData <- reactive({
        req(con)  
        
        # Nested progress bar for filtering data
        withProgress(message = 'Filtering data...', value = 0, {
          
          # Increment the progress bar during each step
          incProgress(0.2, detail = "Building query...")
          
          # Determine the antibody condition based on user selection
          if (input$antibody_range == "Either") {
            antibody_condition <- "N_antibody >= 1"
          } else if (input$antibody_range == "Yes") {
            antibody_condition <- "N_antibody = 1"
          } else {
            antibody_condition <- "N_antibody > 1"
          }
          
          # Determine the CDR3 condition based on user selection
          if (input$cdr3_value == "Either") {
            cdr3_condition <- "CDR3 IN (0,1)"
          } else if (input$cdr3_value == "Yes") {
            cdr3_condition <- "CDR3 = 1"
          } else {
            cdr3_condition <- "CDR3 = 0"
          }
          
          # Add N_patient condition based on the slider range
          n_patient_condition <- sprintf("N_patient BETWEEN %d AND %d", input$n_patient_range[1], input$n_patient_range[2])
          # Add N_antibody condition based on the slider range
          n_antibody_condition <- sprintf("N_antibody BETWEEN %d AND %d", input$n_antibody_range[1], input$n_antibody_range[2])
          
          # Add BSource condition (if All is selected, no filter applied)
          bsource_condition <- if (input$bsource == "All") "1=1" else sprintf("BSource = '%s'", input$bsource)
          
          # Add BType condition (if All is selected, no filter applied)
          btype_condition <- if (input$btype == "All") "1=1" else sprintf("BType = '%s'", input$btype)
          
          # Add Isotype condition (if All is selected, no filter applied)
          isotype_condition <- if (input$isotype == "All") "1=1" else sprintf("Isotype = '%s'", input$isotype)
          
          
          # Increment progress
          incProgress(0.4, detail = "Executing query...")
          # SQL query to filter the data based on all conditions
          query <- sprintf("SELECT Sequence, BSource, BType, Isotype, N_patient, N_antibody, Length_aa, CDR3 
                          FROM DATDB WHERE %s AND %s AND %s AND %s AND %s AND %s AND %s",
                           antibody_condition,
                           cdr3_condition,
                           n_patient_condition,
                           n_antibody_condition,
                           bsource_condition,
                           btype_condition,
                           isotype_condition)
          
          # Execute query
          data <- tryCatch({
            dbGetQuery(con, query)
          }, error = function(e) {
            sendSweetAlert(session, title = "Error", text = "Error filtering data from the database.", type = "error")
            NULL
          })
          
          
          # Increment progress
          incProgress(0.8, detail = "Processing data...")
          # Return the filtered data
          data
          
        })
      })
      
      
      incProgress(0.8, detail = "Rendering plots and summary...") # Progress bar
      
      # Pie chart 1 (plotly) output for BSource and BType
      output$piechart1 <- renderPlotly({
        req(filteredData())  
        
        data <- filteredData()  
        
        if (is.null(data) || nrow(data) == 0) {
          return(NULL)
        }
        
        # Read metadata (assuming this loads your metadata correctly)
        metadata <- fread(file.path("www",   "metadata1.csv"))
        
        # Define the BSource and BType levels
        BSource <- unique(metadata$BSource)
        BType <- unique(metadata$BType)
        
        # Group and count sequences for each hierarchical level (BSource -> BType)
        data2 <- data %>%
          # Calculate the count of sequences for each BSource
          group_by(BSource) %>%
          summarise(N_BSource = n(), .groups = 'drop') %>%
          
          # Left join the count of sequences for each BType within each BSource
          left_join(
            data %>%
              group_by(BSource, BType) %>%
              summarise(N_BType = n(), .groups = 'drop'),
            by = "BSource"
          )
        
        # Create labels and parents for the sunburst plot
        labels <- c(
          input$disease,  # Root node
          unique(data2$BSource),  # Level 1: BSource
          paste(data2$BSource, data2$BType, sep = ":")  # Level 2: BType
        )
        
        parents <- c(
          "",  # Root node has no parent
          rep(input$disease, length(unique(data2$BSource))),  # Level 1: Parent is Disease
          data2$BSource  # Level 2: Parent is BSource
        )
        
        # Combine the values for each level
        values <- c(
          sum(data2$N_BType),  # Root: sum of all sequences
          as.numeric(unique(data2[,1:2])[[2]]),
          data2$N_BType  # Level 2: sum by BType
        )
        # Ensure values are numeric
        values <- as.numeric(values)
        
        # Handle dynamic color palette ensuring a minimum of 3 colors
        n_colors <- length(labels) - 1  # Subtracting 1 to account for the root node "Disease"
        n_colors <- max(n_colors, 3)  # Ensure at least 3 colors
        
        # Create color palette using Pastel2 (maximum 8 colors) and dynamically extend if needed
        if (n_colors <= 8) {
          color_palette <- brewer.pal(n_colors, "Pastel2")
        } else {
          color_palette <- colorRampPalette(brewer.pal(8, "Pastel2"))(n_colors)
        }
        
        # Create the sunburst plot with BSource and BType levels and apply the color palette
        fig <- plot_ly(
          labels = labels,
          parents = parents,
          values = values,  # The size of each node
          type = 'sunburst',
          marker = list(colors = c("rgba(0, 0, 0, 0)", color_palette)),  # Adding white for the root node and the generated color palette
          branchvalues = 'total',
          textinfo = 'label+percent entry',  # Show label and percentage
          insidetextorientation = 'tangential'  # Ensure labels are shown inside the plot in a readable way
        ) %>%
          layout(
            title = list(
              text = "<b>Percentage of peptides by <br> source and type of B-cells<b><br>",
              x = 0.5,  # Center the title horizontally
              y = 0.95,  # Adjust vertical position
              font = list(family = "Arial", size = 18),
              yanchor = "bottom"
            ),
            margin = list(t = 50, l = 0, r = 0, b = 0),  # Adjust margins for proper spacing
            font = list(family = "Arial", size = 14),  # Global font settings
            plot_bgcolor = "rgba(0, 0, 0, 0)",  # Transparent plot background
            paper_bgcolor = "rgba(0, 0, 0, 0)",  # Transparent page background
            showlegend = FALSE,  # Hide the legend for a cleaner look
            autosize = TRUE  # Let the plot adjust automatically
          )
        
        # Display the plot
        fig
        
      })
      
      
      
      
      # Pie chart 2 (plotly) output for Isotype and uniqueness
      output$piechart2 <- renderPlotly({
        req(filteredData())  
        
        data <- filteredData()  
        
        if (is.null(data) || nrow(data) == 0) {
          return(NULL)
        }
        
        # Prepare data for plotly 
        isotypes <- c("IGHA", "IGHD", "IGHE", "IGHG", "IGHM", "Bulk")
        uniqueness <- c("Unique", "non-Unique")
        
        isotype_counts <- data %>%
          filter(Isotype %in% isotypes) %>%
          group_by(Isotype) %>%
          summarise(N_peptide_Isotype = n())  
        
        uniqueness_counts <- data %>%
          filter(Isotype %in% isotypes) %>%
          mutate(Uniqueness = ifelse(N_antibody == 1, "Unique", "non-Unique")) %>%
          group_by(Isotype, Uniqueness) %>%
          summarise(N_peptide_Uniqueness = n())  
        
        complete_grid <- expand.grid(Isotype = isotypes, Uniqueness = uniqueness)
        
        uniqueness_counts_complete <- complete_grid %>%
          left_join(uniqueness_counts, by = c("Isotype", "Uniqueness")) %>%
          mutate(N_peptide_Uniqueness = replace_na(N_peptide_Uniqueness, 0))  
        
        data2 <- uniqueness_counts_complete %>%
          left_join(isotype_counts, by = "Isotype") %>%
          mutate(N_peptide_Isotype = replace_na(N_peptide_Isotype, 0))  
        
        labels <- c(input$disease,                    
                    data2$Isotype[1:6],              
                    paste(data2$Isotype, data2$Uniqueness, sep = ":"))  
        
        parents <- c("",                                
                     rep(input$disease, 6),              
                     data2$Isotype)                     
        
        values <- c(sum(data2$N_peptide_Uniqueness),    
                    data2$N_peptide_Isotype[1:6],       
                    data2$N_peptide_Uniqueness)         
        
        isotype_colors <- brewer.pal(6, "Pastel1")  
        
        uniqueness_colors <- sapply(isotype_colors, function(col) adjustcolor(col, alpha.f = 0.6))
        
        custom_colors <- c("rgba(0, 0, 0, 0)", isotype_colors, uniqueness_colors)  
        
        fig <- plot_ly(
          labels = labels,     
          parents = parents,   
          values = values,     
          type = 'sunburst',   
          branchvalues = 'total',  
          marker = list(colors = custom_colors)  ,
          textinfo = 'label+percent entry',  # Show label and percentage
          insidetextorientation = 'tangential'  # Ensure labels are shown inside the plot in a readable way
        ) %>%
          layout(
            title = list(
              text = "<b>Percentage of peptides by <br> antibody isotypes and peptide uniqueness <b> <br>",  # Set your title here
              x = 0.5,  # Center the title horizontally
              y = 0.95,  # Adjust vertical position
              font = list(family = "Arial", size = 18),
              yanchor = "bottom"
            ),
            margin = list(t = 50, l = 0, r = 0, b = 0),  # Adjust margins for proper spacing
            font = list(family = "Arial", size = 14),  # Global font settings
            plot_bgcolor = "rgba(0, 0, 0, 0)",  # Transparent plot background
            paper_bgcolor = "rgba(0, 0, 0, 0)",  # Transparent page background
            showlegend = FALSE,  # Hide the legend for a cleaner look
            autosize = TRUE  # Let the plot adjust automatically
          )
        
        # Show the plot
        fig
        
      })
      
      # Output: Render the filtered data summary
      output$filteredSummary <- renderUI({
        req(filteredData())  
        
        data <- filteredData()  
        
        if (is.null(data) || nrow(data) == 0) {
          return(h4("Filtered data is empty"))
        }
        
        # Summary calculations
        num_sequences <- nrow(data)
        unique_peptides <- sum(data$N_antibody == 1, na.rm = TRUE)  
        non_unique_peptides <- sum(data$N_antibody > 1, na.rm = TRUE)  
        cdr3_count_1 <- sum(data$CDR3 == 1, na.rm = TRUE)
        cdr3_count_0 <- sum(data$CDR3 == 0, na.rm = TRUE)
        
        # Display summary
        return(tagList(
          h4(HTML("<strong style='font-size:16px;'>Summary:</strong>")),
          p(HTML(paste("<span style='font-size:16px;'>Number of peptides:</span> <strong style='font-size:16px;'>", formatC(num_sequences, format = "d", big.mark = ","), "</strong>"))),
          p(HTML(paste("<span style='font-size:16px;'>Unique peptides:</span> <strong style='font-size:16px;'>", formatC(unique_peptides, format = "d", big.mark = ","), "</strong>"))),
          p(HTML(paste("<span style='font-size:16px;'>Non-unique peptides:</span> <strong style='font-size:16px;'>", formatC(non_unique_peptides, format = "d", big.mark = ","), "</strong>"))),
          p(HTML(paste("<span style='font-size:16px;'>In CDR3 region:</span> <strong style='font-size:16px;'>", formatC(cdr3_count_1, format = "d", big.mark = ","), "</strong>"))),
          p(HTML(paste("<span style='font-size:16px;'>Not in CDR3 region:</span> <strong style='font-size:16px;'>", formatC(cdr3_count_0, format = "d", big.mark = ","), "</strong>")))
        ))
        
      })
      
      # clean memory
      gc()
      # Complete progress bar
      incProgress(1, detail = "Complete")
    })
    
    
    
  })
  
  
  # Dynamically show/hide the "Number of antibodies" slider based on "Unique peptides?" selection
  observeEvent(input$antibody_range, {
    if (input$antibody_range == "No" || input$antibody_range == "Either") {
      shinyjs::show("n_antibody_range")
    } else {
      shinyjs::hide("n_antibody_range")
    }
  })
  
  
  # Dynamically set the max of the N_patient slider based on the max value in the data
  observe({
    if (input$disease == "None") {
      max_n_patient <- 10
    } else  { max_n_patient <- max(db_list$max_patient[db_list$Disease == input$disease], na.rm = TRUE) }
    updateSliderInput(session, "n_patient_range", 
                      min = 1, 
                      max = max_n_patient, 
                      value = c(1, max_n_patient))  # Set max dynamically
  })
  
  
  # Dynamically set the max of the N_antibody slider based on the max value in the data
  observe({
    if (input$disease == "None") {
      max_n_antibody <- 10
    } else  { max_n_antibody <- max(db_list$max_antibody[db_list$Disease == input$disease], na.rm = TRUE) }
    updateSliderInput(session, "n_antibody_range", 
                      min = 1, 
                      max = max_n_antibody, 
                      value = c(1, max_n_antibody))  # Set max dynamically
  })
  
  
  
  # Download handler for downloading the filtered data as CSV
  output$downloadDataCSV <- downloadHandler(
    filename = function() {
      paste(input$disease, "_filtered_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      data <- filteredData()
      fwrite(data, file)
    }
  )
  
  # Download handler for downloading the filtered data as FASTA using `Biostrings::writeXStringSet`
  output$downloadDataFASTA <- downloadHandler(
    filename = function() {
      paste(input$disease, "_", Sys.Date(), ".fasta", sep = "")
    },
    content = function(file) {
      # Start progress bar
      withProgress(message = "Preparing FASTA file...", value = 0, {
        
        # Simulate progress (this is optional; adjust depending on your app's speed)
        incProgress(0.2, detail = "Filtering data...")
        
        # Get the filtered data
        data <- filteredData()
        
        # Simulate progress
        incProgress(0.4, detail = "Creating sequences...")
        
        # Create headers in OAS/UniProt-like format
        headers <- paste("OAS", input$disease, data$BSource, data$BType, data$Isotype, 
                         paste0("Nab_", data$N_antibody), sep = "|")
        
        # Create AAStringSet object (for protein sequences)
        fasta_data <- AAStringSet(data$Sequence)
        
        # Assign headers to the sequence names
        names(fasta_data) <- headers
        
        # Simulate progress
        incProgress(0.6, detail = "Writing FASTA file...")
        
        # Write the FASTA file using Biostrings
        writeXStringSet(fasta_data, filepath = file, format = "fasta")
        
        # Clean up memory
        rm(data, headers,fasta_data)
        gc()  # Trigger garbage collection to free up memory
        
        # Complete progress
        incProgress(1, detail = "Done.")
      })
    }
  )
  
  
  # Download handler for downloading the filtered data as FASTA with decoy sequences
  output$downloadDataFASTAdecoy <- downloadHandler(
    filename = function() {
      paste(input$disease, "_with_decoys_", Sys.Date(), ".fasta", sep = "")
    },
    content = function(file) {
      # Start progress bar
      withProgress(message = "Preparing FASTA file...", value = 0, {
        
        # Simulate progress (optional)
        incProgress(0.2, detail = "Filtering data...")
        
        # Get the filtered data
        data <- filteredData()
        
        # Simulate progress
        incProgress(0.4, detail = "Creating sequences...")
        
        # Create headers in OAS/UniProt-like format for original sequences
        headers <- paste("OAS", input$disease, data$BSource, data$BType, data$Isotype, 
                         paste0("Nab_", data$N_antibody), sep = "|")
        
        # Create AAStringSet object for original sequences (for protein sequences)
        fasta_data <- AAStringSet(data$Sequence)
        
        # Create reversed sequences
        reversed_sequences <- reverse(fasta_data)
        
        # Create headers for reversed sequences, prefixed with "rev_"
        reversed_headers <- paste("rev_OAS", input$disease, data$BSource, data$BType, data$Isotype, 
                                  paste0("Nab_", data$N_antibody), sep = "|")
        
        # Combine original and reversed sequences
        all_sequences <- c(fasta_data, reversed_sequences)
        all_headers <- c(headers, reversed_headers)
        
        # Assign headers to the sequence names
        names(all_sequences) <- all_headers
        
        # Simulate progress
        incProgress(0.6, detail = "Writing FASTA file...")
        
        # Write the combined sequences to a temporary FASTA file
        temp_fasta <- tempfile(fileext = ".fasta")
        writeXStringSet(all_sequences, filepath = temp_fasta, format = "fasta")
        
        # Append decoy files to the temp_fasta
        file.append(temp_fasta, "protDB/decoy_UniProt_SP_Human_2024_09_19.fasta")
        file.append(temp_fasta, "protDB/decoy_cRAP.fasta")
        
        # Copy the combined FASTA file to the final output location
        file.copy(temp_fasta, file, overwrite = TRUE)
        
        # Clean up temporary file
        unlink(temp_fasta)
        
        # Clean up memory
        rm(data, fasta_data, reversed_sequences, all_sequences, headers, reversed_headers)
        gc()  # Trigger garbage collection to free up memory
        
        # Complete progress
        incProgress(1, detail = "Done.")
      })
    }
  )
  
  
  
  
  # Info messages using sendSweetAlert
  observeEvent(input$info_antibody, {
    sendSweetAlert(
      session,
      title = "Unique peptide?",
      text = "This option allows you to filter data based on whether the peptide is unique. \n 
      Unique peptide is the peptide appear in only one antibody.",
      type = "info"
    )
  })
  
  observeEvent(input$info_N_antibody, {
    sendSweetAlert(
      session,
      title = "Number of antibodies",
      text = "Filter peptides based on the number of antibodies the peptides appear in. A larger number means these peptides are more popular.",
      type = "info"
    )
  })
  
  observeEvent(input$info_cdr3, {
    sendSweetAlert(
      session,
      title = "Peptides in CDR3 region?",
      text = "This option allows you to filter data based on whether the sequence is present in the CDR3 region.",
      type = "info"
    )
  })
  
  observeEvent(input$info_n_patient, {
    sendSweetAlert(
      session,
      title = "Number of patients",
      text = "Filter peptides based on the number of patients the peptides appear in. A larger number means these peptides are more popular.",
      type = "info"
    )
  })
  
  # Server logic for sendSweetAlert
  observeEvent(input$info_fasta_decoy, {
    sendSweetAlert(
      session = session,
      title = "Download FASTA with decoys",
      text = "This option will add reversed sequences as decoys with the prefix 'rev_' to the sequence names. It will also add contaminants (cRAP database from The Global Proteome Machine) and UniProt proteins (human, reviewed, with isoforms, date: 19/09/2024) to the peptide data.",
      type = "info"
    )
  })
  
  
}

