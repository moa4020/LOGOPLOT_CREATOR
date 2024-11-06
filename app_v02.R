# Load necessary libraries
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("Biostrings")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
  install.packages("ggseqlogo")
}
if (!requireNamespace("msa", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("msa")
}
if (!requireNamespace("shiny", quietly = TRUE)) {
  install.packages("shiny")
}

library(shiny)
library(Biostrings)
library(ggplot2)
library(ggseqlogo)
library(msa)

###############################################################################################
# Define Function Logic
###############################################################################################

# Function to clean the FASTA file and filter out invalid sequences
clean_fasta <- function(fasta_file, seq_type = "DNA") {
  # Load sequences as a BStringSet
  sequences <- readBStringSet(fasta_file)
  
  # Define nonsense sequences to filter out
  nonsense_patterns <- c("#", "/", "?")
  
  # Filter out sequences containing any of the nonsense patterns
  valid_sequences <- lapply(sequences, function(seq) {
    seq_str <- as.character(seq)
    
    # Use vcountPattern to check if any nonsense patterns are present in the sequence
    if (any(sapply(nonsense_patterns, function(pattern) vcountPattern(pattern, seq_str) > 0))) {
      return(NULL) # Return NULL if nonsense pattern is found
    } else {
      return(seq) # Return sequence if it's valid
    }
  })
  
  # Remove NULL entries from the list (sequences containing nonsense patterns)
  valid_sequences <- valid_sequences[!sapply(valid_sequences, is.null)]
  
  # Convert sequences to the specified sequence type
  valid_sequences <- switch(seq_type,
                            "DNA" = as(valid_sequences, "DNAStringSet"),
                            "RNA" = as(valid_sequences, "RNAStringSet"),
                            "AA"  = as(valid_sequences, "AAStringSet"),
                            stop("Invalid sequence type. Choose 'DNA', 'RNA', or 'AA'."))
  
  return(valid_sequences)
}

# Function to load the FASTA file and extract features
load_fasta_with_features <- function(fasta_file, seq_type = "DNA") {
  # Perform cleaning step to remove invalid sequences
  valid_sequences <- clean_fasta(fasta_file, seq_type)
  
  # Initialize vectors to store extracted features
  DonorSubtype <- character(length(valid_sequences))
  DonorID <- character(length(valid_sequences))
  GenomeType <- character(length(valid_sequences))
  GenomeSubtype <- character(length(valid_sequences))
  SeqID <- character(length(valid_sequences))
  SeqData <- character(length(valid_sequences))
  
  # Loop through each sequence and extract features from the header
  for (i in seq_along(valid_sequences)) {
    header <- names(valid_sequences)[i]
    header_parts <- strsplit(header, "\\.")[[1]]
    
    # Check if the header contains exactly 5 parts
    if (length(header_parts) != 5) {
      stop("Each header must contain exactly 5 parts separated by periods ('.'). Format: DonorSubtype.DonorID.GenomeType.GenomeSubtype.SeqID.")
    }
    
    # Assign extracted parts to respective vectors
    DonorSubtype[i] <- header_parts[1]
    DonorID[i] <- header_parts[2]
    GenomeType[i] <- header_parts[3]
    GenomeSubtype[i] <- header_parts[4]
    SeqID[i] <- header_parts[5]
    SeqData[i] <- as.character(valid_sequences[[i]])
  }
  
  # Combine extracted features into a data frame
  fasta_df <- data.frame(
    DonorSubtype = DonorSubtype,
    DonorID = DonorID,
    GenomeType = GenomeType,
    GenomeSubtype = GenomeSubtype,
    SeqID = SeqID,
    Sequence = SeqData,
    stringsAsFactors = FALSE
  )
  
  return(fasta_df)
}

# Function to split fasta file and create logoplots with customization
create_logoplots <- function(fasta_file, 
                             seq_type = "DNA", 
                             split_by = c("DonorSubtype", "DonorID", "GenomeType", "GenomeSubtype"), 
                             filter_by = list(), 
                             output_file = "logoplots.pdf",
                             width = 10,
                             height = 3,
                             annotation_size = 4,
                             output_status) {
  
  output_status("Starting to process the FASTA file...")
  
  # Load sequences using the custom function
  fasta_df <- load_fasta_with_features(fasta_file, seq_type)
  
  # Apply filters based on criteria specified in filter_by
  for (criterion in names(filter_by)) {
    if (criterion %in% colnames(fasta_df)) {
      values <- filter_by[[criterion]]
      fasta_df <- fasta_df[fasta_df[[criterion]] %in% values, ]
    }
  }
  
  # Ensure split_by only contains valid column names
  split_by <- intersect(split_by, c("DonorSubtype", "DonorID", "GenomeType", "GenomeSubtype"))
  if (length(split_by) == 0) {
    stop("Please specify at least one valid column to split by.")
  }
  
  # Generate the splitting factor based on chosen attributes
  split_column <- paste(split_by, collapse = ".")
  fasta_df[[split_column]] <- apply(fasta_df[split_by], 1, paste, collapse = ".")
  
  # Split data by the specified column(s)
  split_groups <- split(fasta_df, fasta_df[[split_column]])
  
  output_status("Creating sequence logos...")
  
  # Create a PDF file to save the plots
  pdf(output_file, width = width, height = height)
  
  # Loop through each group and create a sequence logo plot
  for (group_name in names(split_groups)) {
    group_data <- split_groups[[group_name]]
    sample_size <- nrow(group_data)
    
    # Remove gaps
    group_data$Sequence <- gsub("-", "", group_data$Sequence)
    
    # Convert sequences to a format suitable for ggseqlogo
    seq_list <- BStringSet(group_data$Sequence)
    seq_list <- switch(seq_type,
                       "DNA" = as(seq_list, "DNAStringSet"),
                       "RNA" = as(seq_list, "RNAStringSet"),
                       "AA"  = as(seq_list, "AAStringSet"),
                       stop("Invalid sequence type. Choose 'DNA', 'RNA', or 'AA'."))
    
    if (length(seq_list) > 1) {
      msa_seqs <- msa(seq_list, method = "Muscle")
    } else {
      msa_seqs <- seq_list
    }
    
    msa_seqs <- BStringSet(msa_seqs)
    
    plot <- ggseqlogo(paste0(msa_seqs), method = "prob") +
      ggtitle(paste("Logo Plot for:", group_name)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      xlab(paste("n =", sample_size))
    
    print(plot)
  }
  
  # Close the PDF device
  dev.off()
  
  output_status(paste("Logoplots saved to", output_file))
}

###############################################################################################
# Define UI
###############################################################################################

ui <- fluidPage(
  titlePanel("FASTA Logo Plot Generator"),
  sidebarLayout(
    sidebarPanel(
      fileInput("fasta_file", "Upload FASTA File", accept = c(".fasta", ".fa")),
      selectInput("seq_type", "Sequence Type", choices = c("DNA", "RNA", "AA"), selected = "DNA"),
      checkboxGroupInput("split_by", "Split by", choices = c("DonorSubtype", "DonorID", "GenomeType", "GenomeSubtype"), selected = "DonorID"),
      numericInput("height", "Plot Height", value = 3, min = 1),
      numericInput("width", "Plot Width", value = 10, min = 1),
      numericInput("annotation_size", "Annotation Size", value = 4, min = 1, step = 0.5),
      textInput("output_file", "Output File Name", value = "logoplots.pdf"),
      actionButton("generate_plot", "Generate Logo Plots")
    ),
    mainPanel(
      textOutput("status"),
      downloadButton("download_plot", "Download Logo Plots PDF")
    )
  )
)

###############################################################################################
# Define Server Logic
###############################################################################################

server <- function(input, output, session) {
  
  # Reactive value for status updates
  status_update <- reactiveVal("Waiting for input...")
  
  # Render the status message
  output$status <- renderText({
    status_update()
  })
  
  observeEvent(input$generate_plot, {
    status_update("Starting to process the FASTA file...")
    
    req(input$fasta_file) # Ensure a file is uploaded
    
    # Define parameters for the create_logoplots function
    fasta_file <- input$fasta_file$datapath
    seq_type <- input$seq_type
    split_by <- input$split_by
    height <- input$height
    width <- input$width
    annotation_size <- input$annotation_size
    output_file <- input$output_file
    
    # Check if at least one split option is selected
    if (length(split_by) == 0) {
      output$status <- renderText("Please select at least one split option.")
      return()
    }
    
    tryCatch({
      # Call the function to create logoplots
      create_logoplots(
        fasta_file = fasta_file,
        seq_type = seq_type,
        split_by = split_by,
        filter_by = list(),  # Optionally add filters here
        output_file = output_file,
        height = height,
        width = width,
        annotation_size = annotation_size,
        output_status = status_update  # Pass status update function
      )
      
      status_update("Logoplots generated successfully!")
      
      # Update download button to point to generated file
      output$download_plot <- downloadHandler(
        filename = function() { output_file },
        content = function(file) {
          file.copy(output_file, file)
        }
      )
    }, error = function(e) {
      status_update(paste("Error:", e$message))
    })
  })
}

# Run the app
shinyApp(ui, server)