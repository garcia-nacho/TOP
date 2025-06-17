library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(dplyr)
library(tidyr)


input.folder<-"/media/nacho/Data/temp/toptest/TopValidation/top_progress/"



# Function to run shell commands
run_command <- function(command) {
  result <- system(command, intern = TRUE)
  if (length(result) == 0) return("No processes running.")
  return(result)
}

# Function to detect where Nextflow is running
detect_pipeline_disk <- function() {
  path <- system("df -h $(pwd) | awk 'NR==2 {print $6}'", intern = TRUE)
  return(trimws(path))  # Clean whitespace
}


# UI: Shiny Dashboard Layout
ui <- dashboardPage(
  dashboardHeader(title = "TOP Radar"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Processes", tabName = "processes", icon = icon("tasks")),
      menuItem("System Usage", tabName = "system", icon = icon("chart-line")),
      menuItem("Status", tabName = "status", icon = icon("spinner"))
    )
  ),
  dashboardBody(
    tabItems(
      # Tab 1: Process Monitoring
      tabItem(tabName = "processes",
              fluidRow(
                #box(title = "Nextflow Processes", width = 6, DTOutput("nextflow_table")),
                box(title = "SPAdes Processes", width = 6, DTOutput("spades_table")),
                box(title = "Trimmomatic Processes", width = 6, DTOutput("trim_table")),
                box(title = "Mapping Processes", width = 6, DTOutput("bt2_table"))
                
              )
      ),
      
      # Tab 2: CPU & Memory Usage Graphs
      tabItem(tabName = "system",
              fluidRow(
                box(title = "CPU Usage (%)", width = 3, plotlyOutput("cpu_plot")),
                box(title = "Memory Usage (GB)", width = 3, plotlyOutput("memory_plot")),
                box(title = "Disk Usage (GB)", width = 3, plotlyOutput("disk_plot"))
              )
      ),
      
      tabItem(tabName = "status",
              fluidRow(
                box(title = "Status", width = 12, plotlyOutput("status_plot", height = "3600px"))
                
              )
      )
    ),
    
    
    
    
    tags$script("setInterval(function(){Shiny.onInputChange('refresh', Math.random());}, 10000);")
  )
)

# Server: Collects and updates process information
server <- function(input, output, session) {
  
  observeEvent(input$refresh, {
    
    # Get Nextflow and SPAdes processes

    spades_processes <- run_command("ps aux | grep '[s]pades.py'")
    spades_processes<-gsub(".*/bin/spades.py","spades.py",spades_processes)
    
    
    
    trimmomatic_processes <- run_command("ps aux | grep '[t]rimmomatic.jar'")
    trimmomatic_processes<-gsub(".*trimmomatic.jar","trimmomatic.jar",trimmomatic_processes )
    
    mapping_process<- run_command("ps aux | grep '[b]owtie2 -p'")
    mapping_process<-gsub(".*/bin/bowtie2","bowtie2",mapping_process)
    mapping_process<-gsub(".*/bin/bowtie2","bowtie2",mapping_process)
    mapping_process<-gsub("-2 .*/Reverse/","-2 ",mapping_process)
    mapping_process<-gsub("-1 .*/Forward/","-1 ",mapping_process)
    mapping_process<-gsub("-2 .*/Reverse/","-2 ",mapping_process)
    mapping_process<-gsub("--no-unal.*tempRefGenome -1", "-1", mapping_process)
    mapping_process<-gsub("-S temp.sam","",mapping_process)
    
    
    #dfsstat<-read.csv(paste(input.folder,"run_logs.tsv",sep = ""), header = FALSE)
    
    dfsstat<-system( paste("cat ", input.folder, "*", sep = ""), intern = TRUE)
    dfsstat<-base::strsplit(dfsstat,",")
    dfsstat<-as.data.frame(base::do.call(rbind,dfsstat))
    
    
    colnames(dfsstat)<-c("Sample", "Process","Status")
    dfsstat$Sample<-gsub("_.*", "", dfsstat$Sample)
    
    dupind<- which(duplicated(dfsstat))
    if(length(dupind)>0) dfsstat<-dfsstat[-dupind,]
    
    pendingtable<-as.data.frame(expand.grid(unique(dfsstat$Sample), c("Prefilter", "Trimming","KrakenRaw","Assembly","MLST","KrakenClean",
                                                                      "Mapping", "Abricate","NGStar", "Hicap", "Hicgmlst","HinfPBP3","STX_Fastq", "STX_Contigs", "EMMTypper",
                                                                       "MeningoTyper" ,"NGmaster", "Seqsero","Tartrate", "TBPipeline","BPEprofiler" , "Diphtoscan","Seroba",
                                                                      "TBProfiler", "BohlinsECPipeline","BohlinsECPipelineContigs" , "AMRFinderPlus")))
    colnames(pendingtable)<-c("Sample", "Process")
    pendingtable$Status<-"Pending"
    
    if(length(which(paste(pendingtable$Sample, pendingtable$Process) %in% paste(dfsstat$Sample, dfsstat$Process)))>0){
      pendingtable<-pendingtable[-which(paste(pendingtable$Sample, pendingtable$Process) %in% paste(dfsstat$Sample, dfsstat$Process)),]
    }
    dfsstat<-rbind(dfsstat, pendingtable)
    
    dfsstat$Process<-factor(dfsstat$Process, levels= c("Prefilter", "Trimming","KrakenRaw","Assembly","KrakenClean","MLST",
                                                        "Mapping", "Abricate", "AMRFinderPlus","NGStar", "Hicap", "HinfPBP3","Hicgmlst","STX_Fastq", "STX_Contigs", "EMMTypper","Seroba",
                                                        "MeningoTyper" ,"NGmaster", "Seqsero","Tartrate", "TBPipeline","BPEprofiler" , "Diphtoscan",
                                                        "TBProfiler", "BohlinsECPipeline","BohlinsECPipelineContigs"  ))
    
    dfsstat$Status<-factor(dfsstat$Status, levels = c("Completed", "Pending", "Failed", "Skipped"))
    
    #updatetable
    
    # Convert process output to tables
    output$nextflow_table <- renderDT({
      nextflow_df <- data.frame(Process = nextflow_processes, stringsAsFactors = FALSE)
      datatable(nextflow_df, options = list(dom = 't'))
    })
    
    output$spades_table <- renderDT({
      spades_df <- data.frame(Process = spades_processes, stringsAsFactors = FALSE)
      
      datatable(spades_df, options = list(dom = 't'))
    })
    
    
    output$trim_table <- renderDT({
      trimm_df <- data.frame(Process = trimmomatic_processes, stringsAsFactors = FALSE)
      datatable(trimm_df, options = list(dom = 't'))
    })
    
    output$bt2_table <- renderDT({
      bt2_df <- data.frame(Process = mapping_process, stringsAsFactors = FALSE)
      datatable(bt2_df, options = list(dom = 't'))
    })
    
    output$progress_table <- renderDT({
      datatable(dfsstat, options = list(dom = 't'))
    })
    
    
    
    # CPU Usage (Extract CPU %)
    cpu_info <- run_command("top -bn1 | grep 'Cpu(s)' | awk '{print 100 - $8}'")
    cpu_value <- as.numeric(cpu_info)

    
    # Memory Usage (Extract total & used memory)
    mem_info <- run_command("free -m | awk 'NR==2{print $3/1024, $2/1024}'")
    mem_values <- as.numeric(unlist(strsplit(mem_info, " ")))
    
    # Disk Usage (Extract used and total disk space)
    pipeline_disk <- detect_pipeline_disk()
    disk_info <- run_command(paste("df -h ", pipeline_disk, " | awk 'NR==2 {print $3, $2}'",sep = ""))
    
    #disk_info <- run_command("df -h / | awk 'NR==2 {print $3, $2}'")  # Root partition (/)
    disk_values <- as.numeric(unlist(strsplit(gsub("[A-Z]","",gsub(",",".",disk_info)), " ")))  # Convert to numeric values
    
    
    # Handle empty result case
    if (length(disk_values) != 2) {
      disk_values <- c(0, 1)  # Avoid errors if parsing fails
    }
    
    output$disk_plot <- renderPlotly({
      plot_ly(
        x = c("Used", "Free"),
        y = c(disk_values[1], disk_values[2] - disk_values[1]),
        type = "bar",
        marker = list(color = c("red", "green"))
      ) %>%
        layout(title = "Disk Usage", yaxis = list(title = "GB"))
    })
    
    # CPU Usage Plot
    output$cpu_plot <- renderPlotly({
      plot_ly(x = c("Used", "Free"), y = c(cpu_value, 100 - cpu_value), type = "bar", 
              marker = list(color = c("red", "green"))) %>%
        layout(title = "CPU Usage (%)", yaxis = list(title = "Percentage"))
    })
    
    # output$status_plot <- renderPlotly({
    #   status_levels <- c("Completed", "Pending", "Failed", "Skipped")
    #   dfsstat$Status <- factor(dfsstat$Status, levels = status_levels)
    #   
    #   # Create ggplot heatmap
    #   p <- ggplot(dfsstat) +
    #     geom_tile(aes(Process, Sample, fill = Status), color = "black") +
    #     scale_fill_manual(values = c(
    #       "Completed" = "green",
    #       "Pending" = "grey",
    #       "Failed" = "red",
    #       "Skipped" = "blue"
    #     )) +  
    #     guides(fill = guide_legend(override.aes = list(alpha = 1))) +  
    #     theme_minimal()+
    #     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    #   
    #   # Convert ggplot to plotly
    #   ggplotly(p)
    #   
    # })
    
    output$status_plot <- renderPlotly({
      # Suppose you code statuses as integers
      status_map <- c("Completed" = 1, "Pending" = 2, "Failed" = 3, "Skipped" = 4)
      
      dfsstat_num <- dfsstat %>%
        mutate(StatusNum = status_map[as.character(Status)])
      
      # Convert to wide form: each 'Sample' becomes a row, each 'Process' a column
      heat_data <- dfsstat_num %>%
        select(Sample, Process, StatusNum) %>%
        pivot_wider(names_from = Process, values_from = StatusNum)
      

      
      
      # The first column is "Sample", so we remove that from the matrix:
      samples <- heat_data$Sample
      mat <- as.matrix(heat_data[ , c("Prefilter", "Trimming","KrakenRaw","Assembly","KrakenClean","MLST",
                                       "Mapping", "Abricate", "AMRFinderPlus","NGStar", "Hicap", "HinfPBP3","STX_Fastq", "STX_Contigs", "EMMTypper","Seroba",
                                       "MeningoTyper" ,"NGmaster", "Seqsero","Tartrate", "TBPipeline","BPEprofiler" , "Diphtoscan",
                                       "TBProfiler", "BohlinsECPipeline","BohlinsECPipelineContigs"  )])
      
      
      
      # Make a plotly heatmap
      plot_ly(
        x = colnames(mat),
        y = samples,
        z = mat,
        type = "heatmap",
        # Define your colors in the order of your numeric codes (1 -> green, 2 -> gray, etc.)
        colors = c("green", "grey", "red", "blue"),
        showscale = FALSE  # <-- Hide the color scale legend
      ) %>%
        layout(
          # Reverse y-axis so the first sample is at the top
          yaxis = list(autorange = "reversed"),
          autozise=TRUE,
          # Make sure x-axis is treated as categorical
          xaxis = list(type = "category"),
          # Optionally add some margin so axis labels don't get cut off
          margin = list(l = 120, r = 20, b = 100, t = 50)
        )
    })
    
    
    # Memory Usage Plot
    output$memory_plot <- renderPlotly({
      plot_ly(x = c("Used", "Free"), y = c(mem_values[1], mem_values[2] - mem_values[1]), type = "bar",
              marker = list(color = c("blue", "gray"))) %>%
        layout(title = "Memory Usage (GB)", yaxis = list(title = "GB"))
    })
  })
}

# Run the app
shinyApp(ui = ui, server = server)
