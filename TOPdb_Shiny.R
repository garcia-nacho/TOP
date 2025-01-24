
#Library loading
library(shiny)
library(ggplot2)


#Data loading
pdf(NULL)

df<-read.csv("/media/nacho/Data/OnGoingProjects/TOPDB/db/TOP_Database_24012025.csv")

#Fix emm
df$emm.type<-gsub("EMM","",df$emm.type)

df[ which(is.na(df),arr.ind = TRUE)]<-"ND"
df$count<-1
sps<-unique(df$rMLST_taxon)
sps<-sps[order(sps)]

options(shiny.sanitize.errors = FALSE)
# Define UI
ui <- fluidPage(
  
  # Application title
  titlePanel("TOP DB"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("agent", "Agent", c("All", unique(sps))),
      selectInput("type", "Plot Type",c("BoxPlot", "LinePlot", "BarPlot")),
      selectInput("x", "X-axis", c("None",colnames(df))),
      selectInput("y", "Y-axis",c("None",colnames(df))),
      checkboxInput("aggregate", "Aggregate Data for X", value = FALSE),
      selectInput("col", "Colour",c("None",colnames(df))),
      actionButton("do", "Plot", align = "center"),
      width = 2
    ),
    
    mainPanel(
      plotOutput("top_plot")
      
    )
  )
)

# Define server 
server <- function(input, output, session) {
  
  filtered_data <- reactive({
    if(input$agent!="All"){
      df[df$rMLST_taxon == input$agent, ]  
    }else{
      df
    }
    
  })

  observe({
    # Get unique subcategories from the filtered data
    subcategories <- colnames(filtered_data)
    
    # Update the SubCategory input based on the filtered data
    updateSelectInput(
      session,
      inputId = "x",
      choices = subcategories,
      selected = subcategories[1]  # Default to the first option
    )
  })
  
  
  
  plot_data <- eventReactive(input$do, {
    data <- filtered_data()

    # Perform aggregation if required
    if (input$aggregate && input$x != "None") {
      group_var <- input$x
      
      if(input$col != "None") group_var<-c(group_var, input$col)

      # Ensure Y is numeric for aggregation
      if (input$y != "None" && is.numeric(data[[input$y]])) {
        # data <- data %>%
        #   group_by(across(all_of(group_var))) %>%
        #   summarize(
        #     agg_y = mean(.data[[input$y]], na.rm = TRUE), .groups = "drop"
        #   )
        
       data<- as.data.frame(aggregate(as.formula(paste( "count ~", paste(group_var, collapse = " + ") )), data, length))
        
      } else {
        # If no Y is provided or Y is not numeric, count occurrences
        data<- as.data.frame( aggregate(as.formula(paste( "count ~", paste(group_var, collapse = " + ") )), data, length))

      }
    }
    
    data
  })
  
  
  output$top_plot <- renderPlot({
    # Further subset data by selected SubCategory
   
    
    selected_x <- input$x
    selected_y <- if (input$aggregate) colnames(plot_data())[2] else input$y
    if (input$aggregate) selected_y <- "count"
    selected_col <- input$col
    faceting<-input$facet
    
    print(colnames(plot_data()))
    
    if (selected_x == "None" || selected_y == "None") {
      return(NULL)  # Don't render plot if X or Y are not selected
    }
    
    # Get selected columns for the plot

    p <- ggplot(plot_data(), aes_string(x = selected_x, y = selected_y, color = if (selected_col != "None") selected_col else NULL,
                                        fill = if (selected_col != "None") selected_col else NULL))
    if (input$type == "BoxPlot") {
      p <- p + geom_boxplot()
    } else if (input$type == "LinePlot") {
      p <- p + geom_line()
    } else if (input$type == "BarPlot") {
      p <- p + geom_bar(stat = "identity", position = "stack")
      
      
      #p<- p + geom_col(na.rm = TRUE)
    }

    
    # Add theme
    p  + theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
  })
  


  
}
# Run the application 
shinyApp(ui = ui, server = server)