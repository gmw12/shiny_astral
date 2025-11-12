server <- function(input, output, session) {
  
  # Define common DT options
  DT_OPTIONS <- list(
    order = list(list(1, 'desc')),
    scrollX = TRUE,
    scrollY = "500px",
    scrollCollapse = TRUE,
    fixedColumns = list(leftColumns = 2)
  )
  
  data_dir <<- "/mnt/RawData/4363/OA10034_Astral/SPD30_tsv"
  excel_file_name <<- "df_hela_OA10034.xlsx"
  
  observe({
    if (input$astral_id == 1) {
      data_dir <<- "/mnt/RawData/4363/OA10034_Astral/SPD30_tsv"
      excel_file_name <<- "df_hela_OA10034.xlsx"
    } else if (input$astral_id == 2){
      data_dir <<- "/mnt/RawData/4363/OA10222_Astral/tsv_export"
      excel_file_name <<- "df_hela_OA10222.xlsx"
    }
  })
  
  
  volumes <- c(dc = data_dir, dd = '/mnt/h_black2', h1 = '/mnt/h_black1', h2 = '/mnt/h_black2', wd = '.', Home = fs::path_home(), getVolumes()())
  shinyFileChoose(input, 'sfb_hela_tsv_file', session = session, roots = volumes, filetypes = c('', 'tsv'))
  
  
  #load hela data from excel
  if (file.exists(str_c(data_dir, "/", excel_file_name))) {
    df_hela <<- readxl::read_excel(str_c(data_dir, "/", excel_file_name))
    df_hela <<- as.data.frame(df_hela)
  } else {
    shinyalert("Error! Need to create starting dataframe manually!", type = "error")
  }
  
  comment_dates <- sort(df_hela$Date, decreasing = TRUE)
  updateSelectInput(session, "comment_date", choices = comment_dates, selected = comment_dates[1])
  
  
  
  
  #reactive value to store dataframe    
  df_rv <- reactiveVal(df_hela)
  
  
  observeEvent(input$sfb_hela_tsv_file, { 
    cat(file = stderr(), "sfb_hela_tsv button clicked...", "\n")
    
    if (is.list(input$sfb_hela_tsv_file)) {
      hela_tsv_sbf <- parseFilePaths(volumes, input$sfb_hela_tsv_file)
      output$hela_tsv_file_name <- renderText({basename(hela_tsv_sbf$datapath)})
      
      extract_data(hela_tsv_sbf$datapath)
      format_save_data()
      updateSelectInput(session, "comment_date", choices = comment_dates, selected = comment_dates[1])
      
      df_rv(df_hela)
    }  
    
  })
  
  
  output$mytable <- renderDT({
    df_rv()
  }, options = DT_OPTIONS, extensions = 'FixedColumns')
  
  
  # Configuration for standard plots
  STANDARD_PLOTS <- list(
    list(output = "q1Plot", brush = "q1_brush", brush_table = "q1_brush_table", 
         column = "Sum_First_Quartile"),
    list(output = "q2Plot", brush = "q2_brush", brush_table = "q2_brush_table", 
         column = "Sum_Second_Quartile"),
    list(output = "q3Plot", brush = "q3_brush", brush_table = "q3_brush_table", 
         column = "Sum_Third_Quartile"),
    list(output = "q4Plot", brush = "q4_brush", brush_table = "q4_brush_table", 
         column = "Sum_Last_Quartile"),
    list(output = "proteinPlot", brush = "protein_brush", brush_table = "protein_brush_table", 
         column = "Proteins"),
    list(output = "peptidePlot", brush = "peptide_brush", brush_table = "peptide_brush_table", 
         column = "Peptides"),
    list(output = "precursorPlot", brush = "precursor_brush", brush_table = "precursor_brush_table", 
         column = "Precursors"),
    list(output = "ideal_peakwidthPlot", brush = "ideal_peakwidth_brush", brush_table = "ideal_peakwidth_brush_table", 
         column = "Ratio_ideal"),
    list(output = "wide_peakwidthPlot", brush = "wide_peakwidth_brush", brush_table = "wide_peakwidth_brush_table", 
         column = "Ratio_wide"),
    list(output = "narrow_peakwidthPlot", brush = "narrow_peakwidth_brush", brush_table = "narrow_peakwidth_brush_table", 
         column = "Ratio_narrow")
  )
  
  # Render standard plots dynamically
  observe({
    for (plot_config in STANDARD_PLOTS) {
      local({
        config <- plot_config
        
        output[[config$output]] <- renderPlot({
          render_dynamic_plot(df_rv(), config$column, input$plot_type)
        })
        
        output[[config$brush_table]] <- renderDT({
          create_brush(df_rv(), "Date", config$column, input[[config$brush]])
        }, options = DT_OPTIONS, extensions = 'FixedColumns')
      })
    }
  })
  
  
  # Get numeric columns for dropdowns
  numeric_columns <- reactive({
    df <- df_rv()
    # Get all column names except Date and Comments
    cols <- colnames(df)
    cols <- cols[!cols %in% c("Date", "Comments")]
    return(cols)
  })
  
  # Update dropdown choices when data changes
  observe({
    cols <- numeric_columns()
    
    # Update Custom Metrics dropdown
    updateSelectInput(session, "custom_metric_select", 
                      choices = cols, 
                      selected = cols[1])
    
    # Update Compare Metrics dropdowns
    if(length(cols) >= 2) {
      updateSelectInput(session, "compare_metric_left", 
                        choices = cols, 
                        selected = cols[1])
      updateSelectInput(session, "compare_metric_right", 
                        choices = cols, 
                        selected = cols[2])
    } else if(length(cols) == 1) {
      updateSelectInput(session, "compare_metric_left", 
                        choices = cols, 
                        selected = cols[1])
      updateSelectInput(session, "compare_metric_right", 
                        choices = cols, 
                        selected = cols[1])
    }
  })
  
  # Custom Metrics tab plot
  output$custom_metric_plot <- renderPlot({
    req(input$custom_metric_select)
    render_dynamic_plot(df_rv(), input$custom_metric_select, input$plot_type)
  })
  
  # Custom Metrics brush table
  output$custom_metric_brush_table <- renderDT({
    req(input$custom_metric_select)
    create_brush(df_rv(), "Date", input$custom_metric_select, input$custom_metric_brush)
  }, options = DT_OPTIONS, extensions = 'FixedColumns')
  
  # Compare Metrics - Left plot
  output$compare_plot_left <- renderPlot({
    req(input$compare_metric_left)
    render_dynamic_plot(df_rv(), input$compare_metric_left, input$plot_type)
  })
  
  # Compare Metrics - Right plot
  output$compare_plot_right <- renderPlot({
    req(input$compare_metric_right)
    render_dynamic_plot(df_rv(), input$compare_metric_right, input$plot_type)
  })
  
  
  observeEvent(input$add_comment, {
    cat(file = stderr(), "add_comment clicked...", "\n") 
    updated_table <- df_rv()
    selected_date <- input$comment_date
    comment_text <- input$comment_text
    updated_table$Comments[updated_table$Date == selected_date] <- comment_text
    df_rv(updated_table)
    df_hela <<- updated_table
    write_xlsx(updated_table, str_c(data_dir, "/", excel_file_name))
  })
  
  
  
}