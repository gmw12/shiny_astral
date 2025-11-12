

server <- function(input, output, session) {
  
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
  
  
  #observeEvent(input$add_row, {
  
  observeEvent(input$sfb_hela_tsv_file, { 
    #requires file to work
    #req(input$file)
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
  

  #output$mytable <- renderDT({
  #  df_rv()
  #}, options = list(order = list(list(1, 'desc'))))
  
  
  output$mytable <- renderDT({
    df_rv()
  }, options = list(
    order = list(list(1, 'desc')),
    scrollX = TRUE,
    scrollY = "500px",
    scrollCollapse = TRUE,
    fixedColumns = list(leftColumns = 2)  # Freeze first column
  ), extensions = 'FixedColumns')  # Enable the extension
  
  
  
  observe({
    
    plot_list <- c("q1Plot", "q2Plot", "q3Plot", "q4Plot", 
                   "proteinPlot", "peptidePlot", "precursorPlot", 
                   "ideal_peakwidthPlot", "wide_peakwidthPlot", "narrow_peakwidthPlot")
    
    brush_list <- c("q1_brush", "q2_brush", "q3_brush", "q4_brush", 
                    "protein_brush", "peptide_brush", "precursor_brush", 
                    "ideal_peakwidth_brush", "wide_peakwidth_brush", "narrow_peakwidth_brush")
    
    brush_title_list <- c("q1_brush_table", "q2_brush_table", "q3_brush_table", "q4_brush_table", 
                          "protein_brush_table", "peptide_brush_table", "precursor_brush_table", 
                          "ideal_peakwidth_brush_table", "wide_peakwidth_brush_table", "narrow_peakwidth_brush_table")
    
    column_list <- c("Sum_First_Quartile", "Sum_Second_Quartile", "Sum_Third_Quartile", "Sum_Last_Quartile", 
                     "Proteins", "Peptides", "Precursors", 
                     "Ratio_ideal", "Ratio_wide", "Ratio_narrow")
    
    title_list <- c("Sum of First Quartile", "Sum of Second Quartile", "Sum of Third Quartile", "Sum of Last Quartile", 
                    "Proteins", "Peptides", "Precursors", 
                    "Ratio_ideal", "Ratio_wide", "Ratio_narrow")
    
    x_list <- c("Date", "Date", "Date", "Date", 
                "Date", "Date", "Date", 
                "Date", "Date", "Date")
    
    y_list <- c("Q1", "Q2", "Q3", "Q4", 
                "Proteins", "Peptides", "Precursors", 
                "Ratio_ideal", "Ratio_wide", "Ratio_narrow")
    
    
    if (input$plot_type == 1) {
      
      options_list <- list(
        order = list(list(1, 'desc')),
        scrollX = TRUE,
        scrollY = "500px",
        scrollCollapse = TRUE,
        fixedColumns = list(leftColumns = 2)  # Freeze first column
      )

      for (i in 1:length(plot_list)) {
        local({
          my_i <- i
          output[[plot_list[my_i]]] <- renderPlot({create_line_plot(df_rv(), column_list[my_i], title_list[my_i] , x_list[my_i], y_list[my_i])})
          
          output[[brush_title_list[my_i]]] <- renderDT({
            create_brush(df_rv(), x_list[my_i], column_list[my_i], input[[brush_list[my_i]]])
          }, options = options_list, extensions = 'FixedColumns')
        })
      }
      
      
    } else if (input$plot_type == 2) {
      
      for (i in 1:length(plot_list)) {
        local({
          my_i <- i
          output[[plot_list[my_i]]] <- renderPlot({create_bar_plot(df_rv(), column_list[my_i], title_list[my_i] , x_list[my_i], "Total")})
        })
      }
    }
    
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




