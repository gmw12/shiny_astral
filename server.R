

server <- function(input, output, session) {
  
  data_dir <<- "/mnt/RawData/4363/OA10034_Astral/tsv_export"
  excel_file_name <<- "df_hela_OA10034.xlsx"
  
  observe({
    if (input$astral_id == 1) {
      data_dir <<- "/mnt/RawData/4363/OA10034_Astral/tsv_export"
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
  

  output$mytable <- renderDT({
    df_rv()
  }, options = list(order = list(list(1, 'desc'))))
  
  
  observe({
    if (input$plot_type == 1) {
      
      output$q4Plot <- renderPlot({create_line_plot(df_rv(), "Sum_Last_Quartile", "Sum of Last Quartile", "Date", "Q4")})
      output$q4_brush_table <- renderTable({create_brush(df_rv(), "Date", "Sum_Last_Quartile", input$q4_brush)})
      
      output$proteinPlot <- renderPlot({create_line_plot(df_rv(), "Proteins", "Proteins", "Date", "Proteins")})
      output$protein_brush_table <- renderTable({create_brush(df_rv(), "Date", "Proteins", input$protein_brush)})  
      
      output$peptidePlot <- renderPlot({create_line_plot(df_rv(), "Peptides", "Peptides", "Date", "Peptides")})
      output$peptide_brush_table <- renderTable({create_brush(df_rv(), "Date", "Peptides", input$peptide_brush)})  
      
      output$precursorPlot <- renderPlot({create_line_plot(df_rv(), "Precursors", "Precursors", "Date", "Precursors")})
      output$precursor_brush_table <- renderTable({create_brush(df_rv(), "Date", "Precursors", input$precursor_brush)})  
      
      output$ideal_peakwidthPlot <- renderPlot({create_line_plot(df_rv(), "Ratio_ideal", "Ratio_ideal", "Date", "Ratio_ideal")})
      output$ideal_peakwidth_brush_table <- renderTable({create_brush(df_rv(), "Date", "Ratio_ideal", input$ideal_peakwidth_brush)})  
      
      output$wide_peakwidthPlot <- renderPlot({create_line_plot(df_rv(), "Ratio_wide", "Ratio_wide", "Date", "Ratio_wide")})
      output$wide_peakwidth_brush_table <- renderTable({create_brush(df_rv(), "Date", "Ratio_wide", input$wide_peakwidth_brush)})  
      
      output$narrow_peakwidthPlot <- renderPlot({create_line_plot(df_rv(), "Ratio_narrow", "Ratio_narrow", "Date", "Ratio_narrow")})
      output$narrow_peakwidth_brush_table <- renderTable({create_brush(df_rv(), "Date", "Ratio_narrow", input$narrow_peakwidth_brush)})  
      
    } else if (input$plot_type == 2) {
      
      output$q4Plot <- renderPlot({create_bar_plot(df_rv(), "Sum_Last_Quartile", "Q4", "Date", "Total")})
      
      output$proteinPlot <- renderPlot({create_bar_plot(df_rv(), "Proteins", "Proteins", "Date", "Total")})
      
      output$peptidePlot <- renderPlot({create_bar_plot(df_rv(), "Peptides", "Peptides", "Date", "Total")})
      
      output$precursorPlot <- renderPlot({create_bar_plot(df_rv(), "Precursors", "Precursors", "Date", "Total")})
      
      output$ideal_peakwidthPlot <- renderPlot({create_bar_plot(df_rv(), "Ratio_ideal", "Peak Width", "Date", "Ratio")})
      
      output$wide_peakwidthPlot <- renderPlot({create_bar_plot(df_rv(), "Ratio_wide", "Peak Width", "Date", "Ratio")})
      
      output$narrow_peakwidthPlot <- renderPlot({create_bar_plot(df_rv(), "Ratio_narrow", "Peak Width", "Date", "Ratio")})
      
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




