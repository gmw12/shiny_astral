library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinydashboard)
library(shinyFiles)
library(shinyalert)

library(ggplot2)
library(plotly)

library(writexl)
library(readxl)

library(data.table)
library(DT)

library(stringr)
library(dplyr)


#  initial_excel("OA10034")
#  initial_excel("OA10222")


initial_excel <- function(astral_id){
  
  data_dir <- str_c("/mnt/RawData/4363/", astral_id, "_Astral/tsv_export")
  
  # get list of files in data_dir
  list_files <- list.files(data_dir, pattern = "\\.tsv$", full.names = TRUE)
  
  df_hela <- data.frame(matrix(ncol=9, nrow=0))
  
  for (file in list_files) {
    #for testing file<-list_files[1]
    df_hela <- extract_data(df_hela, file)
  }
  
  df_hela <<- format_save_data(df_hela)
}



rollup_sum <- function(df){
  cat(file = stderr(), "function rollup_sum...", "\n")
  
  protein_df <- df |> dplyr::group_by(PG.ProteinAccessions) |> dplyr::summarise_all(list(sum))
  protein_df <- data.frame(dplyr::ungroup(protein_df))
  
  #save(df, file="test_df"); save(protein_df, file="test_protein_df")
  #. load(file="test_df"); load(file="test_protein_df")
  
  cat(file = stderr(), "function rollup_sum...end", "\n\n")
  return(protein_df)
}


extract_data <- function(file){
 
  df <- fread(file, header=TRUE)
  
  #remove rows that don't have totalquantity
  df <- df[!is.na(df$`EG.TotalQuantity (Settings)`),]
  
  #extract date from filename
  date <- df$R.FileName[1]
  #find last "_" in date string
  last_underscore <- max(gregexpr("_", date)[[1]])
  date <- substr(date, last_underscore + 1, nchar(date))  # remove ".tsv"
  #reformat date
  date <- paste0(substr(date, 5, 6), substr(date, 1, 2), substr(date, 3, 4)  )
  
  
  #protein, peptide, precursor stats
  protein_count <- length(unique(df$PG.ProteinAccessions))
  peptide_count <- length(unique(df$EG.ModifiedSequence))
  precursor_count <- length(unique(df$EG.PrecursorId))
  
  #sum of last quartile of total quantity
  df$`EG.TotalQuantity (Settings)`<- as.numeric(df$`EG.TotalQuantity (Settings)`)
  q3_range = quantile(df$`EG.TotalQuantity (Settings)`, probs = 0.75, na.rm = TRUE)
  Sum_Last_Quartile <- sum(df$`EG.TotalQuantity (Settings)`[df$`EG.TotalQuantity (Settings)` >= q3_range], na.rm = TRUE)
  Sum_Last_Quartile <- round(Sum_Last_Quartile, digits = 0)
  
  
  # Peak width ratio for "good peaks"
  range_min <- 8/60  
  range_max <- 20/60
  
  # Calculate the ratio
  ratio_pw_ideal <- sum(df$EG.PeakWidth >= range_min & df$EG.PeakWidth <= range_max, na.rm = TRUE) / 
    sum(!is.na(df$EG.PeakWidth))
  
  ratio_pw_ideal <- round((ratio_pw_ideal * 100), digits=1)
  
  
  # Peak width ratio for "wide peaks"
  range_min <- 20/60  
  
  # Calculate the ratio
  ratio_pw_wide <- sum(df$EG.PeakWidth >= range_min, na.rm = TRUE) / sum(!is.na(df$EG.PeakWidth))
  
  ratio_pw_wide <- round((ratio_pw_wide * 100), digits=2)
  
  
  # Peak width ratio for "narrow peaks"
  range_max <- 8/60
  
  # Calculate the ratio
  ratio_pw_narrow <- sum(df$EG.PeakWidth <= range_max, na.rm = TRUE) / sum(!is.na(df$EG.PeakWidth))
  
  ratio_pw_narrow <- round((ratio_pw_narrow * 100), digits=1)
  
  empty_comment <- ""
  
  df_hela <<- rbind(df_hela, c(date, protein_count, peptide_count, precursor_count, Sum_Last_Quartile, ratio_pw_ideal, ratio_pw_wide, ratio_pw_narrow, empty_comment))
  
  return()
}


format_save_data <- function(){

  #set column names
  colnames(df_hela) <- c('Date', "Proteins", "Peptides", "Precursors", "Sum_Last_Quartile", "Ratio_ideal", "Ratio_wide", "Ratio_narrow", "Comments")
  
  #converting all columns to numeric
  df_hela <- df_hela[order(as.Date(df_hela$Date, format = "%y%m%d")), ]
  df_hela$Sum_Last_Quartile <- as.numeric(df_hela$Sum_Last_Quartile)
  df_hela$Proteins <- as.numeric(df_hela$Proteins)
  df_hela$Precursors <- as.numeric(df_hela$Precursors)
  df_hela$Peptides <- as.numeric(df_hela$Peptides)
  df_hela$Ratio_ideal <- as.numeric(df_hela$Ratio_ideal)
  df_hela$Ratio_wide <- as.numeric(df_hela$Ratio_wide)
  df_hela$Ratio_narrow <- as.numeric(df_hela$Ratio_narrow)
  
  df_hela <<- df_hela |> distinct()
   
  #save to excel
  write_xlsx(df_hela, str_c(data_dir, "/", excel_file_name))
  
  return()
}


create_line_plot <- function(df, data_column, title, x_lab, y_lab) {
  # Add a column to indicate if there are comments
  df$has_comment <- !is.na(df$Comments) & df$Comments != ""
  
  g <- ggplot(df, aes(x = Date, y = .data[[data_column]])) +
    geom_line(group = 1, color = "blue", linewidth = 1.2) +
    geom_point(aes(shape = has_comment, color = has_comment), size = 4) + 
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) +  # circle vs triangle
    scale_color_manual(values = c("FALSE" = "red", "TRUE" = "green")) +
    coord_cartesian(ylim = c(0, max(df[[data_column]], na.rm = TRUE) * 1.05)) +
    labs(title = title, x = NULL, y = y_lab) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(5.5, 5.5, 20, 5.5, "pt"),
          legend.position = "none")  # No legend
  
  return(g)
}


create_brush <- function(df, x_var, y_var, input_brush) {
  n = nrow(brushedPoints(df, xvar= x_var, yvar=y_var, brush = input_brush))
  if(n==0)
    return()
  else
    return(brushedPoints(df, xvar= x_var, yvar=y_var, brush = input_brush))
}


create_bar_plot <- function(df, data_column, title, x_lab, y_lab) {
  g <- ggplot(df, aes(x = Date, y = .data[[data_column]], fill = Date)) +
    stat_summary(fun = mean, geom = "bar") +
    coord_cartesian(ylim = c(0, (max(df[[data_column]]) * 1.1))) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          plot.title = element_text(color = "blue", size = 16, vjust = 0.5, hjust = 0.5),
          legend.position = "none",
          plot.margin = margin(5.5, 5.5, 20, 5.5, "pt")) +  # Added bottom margin
    ggtitle(title)
  
  return(g)
}



