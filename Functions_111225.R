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

# astral_id = "OA10034"; data_dir <- "/mnt/RawData/4363/OA10034_Astral/SPD30_tsv"; excel_file_name <- "df_hela_OA10034.xlsx"
#  df_hela <- initial_excel("OA10034")
#  initial_excel("OA10222")


initial_excel <- function(astral_id){
  
  data_dir <- str_c("/mnt/RawData/4363/", astral_id, "_Astral/SPD30_tsv")
  
  # get list of files in data_dir
  list_files <- list.files(data_dir, pattern = "\\.tsv$", full.names = TRUE)
  
  df_hela <- data.frame(matrix(ncol=51, nrow=0))
  
  for (file in list_files) {
    #for testing file<-list_files[1]
    df_hela <- extract_data(df_hela, file)
  }
  
  test_df_hela <<- df_hela
  format_save_data(df_hela)
  
  return(df_hela)
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


extract_data <- function(df_hela, file){
  
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
  q1_range = quantile(df$`EG.TotalQuantity (Settings)`, probs = 0.25, na.rm = TRUE)
  q2_range = quantile(df$`EG.TotalQuantity (Settings)`, probs = 0.50, na.rm = TRUE)
  q3_range = quantile(df$`EG.TotalQuantity (Settings)`, probs = 0.75, na.rm = TRUE)
  Sum_First_Quartile <- sum(df$`EG.TotalQuantity (Settings)`[df$`EG.TotalQuantity (Settings)` <= q1_range], na.rm = TRUE)
  Sum_First_Quartile <- round(Sum_First_Quartile, digits = 0)
  Sum_Second_Quartile <- sum(df$`EG.TotalQuantity (Settings)`[df$`EG.TotalQuantity (Settings)` >= q1_range & df$`EG.TotalQuantity (Settings)` <= q2_range ], na.rm = TRUE)
  Sum_Second_Quartile <- round(Sum_Second_Quartile, digits = 0)
  Sum_Third_Quartile <- sum(df$`EG.TotalQuantity (Settings)`[df$`EG.TotalQuantity (Settings)` >= q2_range & df$`EG.TotalQuantity (Settings)` <= q3_range ], na.rm = TRUE)
  Sum_Third_Quartile <- round(Sum_Third_Quartile, digits = 0)
  Sum_Last_Quartile <- sum(df$`EG.TotalQuantity (Settings)`[df$`EG.TotalQuantity (Settings)` >= q3_range], na.rm = TRUE)
  Sum_Last_Quartile <- round(Sum_Last_Quartile, digits = 0)
  
  # Peak width ratio for "good peaks"
  range_min <- 8/60  
  range_max <- 30/60
  
  # Calculate the ratio
  ratio_pw_ideal <- sum(df$EG.PeakWidth >= range_min & df$EG.PeakWidth <= range_max, na.rm = TRUE) / 
    sum(!is.na(df$EG.PeakWidth))
  
  ratio_pw_ideal <- round((ratio_pw_ideal * 100), digits=1)
  
  
  # Peak width ratio for "wide peaks"
  range_min <- 30/60  
  
  # Calculate the ratio
  ratio_pw_wide <- sum(df$EG.PeakWidth >= range_min, na.rm = TRUE) / sum(!is.na(df$EG.PeakWidth))
  
  ratio_pw_wide <- round((ratio_pw_wide * 100), digits=2)
  
  
  # Peak width ratio for "narrow peaks"
  range_max <- 8/60
  
  # Calculate the ratio
  ratio_pw_narrow <- sum(df$EG.PeakWidth <= range_max, na.rm = TRUE) / sum(!is.na(df$EG.PeakWidth))
  
  ratio_pw_narrow <- round((ratio_pw_narrow * 100), digits=1)
  
  #calculate median
  quant_median = median(df$`EG.TotalQuantity (Settings)`, na.rm = TRUE)
  quant_median <- round(quant_median, digits = 0)
  
  
  #count values in df$`EG.TotalQuantity (Settings)`
  count_values_0_50 <- sum(df$`EG.TotalQuantity (Settings)` >= 0 & df$`EG.TotalQuantity (Settings)` <= 50, na.rm = TRUE)
  count_values_50_100 <- sum(df$`EG.TotalQuantity (Settings)` >= 50 & df$`EG.TotalQuantity (Settings)` <= 100, na.rm = TRUE)
  count_values_100_200 <- sum(df$`EG.TotalQuantity (Settings)` >= 100 & df$`EG.TotalQuantity (Settings)` <= 200, na.rm = TRUE)
  count_values_200_500 <- sum(df$`EG.TotalQuantity (Settings)` >= 200 & df$`EG.TotalQuantity (Settings)` <= 500, na.rm = TRUE)
  count_values_500_1000 <- sum(df$`EG.TotalQuantity (Settings)` >= 500 & df$`EG.TotalQuantity (Settings)` <= 1000, na.rm = TRUE)
  count_values_1000_2000 <- sum(df$`EG.TotalQuantity (Settings)` >= 1000 & df$`EG.TotalQuantity (Settings)` <= 2000, na.rm = TRUE)
  count_values_2000_5000 <- sum(df$`EG.TotalQuantity (Settings)` >= 2000 & df$`EG.TotalQuantity (Settings)` <= 5000, na.rm = TRUE)
  count_values_5000_10000 <- sum(df$`EG.TotalQuantity (Settings)` >= 5000 & df$`EG.TotalQuantity (Settings)` <= 10000, na.rm = TRUE)
  count_values_10000 <- sum(df$`EG.TotalQuantity (Settings)` >= 10000, na.rm = TRUE)
  
  
  #quartile values
  quartile_1 = round(quantile(df$`EG.TotalQuantity (Settings)`, probs = 0.25, na.rm = TRUE), digits = 0)
  quartile_2 = round(quantile(df$`EG.TotalQuantity (Settings)`, probs = 0.50, na.rm = TRUE), digits = 0)
  quartile_3 = round(quantile(df$`EG.TotalQuantity (Settings)`, probs = 0.75, na.rm = TRUE), digits = 0)
  
  
  #specific proteins
  protein_list <- c("O00410", "O00571", "A3KMH1", "A5YKK6", "A6NHR9")
  
  protein_O00410_count <- sum(grepl("O00410", df$PG.ProteinAccessions))
  protein_O00571_count <- sum(grepl("O00571", df$PG.ProteinAccessions))
  protein_A3KMH1_count <- sum(grepl("A3KMH1", df$PG.ProteinAccessions))
  protein_A5YKK6_count <- sum(grepl("A5YKK6", df$PG.ProteinAccessions))
  protein_A6NHR9_count <- sum(grepl("A6NHR9", df$PG.ProteinAccessions))
  
  protein_O00410 <- round(sum(df[df$PG.ProteinAccessions=="O00410", ]$`EG.TotalQuantity (Settings)`), digits = 0)
  protein_O00571 <- round(sum(df[df$PG.ProteinAccessions=="O00571", ]$`EG.TotalQuantity (Settings)`), digits = 0)
  protein_A3KMH1 <- round(sum(df[df$PG.ProteinAccessions=="A3KMH1", ]$`EG.TotalQuantity (Settings)`), digits = 0)
  protein_A5YKK6 <- round(sum(df[df$PG.ProteinAccessions=="A5YKK6", ]$`EG.TotalQuantity (Settings)`), digits = 0)
  protein_A6NHR9 <- round(sum(df[df$PG.ProteinAccessions=="A6NHR9", ]$`EG.TotalQuantity (Settings)`), digits = 0)
  
  ratio_O00410_A3KMH1 <- round(protein_O00410 / protein_A3KMH1, digits = 2)
  ratio_O00571_A6NHR9 <- round(protein_O00571 / protein_A6NHR9, digits = 2)
  ratio_A5YKK6_A6NHR9 <- round(protein_A5YKK6 / protein_A6NHR9, digits = 2)
  
  #noise
  sn_median = median(df$EG.SignalToNoise, na.rm = TRUE)
  sn_median <- round(sn_median, digits = 0)
  
  sn_count_values_0_5 <- sum(df$EG.SignalToNoise >= 0 & df$EG.SignalToNoise <= 5, na.rm = TRUE)
  sn_count_values_5_10 <- sum(df$EG.SignalToNoise >= 5 & df$EG.SignalToNoise <= 10, na.rm = TRUE)
  sn_count_values_10_20 <- sum(df$EG.SignalToNoise >= 10 & df$EG.SignalToNoise <= 20, na.rm = TRUE)
  sn_count_values_20_50 <- sum(df$EG.SignalToNoise >= 20 & df$EG.SignalToNoise <= 50, na.rm = TRUE)
  sn_count_values_50 <- sum(df$EG.SignalToNoise >= 50, na.rm = TRUE)
  
  
  df$EG.SignalToNoise <- as.numeric(df$EG.SignalToNoise)
  sn_q1_range <- quantile(df$EG.SignalToNoise, probs = 0.25, na.rm = TRUE)[[1]]
  sn_q1_range <- round(sn_q1_range, digits = 1)
  sn_q2_range <- quantile(df$EG.SignalToNoise, probs = 0.50, na.rm = TRUE)[[1]]
  sn_q2_range <- round(sn_q2_range, digits = 1)
  sn_q3_range <- quantile(df$EG.SignalToNoise, probs = 0.75, na.rm = TRUE)[[1]]
  sn_q3_range <- round(sn_q3_range, digits = 1)
  sn_Sum_First_Quartile <- sum(df$EG.SignalToNoise[df$EG.SignalToNoise <= sn_q1_range], na.rm = TRUE)
  sn_Sum_First_Quartile <- round(sn_Sum_First_Quartile, digits = 0)
  sn_Sum_Second_Quartile <- sum(df$EG.SignalToNoise[df$EG.SignalToNoise >= sn_q1_range & df$EG.SignalToNoise <= sn_q2_range ], na.rm = TRUE)
  sn_Sum_Second_Quartile <- round(sn_Sum_Second_Quartile, digits = 0)
  sn_Sum_Third_Quartile <- sum(df$EG.SignalToNoise[df$EG.SignalToNoise >= sn_q2_range & df$EG.SignalToNoise <= sn_q3_range ], na.rm = TRUE)
  sn_Sum_Third_Quartile <- round(sn_Sum_Third_Quartile, digits = 0)
  sn_Sum_Last_Quartile <- sum(df$EG.SignalToNoise[df$EG.SignalToNoise >= sn_q3_range], na.rm = TRUE)
  sn_Sum_Last_Quartile <- round(sn_Sum_Last_Quartile, digits = 0)
  
  
  #add empty comment
  empty_comment <- ""
  
  df_hela <- rbind(df_hela, c(date, protein_count, peptide_count, precursor_count, Sum_First_Quartile, 
                              Sum_Second_Quartile, Sum_Third_Quartile, Sum_Last_Quartile, ratio_pw_ideal, 
                              ratio_pw_wide, ratio_pw_narrow, quant_median, count_values_0_50, count_values_50_100, 
                              count_values_100_200, count_values_200_500, count_values_500_1000, count_values_1000_2000, 
                              count_values_2000_5000, count_values_5000_10000, count_values_10000,
                              quartile_1, quartile_2, quartile_3, 
                              protein_O00410_count, protein_O00571_count, protein_A3KMH1_count, protein_A5YKK6_count, protein_A6NHR9_count,
                              protein_O00410, protein_O00571, protein_A3KMH1, protein_A5YKK6, protein_A6NHR9,
                              ratio_O00410_A3KMH1, ratio_O00571_A6NHR9, ratio_A5YKK6_A6NHR9,
                              sn_median, sn_count_values_0_5, sn_count_values_5_10, sn_count_values_10_20, sn_count_values_20_50, sn_count_values_50,
                              sn_q1_range, sn_q2_range, sn_q3_range,
                              sn_Sum_First_Quartile, sn_Sum_Second_Quartile, sn_Sum_Third_Quartile, sn_Sum_Last_Quartile,
                              empty_comment))
  
  #set column names
  colnames(df_hela) <- c('Date', "Proteins", "Peptides", "Precursors", "Sum_First_Quartile", 
                         "Sum_Second_Quartile", "Sum_Third_Quartile", "Sum_Last_Quartile", "Ratio_ideal", 
                         "Ratio_wide", "Ratio_narrow", "Median", "Count_0_50", "Count_50_100", 
                         "Count_100_200", "Count_200_500", "Count_500_1000", "Count_1000_2000", 
                         "Count_2000_5000", "Count_5000_10000", "Count_10000",
                         "Quartile1", "Quartile2", "Quartile3", 
                         "O00410_Count", "O00571_Count", "A3KMH1_Count", "A5YKK6_Count", "A6NHR9_Count",
                         "O00410_Quantity", "O00571_Quantity", "A3KMH1_Quantity", "A5YKK6_Quantity", "A6NHR9_Quantity",
                         "Ratio_O00410_A3KMH1",  "Ratio_O00571_A6NHR9", "Ratio_A5YKK6_A6NHR9",
                         "SN_Median", "SN_Count_0_5", "SN_Count_5_10", "SN_Count_10_20", "SN_Count_20_50", "SN_Count_50",
                         "SN_Quartile1", "SN_Quartile2", "SN_Quartile3",
                         "SN_Sum_First_Quartile", "SN_Sum_Second_Quartile", "SN_Sum_Third_Quartile", "SN_Sum_Last_Quartile",
                         "Comments")
  
  #converting all columns to numeric
  df_hela <- df_hela[order(as.Date(df_hela$Date, format = "%y%m%d"), decreasing = TRUE), ]
  
  #set all columns in df_hela to numeric except Date and Comments
  df_hela <- df_hela %>%
    mutate(across(-c(Date, Comments), as.numeric))
  
  df_hela <- df_hela |> distinct()
  
  return(df_hela)
}



format_save_data <- function(df_hela){
  
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



