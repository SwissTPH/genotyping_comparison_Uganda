#################################
# Script to generate the datasets from the .xls file
#
# monica.golumbeanu@swisstph.ch
#################################

library(readxl)
library(xlsx)
library(dplyr)
library(stringr)
library(purrr)

# Functions for writing to excel files
.write_block <- function(wb, sheet, y, rowIndex=seq_len(nrow(y)),
                         colIndex=seq_len(ncol(y)), showNA=TRUE)
{
  rows  <- createRow(sheet, rowIndex)      # create rows
  cells <- createCell(rows, colIndex)      # create cells
  
  for (ic in seq_len(ncol(y)))
    mapply(setCellValue, cells[seq_len(nrow(cells)), colIndex[ic]], y[,ic], FALSE, showNA)
  
  # Date and POSIXct classes need to be formatted
  indDT <- which(sapply(y, function(x) inherits(x, "Date")))
  if (length(indDT) > 0) {
    dateFormat <- CellStyle(wb) + DataFormat(getOption("xlsx.date.format"))
    for (ic in indDT){
      lapply(cells[seq_len(nrow(cells)),colIndex[ic]], setCellStyle, dateFormat)
    }
  }
  
  indDT <- which(sapply(y, function(x) inherits(x, "POSIXct")))
  if (length(indDT) > 0) {
    datetimeFormat <- CellStyle(wb) + DataFormat(getOption("xlsx.datetime.format"))
    for (ic in indDT){
      lapply(cells[seq_len(nrow(cells)),colIndex[ic]], setCellStyle, datetimeFormat)
    }
  }
  
}

write.xlsx.custom <- function(x, file, sheetName="Sheet1",
                              col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE)
{
  if (!is.data.frame(x))
    x <- data.frame(x)    # just because the error message is too ugly
  
  iOffset <- jOffset <- 0
  if (col.names)
    iOffset <- 1
  if (row.names)
    jOffset <- 1
  
  if (append && file.exists(file)){
    wb <- loadWorkbook(file)
  } else {
    ext <- gsub(".*\\.(.*)$", "\\1", basename(file))
    wb  <- createWorkbook(type=ext)
  }  
  sheet <- createSheet(wb, sheetName)
  
  noRows <- nrow(x) + iOffset
  noCols <- ncol(x) + jOffset
  if (col.names){
    rows  <- createRow(sheet, 1)                  # create top row
    cells <- createCell(rows, colIndex=1:noCols)  # create cells
    mapply(setCellValue, cells[1,(1+jOffset):noCols], colnames(x))
  }
  if (row.names)             # add rownames to data x                   
    x <- cbind(rownames=rownames(x), x)
  
  if(nrow(x) > 0) {
    colIndex <- seq_len(ncol(x))
    rowIndex <- seq_len(nrow(x)) + iOffset
    
    .write_block(wb, sheet, x, rowIndex, colIndex, showNA)
  }
  saveWorkbook(wb, file)
  
  invisible()
}

# Specification of the files, to adapt according to user system
reinfections_file = "~/GitRepos/bayesian_analysis_Uganda/data/Reinfections.xls"
day_0_file = "~/GitRepos/bayesian_analysis_Uganda/data/Day 0 Data set.xls"
data_folder = "~/GitRepos/bayesian_analysis_Uganda/datasets_no_additional/"

# Checking if the data contains DAY 0 samples
# day_0_data = readxl::read_xls(day_0_file, sheet = "Day 0 MSP data", 
#                               col_names = TRUE)
# part1_name = substr(day_0_data$`Lab Sample ID`, 1, str_locate(day_0_data$`Lab Sample ID`, "D")-1)
# part2_name = substr(day_0_data$`Lab Sample ID`, str_locate(day_0_data$`Lab Sample ID`, "_"), length(day_0_data$`Lab Sample ID`))
# day_0_data$`Lab Sample ID` = paste0(part1_name, part2_name)
# part1_name = substr(day_0_data$`PATIENT ID`, 1, str_locate(day_0_data$`PATIENT ID`, "D")-1)
# part2_name = substr(day_0_data$`PATIENT ID`, str_locate(day_0_data$`PATIENT ID`, "_"), length(day_0_data$`PATIENT ID`))
# day_0_data$`PATIENT ID` = paste0(part1_name, part2_name)

# Extract names of sheets from the excel data file
sheets_all_data = excel_sheets(reinfections_file)

# Extract each sheet in a separate dataset
for (sheet_name in sheets_all_data) {
  dataset_sheet = readxl::read_xls(reinfections_file, sheet = sheet_name, 
                                   col_names = TRUE)
  # Rename the columns
  dataset_sheet = dataset_sheet %>% rename(PatientID = `PATIENT ID`,
                                           Day = DAY,
                                           Site = SITE,
                                           K1_1 = KI_1,
                                           R033_1 = RO33_1,
                                           R033_2 = RO33_2)
  colnames(dataset_sheet) = str_replace(colnames(dataset_sheet), "PFP2_", "PFPK2_")
  colnames(dataset_sheet) = str_replace(colnames(dataset_sheet), "FC 27_", "FC27_")
  colnames(dataset_sheet) = str_replace(colnames(dataset_sheet), "Poly_", "Poly a_")
  
  # Rename the marker columns with "glurp"
  marker_name = substr(sheet_name, str_locate(sheet_name, "MSP")[2]+2, nchar(sheet_name))
  marker_name = paste0(marker_name, "_")
  if(marker_name == "POLY a_") {
    marker_name = "Poly a_"
  }
  print(paste("Marker:", marker_name))
  print("Before replacement:")
  print(colnames(dataset_sheet))
  colnames(dataset_sheet) = str_replace(colnames(dataset_sheet), marker_name, "glurp_")
  colnames(dataset_sheet) = str_replace(colnames(dataset_sheet), "IC3D7_", "3D7_")
  print("After replacement:")
  print(colnames(dataset_sheet))
  
  # Remove the Days names from the IDs
  part1_name = substr(dataset_sheet$PatientID, 1, str_locate(dataset_sheet$PatientID, "D")-1)
  part2_name = substr(dataset_sheet$PatientID, str_locate(dataset_sheet$PatientID, "_"), length(dataset_sheet$PatientID))
  dataset_sheet$PatientID = paste0(part1_name, part2_name)
  
  # Separate the AL and DP samples
  dataset_sheet_AL = dataset_sheet %>% filter(grepl("_AL_", PatientID))
  dataset_sheet_AL[, 4:ncol(dataset_sheet_AL)] = as.data.frame(lapply(dataset_sheet_AL[, 4:ncol(dataset_sheet_AL)], function(x) as.numeric(as.character(x))))
  dataset_sheet_DP = dataset_sheet %>% filter(grepl("_DP_", PatientID))
  dataset_sheet_DP[, 4:ncol(dataset_sheet_DP)] = as.data.frame(lapply(dataset_sheet_DP[, 4:ncol(dataset_sheet_DP)], function(x) as.numeric(as.character(x))))
  
  # Write the two sheets to the file for each dataset
  AL_file = paste0(data_folder, sheet_name, "_AL", ".xlsx")
  DP_file = paste0(data_folder, sheet_name, "_DP", ".xlsx")
  write.xlsx.custom(as.data.frame(dataset_sheet_AL), file = AL_file, 
                    sheetName="RecurrentInfections", row.names = FALSE)
  write.xlsx.custom(as.data.frame(dataset_sheet_DP), file = DP_file, 
             sheetName="RecurrentInfections", row.names = FALSE)
  
  # Build the Additional sheet per drug
  # Split the data frame into groups based on the "site" column
  grouped_df_AL = dataset_sheet_AL %>% group_split(Site)
  grouped_df_DP = dataset_sheet_DP %>% group_split(Site)
  
  # Select the first three rows from each site and combine them into a new data frame
  selected_rows_AL = map(grouped_df_AL, ~slice(.x, 1:3)) %>% bind_rows()
  selected_rows_AL$Day = 0
  selected_rows_AL$PatientID = c(1:nrow(selected_rows_AL))
  selected_rows_DP = map(grouped_df_DP, ~slice(.x, 1:3)) %>% bind_rows()
  selected_rows_DP$Day = 0
  selected_rows_DP$PatientID = c(1:nrow(selected_rows_DP))
  
  selected_rows_AL = selected_rows_AL[FALSE, ]
  selected_rows_DP = selected_rows_DP[FALSE, ]
  
  write.xlsx.custom(as.data.frame(selected_rows_AL), file = AL_file, sheetName="Additional", 
                    row.names = FALSE, append = TRUE)
  write.xlsx.custom(as.data.frame(selected_rows_DP), file = DP_file, sheetName="Additional", 
                    row.names = FALSE, append = TRUE)
}

