#################################
# plot the distributions of the probabilities of recrudescence
#
# 09.07.2024
# monica.golumbeanu@unibas.ch
#################################

# Read in the data and build the data frame
results_folder = "~/genotyping/analysis_Uganda/results_test/"
list_results = list.dirs(path = results_folder, full.names = TRUE)

for (folder_name in list_results) {
  marker_name = str_replace(folder_name, "MSP ", "")
  marker_name = str_replace(marker_name, "_AL", "")
  marker_name = str_replace(marker_name, "_DP", "")
  
  drug = 
}
