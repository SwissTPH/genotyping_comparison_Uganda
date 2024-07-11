##########################
# Process datasets and group alleles based on bin sizes then returns each 
# dataset with the newly defined allele lengths
#
# 10.07.2024
# monica.golumbeanu@unibas.ch
##########################

source("~/GitRepos/STPHrepos/genotyping_comparison_Uganda/analysis_scripts/define_alleles_markers/combine_alleles.R")

folder_results = "~/GitRepos/bayesian_analysis_Uganda/datasets_no_additional/"
folder_output = "~/GitRepos/bayesian_analysis_Uganda/datasets_combined_alleles/"
final_table = NULL
list_data = list.files(folder_results, full.names = TRUE)

# Build the long data table with all the alleles
for (f in list_data) {
  marker_name = str_extract(f, "(?<=MSP ).*?(?=_)")
  data_alleles = read_excel(f, sheet = "RecurrentInfections")
  
  data_alleles_long = data_alleles %>% pivot_longer(-c(PatientID, Day, Site), 
                                                    names_to = "allele_name", 
                                                    values_to = "allele_length")
  data_alleles_long$allele_id = substr(data_alleles_long$allele_name, 1, regexpr("_", data_alleles_long$allele_name)-1)
  data_alleles_long$Drug = str_extract(data_alleles_long$PatientID, "(?<=_)[A-Z]+(?=_)")
  data_alleles_long$Third_marker = marker_name
  glurp_idx = which(data_alleles_long$allele_id == "glurp")
  data_alleles_long[glurp_idx, "allele_id"] = marker_name
  data_alleles_long[glurp_idx, "allele_name"] = gsub("glurp", marker_name, data_alleles_long$allele_name[glurp_idx])
  
  final_table = rbind.data.frame(final_table, data_alleles_long)
}

# Merge alleles for each marker based on corresponding bin sizes
list_markers = unique(final_table$allele_id)
merged_allele_table = NULL
for (marker_name in unique(list_markers)) {
  marker_drug_tab = final_table %>% filter(allele_id == marker_name)
  allele_vector = marker_drug_tab %>% select(allele_length)
  bin_size = marker_bins %>% 
              filter(marker_id == marker_name) %>% 
              select(bin_size)

  merged_alleles = merge_alleles(allele_vector = allele_vector$allele_length, 
                                 bin_size = bin_size$bin_size)
  
  marker_drug_tab$final_allele_length = sapply(marker_drug_tab$allele_length, 
                                               function(x) {
                                                 if (is.na(x)) {
                                                   return(NA)
                                                 } else {
                                                   return(merged_alleles[which.min(abs(merged_alleles - x))])
                                                 }
                                               })
  
  merged_allele_table = rbind.data.frame(merged_allele_table, marker_drug_tab)
  
}

# Generate the dataset files
for (third_marker in unique(merged_allele_table$Third_marker)) {
  for (drug in unique(merged_allele_table$Drug)) {
    # Select dataset by drug and marker 
    output_file = paste0(folder_output, "combined_MSP_", third_marker, "_", drug, ".csv")
    merged_allele_table_dataset = merged_allele_table %>% filter(Drug == drug, Third_marker == third_marker)
    
    # Transform back from long format to wide format
    merged_allele_table_dataset$allele_length = merged_allele_table$allele_id = NULL
    output_table = merged_allele_table_dataset %>% 
                      pivot_wider(names_from = allele_name, 
                                  values_from = final_allele_length)

    write.csv2(output_table, output_file, row.names = FALSE, quote = FALSE)
  }
}




