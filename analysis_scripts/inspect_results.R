#################################
# Inspect results
# created 26.02.2024
# monica.golumbeanu@swisstph.ch
################################

list_datasets = list.files("~/genotyping/analysis_Uganda/results/", full.names = TRUE)
thresh = 0.8

res = NULL
table_all_prob_AL = NULL
table_all_prob_DP = NULL
for (f in list_datasets) {
  print(f)
  proba_recr = read.csv(paste0(f, "/", "probability_of_recrudescence_ALL.csv"), skip = 1, header = FALSE)
  
  hist(proba_recr$V3)
  
  nb_r = as.integer(sum(as.double(proba_recr$V3) > thresh))
  print(nb_r)
  row_res = cbind.data.frame(f, nb_r)
  res = rbind.data.frame(res, row_res)
}
colnames(res) = c("file_name", "n_recrudescences")
