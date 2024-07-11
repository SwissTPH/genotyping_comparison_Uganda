#################################
# Script to calculate the plots for the presentation
# 
#################################

library(readxl)
library(dplyr)
library(tidyr)
library(plotly)
library(ggplot2)
library(tidyverse)
library(ggpubr)

source("~/GitRepos/STPHrepos/genotyping_comparison_Uganda/analysis_scripts/define_alleles_markers/combine_alleles.R")

# Function for processing the data for plotting
process_plot_data = function(plot_data_df) {
  # Make sure that the colnames are the ones expected later for plotting
  colnames(plot_data_df) = c("Haplotype", "Amount", "Frequency")
  
  # Remove rows with NA size
  plot_data_df = plot_data_df[!is.na(plot_data_df$Haplotype), c("Haplotype", "Amount", "Frequency")]
  
  # Process the data frame, define labels
  plot_data_df$Haplotype = factor(plot_data_df$Haplotype, levels = sort(plot_data_df$Haplotype))
  # plot_data_df$Label = paste0("n = ", plot_data_df$Amount, "\n(", round(plot_data_df$Amount/sum(plot_data_df$Amount), digits = 2)*100,"%)")
  plot_data_df$Label = paste0(round(plot_data_df$Amount/sum(plot_data_df$Amount), digits = 2)*100, "%")
  
  plot_data_df = plot_data_df[order(plot_data_df$Amount, decreasing = TRUE), ]
  
  # Create sector labels
  plot_data_df$pct = as.double(round(plot_data_df$Amount/sum(plot_data_df$Amount), 2) * 100)
  plot_data_df[4:nrow(plot_data_df), "Label"] = ""
  # plot_data_df[4:nrow(plot_data_df), "pct"] = ""
  # 
  # plot_data_df$Label[which(plot_data_df$pct <= 4)] = 0 
  # plot_data_df$pct[which(plot_data_df$pct <= 4)] = 0  # Anything less than 4% should be blank
  # plot_data_df$pct = paste0(plot_data_df$pct, "%")
  # plot_data_df$pct[which(plot_data_df$pct == "0%")] = ""
  # plot_data_df$Label[which(plot_data_df$pct == "")] = ""
  
  return(plot_data_df)
}

plot_pie = function(data_df, color_marker, marker_name){
  
  plot_data_df = process_plot_data(data_df)
  title_marker = paste0(marker_name, " (", nrow(plot_data_df), " genotypes)")
  
  # Define the color palettes
  # colfunc_ama1D3 <- colorRampPalette(rev(c("#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84", "#fc8d59")))
  
  # Plotting one pie with plotly
  a = plot_ly(plot_data_df, labels = ~Haplotype, values = ~Amount, type = 'pie',
              textposition = 'outside',
              text = ~Label, #pct
              textinfo = 'text',
              # textinfo = 'Label',
              marker = list(colors = rep(color_marker, nrow(plot_data_df)), #colfunc_ama1D3(nrow(plot_data_df))
                            line = list(color='black', width=0.5))) %>% ##FFFFFF
    layout(title = title_marker,
           xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
           yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
           uniformtext = list(minsize = 12, mode='hide'), 
           showlegend = FALSE,
           width = 600,
           height = 600)
  # legend = list(title=list(text='<b>Haplotypes
  #                        (N=28):</b>'),
  #               orientation = 'h
  # font = list(size=11))) #add this line for Glurp only (fragment sizes are longer, otherwise legend does not fit)
  return(a)
}

plot_pie_ggplot = function(data_df, color_marker, marker_name) {
  plot_data_df = process_plot_data(data_df)
  title_marker = paste0(marker_name, "\n(", nrow(plot_data_df), " genotypes)")
  
  df2 <- plot_data_df %>% 
    mutate(csum = rev(cumsum(rev(Amount))), 
           pos = Amount/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Amount/2, pos),
           pos_inside = pos - Amount/4,
           empty_label = "")
  
  # Add angle and hjust for better label positioning
  plot_data_df <- plot_data_df %>%
    mutate(mid = cumsum(Amount) - Amount / 2,
           angle = 0, # - 360 * (mid / sum(Amount)),
           hjust = ifelse(angle < -90, 0.1, 0))
  # angle = ifelse(angle < -90, angle + 180, angle))
  
  # plot_data_df <- plot_data_df[order(plot_data_df$Amount), ]
  
  pie_chart = ggplot(plot_data_df, aes(x = "", y = Amount, fill = fct_inorder(as.character(Amount)))) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    guides(fill = guide_legend(title = "Group")) +
    scale_fill_manual(values = rep(color_marker, 100)) +
    scale_y_continuous(breaks = df2$pos, label = df2$empty_label) + #plot_data_df$Label
    geom_text(aes(y = df2$pos, label = Label), size = 7, color = "black", fontface = "bold") +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 17),  
          legend.position = "none", # Removes the legend
          panel.background = element_rect(fill = "white")) +
    ggtitle(title_marker) +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  return(pie_chart)
}

folder_results = "~/GitRepos/bayesian_analysis_Uganda/datasets_no_additional/"
final_table = NULL
list_data = list.files(folder_results, full.names = TRUE)

for (f in list_data) {
  marker_name = str_extract(f, "(?<=MSP ).*?(?=_)")
  data_alleles = read_excel(f, sheet = "RecurrentInfections")
  # Select only day 0:
  data_alleles = data_alleles %>% filter(Day == 0)
  
  data_alleles_long = data_alleles %>% pivot_longer(-c(PatientID, Day, Site), 
                                                    names_to = "allele_name", 
                                                    values_to = "allele_length")
  data_alleles_long$allele_id = substr(data_alleles_long$allele_name, 1, regexpr("_", data_alleles_long$allele_name)-1)
  data_alleles_long$drug = str_extract(data_alleles_long$PatientID, "(?<=_)[A-Z]+(?=_)")
  data_alleles_long = data_alleles_long %>% filter(!is.na(allele_length))
  data_alleles_long$Day = NULL
  data_alleles_long = distinct(data_alleles_long)
  data_alleles_long[which(data_alleles_long$allele_id == "glurp"), "allele_id"] = marker_name
  
  final_table = rbind.data.frame(final_table, data_alleles_long)
  
}

# Final table with unmerged alleles
final_table = distinct(final_table)
list_markers = unique(final_table$allele_id)

# Merge alleles based on bin sizes
merged_allele_table = NULL
for (marker_name in unique(list_markers)) {
  marker_drug_tab = final_table %>% 
    filter(allele_id == marker_name)
  allele_vector = marker_drug_tab %>%
    select(allele_length)
  bin_size = marker_bins %>% 
    filter(marker_id == marker_name) %>% 
    select(bin_size)
  merged_alleles = merge_alleles(allele_vector = allele_vector$allele_length, 
                                 bin_size = bin_size$bin_size)
  
  marker_drug_tab$true_alleles = sapply(marker_drug_tab$allele_length, 
                                        function(x) merged_alleles[which.min(abs(merged_alleles - x))])
  marker_drug_tab$allele_length = marker_drug_tab$allele_name = NULL
  merged_allele_table = rbind.data.frame(merged_allele_table, marker_drug_tab)
}

merged_allele_table = unique(merged_allele_table)
n_patients = length(unique(merged_allele_table$PatientID))

# Figure of diversity across all sites and drug arms
final_table_plot = merged_allele_table %>% 
  group_by(allele_id, true_alleles) %>% 
  summarize(Amount = n())

final_table_plot$Frequency = final_table_plot$Amount/n_patients

colors_amp_seq = c("#fc8d59", "#fc8d59", "#fc8d59", 
                   "#F692EF","#F692EF", "#74a9cf", 
                   "#ef3b2c", "#807dba", "#f768a1", 
                   "grey", "grey", "grey",
                   "#82C0FF", "#41b6c4", "#CD5C5C")
list_markers = unique(final_table_plot$allele_id)
list_markers = c("K1", "MAD20", "R033","3D7", "FC27", "GLURP", 
                 "313", "383", "POLY a", "NULL", "NULL", "NULL", "TA1", "PFPK2", "TA109", "2490")
p_array = vector('list', length(list_markers))

for (i in 1:length(list_markers)) {
  marker_name = list_markers[i]
  print(marker_name)
  final_table_plot1 = final_table_plot %>% filter(allele_id == marker_name) 
  final_table_plot1 = final_table_plot1 %>% ungroup() %>% dplyr::select(true_alleles, Amount, Frequency)
  final_table_plot1 = distinct(final_table_plot1)
  
  if (marker_name != "NULL") {
    p_array[[i]] = plot_pie_ggplot(final_table_plot1, colors_amp_seq[i], marker_name)
  } else {
    p_array[[i]] = NULL
  }
  
}

figure_A = ggarrange(plotlist = p_array, nrow = 3, ncol = 6)

# ggsave("~/genotyping/analysis_Uganda/pie_charts_diversity.pdf",
#        plot = figure_A,  width = 22, height = 10)

# Figure of diversity for each site and drug arm
final_table_plot_drug_site = merged_allele_table %>% 
  group_by(allele_id, drug, Site, true_alleles) %>% 
  summarize(Amount = n())

final_table_plot_pat = unique(merged_allele_table %>% select(PatientID, drug, Site))
final_table_plot_pat = final_table_plot_pat %>% group_by(Site, drug) %>%
                        summarise(patient_count = n()) 
final_table_plot_drug_site = merge(final_table_plot_drug_site, final_table_plot_pat, by = c("Site", "drug"))
final_table_plot_drug_site = final_table_plot_drug_site %>% mutate(Frequency = Amount/patient_count)

colors_amp_seq = c("#fc8d59", "#fc8d59", "#fc8d59", 
                   "#F692EF","#F692EF", "#74a9cf", 
                   "#ef3b2c", "#807dba", "#f768a1", 
                   "grey", "grey", "grey",
                   "#82C0FF", "#41b6c4", "#CD5C5C")
list_markers = unique(final_table_plot_drug_site$allele_id)
list_markers = c("K1", "MAD20", "R033","3D7", "FC27", "GLURP", 
                 "313", "383", "POLY a", "NULL", "NULL", "NULL", "TA1", "PFPK2", "TA109", "2490")

for (drug_name in unique(final_table_plot_drug_site$drug)) {
  for (site in unique(final_table_plot_drug_site$Site)) {
    p_array = NULL
    p_array = vector('list', length(list_markers))
    for (i in 1:length(list_markers)) {
      marker_name = list_markers[i]
      print(marker_name)
      final_table_plot1 = final_table_plot_drug_site %>% filter(allele_id == marker_name &
                                                        drug == drug_name & Site == site) 
      final_table_plot1 = final_table_plot1 %>% ungroup() %>% dplyr::select(true_alleles, Amount, Frequency)
      final_table_plot1 = distinct(final_table_plot1)
      
      if (marker_name != "NULL") {
        p_array[[i]] = plot_pie_ggplot(final_table_plot1, colors_amp_seq[i], marker_name)
      } else {
        p_array[[i]] = NULL
      }
      
    }
    
    figure_A = ggarrange(plotlist = p_array, nrow = 3, ncol = 6)
    figure_A = annotate_figure(figure_A, top = text_grob(paste0("Site ", site, ", ", drug_name, " arm"), 
                                          color = "black", face = "bold", size = 30))
    
    ggsave(paste0("~/genotyping/analysis_Uganda/pie_charts_diversity_", site, "_", drug_name, ".pdf"),
           plot = figure_A,  width = 22, height = 10)
  }
}

# Figure of diversity for each site 
final_table_plot_site = merged_allele_table %>% 
  group_by(allele_id, Site, true_alleles) %>% 
  summarize(Amount = n())

final_table_plot_pat = unique(merged_allele_table %>% select(PatientID, Site))
final_table_plot_pat = final_table_plot_pat %>% group_by(Site) %>%
                        summarise(patient_count = n()) 
final_table_plot_site = merge(final_table_plot_site, final_table_plot_pat, by = c("Site"))
final_table_plot_site = final_table_plot_site %>% mutate(Frequency = Amount/patient_count)

colors_amp_seq = c("#fc8d59", "#fc8d59", "#fc8d59", 
                   "#F692EF","#F692EF", "#74a9cf", 
                   "#ef3b2c", "#807dba", "#f768a1", 
                   "grey", "grey", "grey",
                   "#82C0FF", "#41b6c4", "#CD5C5C")
list_markers = unique(final_table_plot_drug_site$allele_id)
list_markers = c("K1", "MAD20", "R033","3D7", "FC27", "GLURP", 
                 "313", "383", "POLY a", "NULL", "NULL", "NULL", "TA1", "PFPK2", "TA109", "2490")

for (site in unique(final_table_plot_site$Site)) {
  p_array = NULL
  p_array = vector('list', length(list_markers))
  for (i in 1:length(list_markers)) {
    marker_name = list_markers[i]
    print(marker_name)
    final_table_plot1 = final_table_plot_site %>% filter(allele_id == marker_name & Site == site) 
    final_table_plot1 = final_table_plot1 %>% ungroup() %>% dplyr::select(true_alleles, Amount, Frequency)
    final_table_plot1 = distinct(final_table_plot1)
    
    if (marker_name != "NULL") {
      p_array[[i]] = plot_pie_ggplot(final_table_plot1, colors_amp_seq[i], marker_name)
    } else {
      p_array[[i]] = NULL
    }
    
  }
  
  figure_A = ggarrange(plotlist = p_array, nrow = 3, ncol = 6)
  figure_A = annotate_figure(figure_A, top = text_grob(paste0("Site ", site), 
                                                       color = "black", face = "bold", size = 30))
  
  ggsave(paste0("~/genotyping/analysis_Uganda/pie_charts_diversity_", site, ".pdf"),
         plot = figure_A,  width = 22, height = 10)
}



