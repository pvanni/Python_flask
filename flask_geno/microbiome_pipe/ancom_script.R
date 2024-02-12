#install.packages('exactRankTests')
#install.packages('nlme')
#install.packages('dplyr')
#install.packages('ggplot2')
#install.packages('compositions')
#install.packages('readr')
#install.packages('tidyverse')
#install.packages('stringr')

library(exactRankTests)
library(nlme)
library(dplyr)
library(ggplot2)
library(compositions)
library(readr)
library(tidyverse)
library(stringr)

target_variables <- c('Koliikki', 'HA_kombi', 'Atopia_kombi', 'Maito_kombi', 'obesity_vs_control', 'obesity_and_ow_vs_control')
data_types <- c('16s_meco', '16s_6mo', '16s_18mo', 'its_meco', 'its_6mo', 'its_18mo')
lib_cuts <- c(1000, 1000, 1000, 500, 500, 500)

mapply(function(data_type, lib_cut) {
  
  for (target in target_variables){
  
    group = "meco"
    main_var = target
    
    setwd("D:/väikkäri/kolmas/2023/ancom_analyses")
    table_path = str_glue('{data_type}_feature_table.tsv')
    meta_path = str_glue('{data_type}_{main_var}_metadata.tsv')
    output_name = str_glue('{data_type}_{main_var}_results.tsv')
    
    # This code is from https://github.com/FrederickHuangLin/ANCOM for 2 group differential abundance analysis
    otu_data = read_tsv(table_path, skip=1)
    #otu_data = read_tsv(table_path)
    otu_id = otu_data$"#OTU ID"
    otu_data = data.frame(otu_data[,-1], check.names=FALSE)
    rownames(otu_data) = otu_id
    
    meta_data = read_tsv(meta_path)
    
    source("D:/väikkäri/kolmas/scripts/r_files/ancom_v2.1.R")
    
    prepro = feature_table_pre_process(feature_table = otu_data, meta_data = meta_data, 
                                       sample_var = "#SampleID", group_var = NULL, 
                                       out_cut = 0.05, zero_cut = 0.90, lib_cut = lib_cut, neg_lb = FALSE)
    
    feature_table = prepro$feature_table
    meta_data = prepro$meta_data
    struc_zero = prepro$structure_zeros
    
    res = ANCOM(feature_table = feature_table, meta_data = meta_data, 
                struc_zero = struc_zero, main_var = main_var, 
                p_adj_method = "BH", alpha = 0.05, adj_formula = NULL, rand_formula = NULL)
    
    write_tsv(res$out, output_name)
  }
}, data_types, lib_cuts)