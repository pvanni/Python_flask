library(exactRankTests)
library(nlme)
library(dplyr)
library(ggplot2)
library(compositions)
library(readr)
library(tidyverse)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

main_var = args[1]
data_type = args[2]
project_path = args[3]
lib_cut = 1000

# Fix this
setwd("D:/väikkäri/kolmas/2023/ancom_analyses")
# Fix this
table_path = str_glue('{data_type}_feature_table.tsv')
# {project_path}/start_files/metadata.tsv
meta_path = str_glue('{project_path}/start_files/metadata.tsv')
# Fix this
output_name = str_glue('{data_type}_{main_var}_results.tsv')

# This code is from https://github.com/FrederickHuangLin/ANCOM for 2 group differential abundance analysis
otu_data = read_tsv(table_path, skip=1)
#otu_data = read_tsv(table_path)
otu_id = otu_data$"#OTU ID"
otu_data = data.frame(otu_data[,-1], check.names=FALSE)
rownames(otu_data) = otu_id

meta_data = read_tsv(meta_path)

# Fix this source to the flask_geno folder
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