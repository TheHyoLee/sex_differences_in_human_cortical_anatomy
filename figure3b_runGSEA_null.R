# This script is meant to be run on computing cluster to run GSEA for 10,000 spin-permutation-based ranked gene lists

# Parse input argument
args <- commandArgs(trailingOnly=TRUE)
sBS = as.character(args[1])

# Import libraries
library(tidyverse)
library(RColorBrewer)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(readxl)
library(writexl)

# Set paths
basePath <- "/path_to_results/GSEA_R_UKB/"
inPath <- "/path_to_results/input_bs/"
outPath <- "/path_to_results/output_bs/"

# Read the ranked gene list
magicc <- read_excel(paste0(basePath, "CorticalSexDiff_MAGICC_Correlation.xlsx"))

# Extract pathways and ranks (magicc)
magicc_pathways <- lapply(magicc[,8:ncol(magicc)], function(boolean_column) {magicc$Gene_ID[boolean_column]})
names(magicc_pathways) <- colnames(magicc)[8:ncol(magicc)]

magicc_rank_fcv <- scan(paste0(inPath, "MAGICC_R_CV_F_Spin_BS", sBS, ".txt"))
names(magicc_rank_fcv) <- magicc$Gene_ID
magicc_rank_fcv <- sort(magicc_rank_fcv, decreasing = TRUE)

magicc_rank_fsa <- scan(paste0(inPath, "MAGICC_R_SA_F_Spin_BS", sBS, ".txt"))
names(magicc_rank_fsa) <- magicc$Gene_ID
magicc_rank_fsa <- sort(magicc_rank_fsa, decreasing = TRUE)

magicc_rank_fct <- scan(paste0(inPath, "MAGICC_R_CT_F_Spin_BS", sBS, ".txt"))
names(magicc_rank_fct) <- magicc$Gene_ID
magicc_rank_fct <- sort(magicc_rank_fct, decreasing = TRUE)

magicc_rank_mcv <- scan(paste0(inPath, "MAGICC_R_CV_M_Spin_BS", sBS, ".txt"))
names(magicc_rank_mcv) <- magicc$Gene_ID
magicc_rank_mcv <- sort(magicc_rank_mcv, decreasing = TRUE)

magicc_rank_msa <- scan(paste0(inPath, "MAGICC_R_SA_M_Spin_BS", sBS, ".txt"))
names(magicc_rank_msa) <- magicc$Gene_ID
magicc_rank_msa <- sort(magicc_rank_msa, decreasing = TRUE)

magicc_rank_mct <- scan(paste0(inPath, "MAGICC_R_CT_M_Spin_BS", sBS, ".txt"))
names(magicc_rank_mct) <- magicc$Gene_ID
magicc_rank_mct <- sort(magicc_rank_mct, decreasing = TRUE)


# Run fGSEA (magicc)
magicc_fsea_fcv <- fgsea(pathways = magicc_pathways,
                         stats = magicc_rank_fcv,
                         scoreType = 'pos',
                         minSize = 1,
                         maxSize = length(magicc_rank_fcv)-1,
                         nPermSimple = 10000)
write_xlsx(magicc_fsea_fcv, paste0(outPath, "MAGICC_GSEA_F_CV_BS", sBS, ".xlsx"))

magicc_fsea_fsa <- fgsea(pathways = magicc_pathways,
                         stats = magicc_rank_fsa,
                         scoreType = 'pos',
                         minSize = 1,
                         maxSize = length(magicc_rank_fsa)-1,
                         nPermSimple = 10000)
write_xlsx(magicc_fsea_fcv, paste0(outPath, "MAGICC_GSEA_F_SA_BS", sBS, ".xlsx"))

magicc_fsea_fct <- fgsea(pathways = magicc_pathways,
                         stats = magicc_rank_fct,
                         scoreType = 'pos',
                         minSize = 1,
                         maxSize = length(magicc_rank_fct)-1,
                         nPermSimple = 10000)
write_xlsx(magicc_fsea_fcv, paste0(outPath, "MAGICC_GSEA_F_CT_BS", sBS, ".xlsx"))

magicc_fsea_mcv <- fgsea(pathways = magicc_pathways,
                         stats = magicc_rank_mcv,
                         scoreType = 'pos',
                         minSize = 1,
                         maxSize = length(magicc_rank_mcv)-1,
                         nPermSimple = 10000)
write_xlsx(magicc_fsea_fcv, paste0(outPath, "MAGICC_GSEA_M_CV_BS", sBS, ".xlsx"))

magicc_fsea_msa <- fgsea(pathways = magicc_pathways,
                         stats = magicc_rank_msa,
                         scoreType = 'pos',
                         minSize = 1,
                         maxSize = length(magicc_rank_msa)-1,
                         nPermSimple = 10000)
write_xlsx(magicc_fsea_fcv, paste0(outPath, "MAGICC_GSEA_M_SA_BS", sBS, ".xlsx"))


magicc_fsea_mct <- fgsea(pathways = magicc_pathways,
                         stats = magicc_rank_mct,
                         scoreType = 'pos',
                         minSize = 1,
                         maxSize = length(magicc_rank_mct)-1,
                         nPermSimple = 10000)
write_xlsx(magicc_fsea_fcv, paste0(outPath, "MAGICC_GSEA_M_CT_BS", sBS, ".xlsx"))
