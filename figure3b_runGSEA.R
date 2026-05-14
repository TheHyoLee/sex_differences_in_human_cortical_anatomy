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
basePath <- "/path_to_results/"
inPath <- "/path_to_results/input_main/"
outPath <- "/path_to_results/output_main/"

# Read the ranked gene list
magicc <- read_excel(paste0(basePath, "CorticalSexDiff_MAGICC_Correlation.xlsx"))

# Extract pathways and ranks (magicc)
magicc_pathways <- lapply(magicc[,8:ncol(magicc)], function(boolean_column) {magicc$Gene_ID[boolean_column]})
names(magicc_pathways) <- colnames(magicc)[8:ncol(magicc)]

magicc_rank_fcv <- magicc$Female_CV
names(magicc_rank_fcv) <- magicc$Gene_ID
magicc_rank_fcv <- sort(magicc_rank_fcv, decreasing = TRUE)

magicc_rank_fsa <- magicc$Female_SA
names(magicc_rank_fsa) <- magicc$Gene_ID
magicc_rank_fsa <- sort(magicc_rank_fsa, decreasing = TRUE)

magicc_rank_fct <- magicc$Female_CT
names(magicc_rank_fct) <- magicc$Gene_ID
magicc_rank_fct <- sort(magicc_rank_fct, decreasing = TRUE)

magicc_rank_mcv <- magicc$Male_CV
names(magicc_rank_mcv) <- magicc$Gene_ID
magicc_rank_mcv <- sort(magicc_rank_mcv, decreasing = TRUE)

magicc_rank_msa <- magicc$Male_SA
names(magicc_rank_msa) <- magicc$Gene_ID
magicc_rank_msa <- sort(magicc_rank_msa, decreasing = TRUE)

magicc_rank_mct <- magicc$Male_CT
names(magicc_rank_mct) <- magicc$Gene_ID
magicc_rank_mct <- sort(magicc_rank_mct, decreasing = TRUE)

# Run fGSEA (magicc)
magicc_fsea_fcv <- fgsea(pathways = magicc_pathways,
                         stats = magicc_rank_fcv,
                         scoreType = 'pos',
                         minSize = 1,
                         maxSize = length(magicc_rank_fcv)-1,
                         nPermSimple = 10000)
write_xlsx(magicc_fsea_fcv, paste0(outPath, "MAGICC_GSEA_F_CV.xlsx"))

magicc_fsea_fsa <- fgsea(pathways = magicc_pathways,
                         stats = magicc_rank_fsa,
                         scoreType = 'pos',
                         minSize = 1,
                         maxSize = length(magicc_rank_fsa)-1,
                         nPermSimple = 10000)
write_xlsx(magicc_fsea_fsa, paste0(outPath, "MAGICC_GSEA_F_SA.xlsx"))

magicc_fsea_fct <- fgsea(pathways = magicc_pathways,
                         stats = magicc_rank_fct,
                         scoreType = 'pos',
                         minSize = 1,
                         maxSize = length(magicc_rank_fct)-1,
                         nPermSimple = 10000)
write_xlsx(magicc_fsea_fct, paste0(outPath, "MAGICC_GSEA_F_CT.xlsx"))

magicc_fsea_mcv <- fgsea(pathways = magicc_pathways,
                         stats = magicc_rank_mcv,
                         scoreType = 'pos',
                         minSize = 1,
                         maxSize = length(magicc_rank_mcv)-1,
                         nPermSimple = 10000)
write_xlsx(magicc_fsea_mcv, paste0(outPath, "MAGICC_GSEA_M_CV.xlsx"))

magicc_fsea_msa <- fgsea(pathways = magicc_pathways,
                         stats = magicc_rank_msa,
                         scoreType = 'pos',
                         minSize = 1,
                         maxSize = length(magicc_rank_msa)-1,
                         nPermSimple = 10000)
write_xlsx(magicc_fsea_msa, paste0(outPath, "MAGICC_GSEA_M_SA.xlsx"))

magicc_fsea_mct <- fgsea(pathways = magicc_pathways,
                         stats = magicc_rank_mct,
                         scoreType = 'pos',
                         minSize = 1,
                         maxSize = length(magicc_rank_mct)-1,
                         nPermSimple = 10000)
write_xlsx(magicc_fsea_mct, paste0(outPath, "MAGICC_GSEA_M_CT.xlsx"))