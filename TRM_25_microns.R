#TRM_25_microns.R#
#Travis Tu; Last Edited: 5/8/2017#
#Folders Needed: Colt_TRM, Recoordinate_Dist_Files, Filtered#

#Install and load required packages
required = c("data.table")
for (lib in required)
{
  if (!require(lib, character.only=TRUE))
  {
    install.packages(lib, repos="http://cran.rstudio.com")
    suppressMessages(library(lib, character.only=TRUE, quietly=TRUE))
  }
}

library(data.table)

setwd("C:/Users/ttu/Documents/TRM/Files")

RDF_all_files <- list.files("C:/Users/ttu/Documents/TRM/Files")

final <- as.data.frame(matrix(0, ncol = 5, nrow = length(RDF_all_files)))
colnames(final) <- c("Patient_ID", "CK_to_CD8+CD103+", "CK_to_CD8+CD103-", "CD8+CD103+_to_CK", "CD8+CD103-_to_CK")
final[, 1] <- RDF_all_files

final_Cancer <- as.data.frame(matrix(0, ncol = 5, nrow = length(RDF_all_files)))
colnames(final_Cancer) <- c("Patient_ID", "CK_to_CD8+CD103+", "CK_to_CD8+CD103-", "CD8+CD103+_to_CK", "CD8+CD103-_to_CK")
final_Cancer[, 1] <- RDF_all_files

final_Stroma <- as.data.frame(matrix(0, ncol = 5, nrow = length(RDF_all_files)))
colnames(final_Stroma) <- c("Patient_ID", "CK_to_CD8+CD103+", "CK_to_CD8+CD103-", "CD8+CD103+_to_CK", "CD8+CD103-_to_CK")
final_Stroma[, 1] <- RDF_all_files

final_Tissue <- as.data.frame(matrix(0, ncol = 5, nrow = length(RDF_all_files)))
colnames(final_Tissue) <- c("Patient_ID", "CK_to_CD8+CD103+", "CK_to_CD8+CD103-", "CD8+CD103+_to_CK", "CD8+CD103-_to_CK")
final_Tissue[, 1] <- RDF_all_files

TA <- read.delim("C:/Users/ttu/Documents/TRM/Tissue_Areas/Colt_TRM_Sum_Areas.txt", header = TRUE)

for(i in 1:length(RDF_all_files)) {
  file <- fread(RDF_all_files[i], header = TRUE)
  
  row_ID <- gsub(".*\\Recoor_", "\\1", RDF_all_files[i])
  row_ID <- gsub("_cell+.*", "\\1", row_ID)
  row_num <- which(grepl(row_ID, TA$Patient_ID))
  
  CK_CD8_CD103 <- file[file$Phenotype == "CK" & file$"Distance to CD8+CD103+" <= 50, ] #CK to CD8+CD103+
  CK_CD8 <- file[file$Phenotype == "CK" & file$"Distance to CD8+CD103-" <= 50, ] #CK to CD8+CD103-
  
  CD8_CD103_CK <- file[file$Phenotype == "CD8+CD103+" & file$"Distance to CK" <= 50, ] #50 pixels = 25 microns
  CD8_CK <- file[file$Phenotype == "CD8+CD103-" & file$"Distance to CK" <= 50, ]
  
  CK_CD8_CD103_Cancer <- nrow(CK_CD8_CD103)/((TA[row_num, 2] * 0.246)/1000000)
  CK_CD8_Cancer <- nrow(CK_CD8)/((TA[row_num, 2] * 0.246)/1000000)
  CD8_CD103_Cancer <- nrow(CD8_CD103_CK)/((TA[row_num, 2] * 0.246)/1000000)
  CD8_CK_Cancer <- nrow(CD8_CK)/((TA[row_num, 2] * 0.246)/1000000)
  
  CK_CD8_CD103_Stroma <- nrow(CK_CD8_CD103)/((TA[row_num, 3] * 0.246)/1000000)
  CK_CD8_Stroma <- nrow(CK_CD8)/((TA[row_num, 3] * 0.246)/1000000)
  CD8_CD103_Stroma <- nrow(CD8_CD103_CK)/((TA[row_num, 3] * 0.246)/1000000)
  CD8_CK_Stroma <- nrow(CD8_CK)/((TA[row_num, 3] * 0.246)/1000000)
  
  CK_CD8_CD103_Tissue <- nrow(CK_CD8_CD103)/((TA[row_num, 4] * 0.246)/1000000)
  CK_CD8_Tissue <- nrow(CK_CD8)/((TA[row_num, 4] * 0.246)/1000000)
  CD8_CD103_Tissue <- nrow(CD8_CD103_CK)/((TA[row_num, 4] * 0.246)/1000000)
  CD8_CK_Tissue <- nrow(CD8_CK)/((TA[row_num, 4] * 0.246)/1000000)
  
  final[i, c(2, 3, 4, 5)] <- c(nrow(CK_CD8_CD103), nrow(CK_CD8), nrow(CD8_CD103_CK), nrow(CD8_CK))
  final_Cancer[i, c(2, 3, 4, 5)] <- c(CK_CD8_CD103_Cancer, CK_CD8_Cancer, CD8_CD103_Cancer, CD8_CK_Cancer)
  final_Stroma[i, c(2, 3, 4, 5)] <- c(CK_CD8_CD103_Stroma, CK_CD8_Stroma, CD8_CD103_Stroma, CD8_CK_Stroma)
  final_Tissue[i, c(2, 3, 4, 5)] <- c(CK_CD8_CD103_Tissue, CK_CD8_Tissue, CD8_CD103_Tissue, CD8_CK_Tissue)
}
fwrite(final, file = "C:/Users/ttu/Documents/TRM/Filtered/TRM_25.txt", sep = "\t", row.names = FALSE)
fwrite(final_Cancer, file = "C:/Users/ttu/Documents/TRM/Filtered/TRM_25_Cancer_Normalized.txt", sep = "\t", row.names = FALSE)
fwrite(final_Stroma, file = "C:/Users/ttu/Documents/TRM/Filtered/TRM_25_Stroma_Normalized.txt", sep = "\t", row.names = FALSE)
fwrite(final_Tissue, file = "C:/Users/ttu/Documents/TRM/Filtered/TRM_25_Tissue_Normalized.txt", sep = "\t", row.names = FALSE)
