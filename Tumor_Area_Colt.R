#Tumor_Area_Colt.R#
#Travis Tu; Last Edited: 11/15/2016#
#Folders: Tumor_Area, Sum

setwd("W:/Desktop/DATA_Transformation_Scripts/Tissue Seg")

all_files_new <- list.files("W:/Desktop/DATA_Transformation_Scripts/Tissue Seg")

final <- as.data.frame(matrix(0, ncol = 7, nrow = length(all_files_new)))
colnames(final) <- c("Patient_ID", "Tumor_Area_Sum", "Stroma_Area_Sum", "Tissue_Area_Sum", 
                     "Tumor_Area_mm2", "Stroma_Area_mm2", "Tissue_Area_mm2")
final[, 1] <- all_files_new

for (i in 1:length(all_files_new)) {
  file_to_sum <- read.delim(all_files_new[i], sep = "\t", header = TRUE)
  
  TC_Index <- which(grepl("Tissue.Category", names(file_to_sum)))
  colnames(file_to_sum)[TC_Index] <- "Tissue_Category"
  Region_Index <- which(grepl("Region.Area..pixels.", names(file_to_sum)))
  colnames(file_to_sum)[Region_Index] <- "Region_Area_Pixels"
  
  Tumor <- file_to_sum[toupper(file_to_sum$Tissue_Category) == "TUMOR", ] #Change your tumor tissue category#
  Stroma <- file_to_sum[toupper(file_to_sum$Tissue_Category) == "STROMA", ] #Change your stroma tissue category#
  
  Tumor_Sum <- sum(Tumor$Region_Area_Pixels)
  Stroma_Sum <- sum(Stroma$Region_Area_Pixels)
  Tissue_Sum <- sum(Tumor_Sum, Stroma_Sum)
  
  Tumor_Sum_mm2 <- (Tumor_Sum * 0.246)/1000000
  Stroma_Sum_mm2 <- (Stroma_Sum * 0.246)/1000000
  Tissue_Sum_mm2 <- sum(Tumor_Sum_mm2, Stroma_Sum_mm2)
  
  final[i, c(2:7)] <- c(Tumor_Sum, Stroma_Sum, Tissue_Sum, Tumor_Sum_mm2, Stroma_Sum_mm2, Tissue_Sum_mm2)
}

tumor_fail <- which(final[, 2] == 0)
stroma_fail <- which(final[, 3] == 0)

if (length(tumor_fail) > 0) {
  message(paste("Tumor Area was 0 for", toString(final[tumor_fail, 1])))
  for (i in 1:length(final[tumor_fail, 1])) {
    tumor2 <- read.delim(final[tumor_fail[i], 1], header = TRUE)
    unique_tc <- as.vector(unique(tumor2$Tissue.Category))
    if (i == 1) {
      tc <- unique_tc
    }
    else {
      tc <- c(tc, unique_tc)
    }
  }
  final_tumor_tc <- toupper(unique(tc))
  message(paste("Valid choices for these files are:", toString(final_tumor_tc)))
}

if (length(stroma_fail) > 0) {
  message(paste("Stroma Area was 0 for", toString(final[stroma_fail, 1])))
  for (i in 1:length(final[stroma_fail, 1])) {
    stroma2 <- read.delim(final[stroma_fail[i], 1], header = TRUE)
    unique_tc <- as.vector(unique(stroma2$Tissue.Category))
    if (i == 1) {
      tc <- unique_tc
    }
    else {
      tc <- c(tc, unique_tc)
    }
  }
  final_stroma_tc <- toupper(unique(tc))
  message(paste("Valid choices for these files are:", toString(final_stroma_tc)))
}

write.table(final, file = "W:/Desktop/DATA_Transformation_Scripts/Sum/Sum_Areas.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE)
