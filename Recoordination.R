#Recoordination.R#
#Travis Tu; Last Edited: 8/22/2018#
#Folders needed: File Merging, Files to merge, Recor

setwd("C:/Users/ttu/Documents/File Merging/Files to merge")
all_files <- list.files("C:/Users/ttu/Documents/File Merging/Files to merge")

for (file in all_files){
  #Imports data from .txt file#
  merged <- read.delim(file, sep = "\t", header = TRUE)
  #Extract stage coordinates from Sample_Name#
  Stage_X <- gsub(".*_\\[+([0-9])", "\\1", merged$Sample.Name)
  Stage_X <- gsub("\\,+.*", "\\1", Stage_X)
  Stage_X <- as.numeric(Stage_X)
  merged$Stage_X_Position <- Stage_X
  
  Stage_Y <- gsub(".*\\,+([0-9])", "\\1", merged$Sample.Name)
  Stage_Y <- gsub("]+.*", "\\1", Stage_Y)
  Stage_Y <- as.numeric(Stage_Y)
  merged$Stage_Y_Position <- Stage_Y
  
  #Get dimensions of the image
  #x_dim <- 1392
  x_dim <- 1352
  #x_dim <- 1340
  #x_dim <- 1344
  
  #y_dim <- 1040
  y_dim <- 1012
  #y_dim <- 1004
  #y_dim <- 1008
  
  #Get difference between centers
  #x_center_diff <- 695.2
  x_center_diff <- 669
  #x_center_diff <- 669
  #x_center_diff <- 669
  
  #y_center_diff <- 519.4
  y_center_diff <- 500
  #y_center_diff <- 500
  #y_center_diff <- 500
  
  #Obtain minimum X and minimum Y (top left corner of image)#
  merged$Min_X <- min(merged$Stage_X_Position)
  merged$Min_Y <- min(merged$Stage_Y_Position)
  
  #Obtain offset for X and Y
  merged$X_Offset <- merged$Min_X
  merged$Y_Offset <- merged$Min_Y
  
  merged$X_Diff <- merged$Stage_X_Position - merged$Min_X
  merged$Cell.X.Position <- merged$X_Offset + merged$Cell.X.Position + x_dim*(merged$X_Diff/x_center_diff)
  
  merged$Y_Diff <- merged$Stage_Y_Position - merged$Min_Y
  merged$Cell.Y.Position <- merged$Y_Offset + merged$Cell.Y.Position + y_dim*(merged$Y_Diff/y_center_diff)
  
  #Output new file with corrected cell coordinates#
  write.table(merged, file = paste("C:/Users/ttu/Documents/File Merging/Recor/Recor_", file, sep = ""), sep = "\t", row.names = FALSE)
}
