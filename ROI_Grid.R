#ROI_Grid.R#
#Travis Tu; Last Edited: 5/31/2017#
#Folders Needed: ROI_Grid, XML, Dataframe#

#Install and load required packages
required = c("RANN", "data.table", "plyr")
for (lib in required)
{
  if (!require(lib, character.only=TRUE))
  {
    install.packages(lib, repos="http://cran.rstudio.com")
    suppressMessages(library(lib, character.only=TRUE, quietly=TRUE))
  }
}

library(data.table)
library(RANN)
library(plyr)

#Binary Search#
search <- function(array, targetValue) {
  min <- 1
  max <- length(array)
  while (min <= max) {
    guess <- floor((min + max)/2)
    if (array[guess] == targetValue) {
      return (guess)
    }
    else if (array[guess] < targetValue) {
      min <- guess + 1
    }
    else {
      max <- guess - 1
    }
  }
}

#Go through XML file to pull out coordinates of high powered fields
#These were selected by the ROI circling tool and may not match with the .im3 files
setwd("C:/Users/ttu/Documents/ROI_Grid/XML")

list_XML <- list.files("C:/Users/ttu/Documents/ROI_Grid/XML")

for (i in 1:length(list_XML)){
  doc <- read.delim(list_XML[i])
  X_Index <- which(grepl("<X>", doc[,1]))
  Y_Index <- which(grepl("<Y>", doc[,1]))
  
  X_list <- c()
  Y_list <- c()
  
  for (j in 1:length(X_Index)) {
    X_Coor <- doc[X_Index[j], 1]
    X_Coor <- gsub(".*<X>", "\\1", X_Coor)
    X_Coor <- gsub("</X>.*", "\\1", X_Coor)
    
    X_Coor <- as.numeric(X_Coor)
    X_Coor <- X_Coor + 335 #Add in X offset here
    X_Coor <- floor(X_Coor)
    X_list[j] <- X_Coor
  }
  
  for (k in 1:length(Y_Index)) {
    Y_Coor <- doc[Y_Index[k], 1]
    Y_Coor <- gsub(".*<Y>", "\\1", Y_Coor)
    Y_Coor <- gsub("</Y>.*", "\\1", Y_Coor)
    
    Y_Coor <- as.numeric(Y_Coor)
    Y_Coor <- Y_Coor + 250 #Add in Y offset here
    Y_Coor <- floor(Y_Coor)
    Y_list[k] <- Y_Coor
  }
  
  X_DF <- as.data.frame(X_list)
  Y_DF <- as.data.frame(Y_list)
  XY_DF <- cbind(X_DF, Y_DF)
  XY_DF$ID <- 1:nrow(XY_DF)
  write.table(XY_DF, file = paste("C:/Users/ttu/Documents/ROI_Grid/Dataframe/", list_XML[i], sep = ""), 
              row.names = FALSE, sep = "\t")
}

#Extract the stage coordinates of the .im3 files from a list of .im3 files
setwd("C:/Users/ttu/Documents/ROI_Grid/Images")

list_Images <- list.files("C:/Users/ttu/Documents/ROI_Grid/Images")

if (length(list_Images) > 0) {
  X_Image_list <- c()
  Y_Image_list <- c()
  
  for (i in 1:length(list_Images)) {
    Image_X_Coor <- gsub(".*_\\[", "\\1", list_Images[i])
    Image_X_Coor <- gsub(",.*", "\\1", Image_X_Coor)
    Image_X_Coor <- as.numeric(Image_X_Coor)
    X_Image_list[i] <- Image_X_Coor
    
    Image_Y_Coor <- gsub(".*,", "\\1", list_Images[i])
    Image_Y_Coor <- gsub("\\].*", "\\1", Image_Y_Coor)
    Image_Y_Coor <- as.numeric(Image_Y_Coor)
    Y_Image_list[i] <- Image_Y_Coor
    
    Image_X_DF <- as.data.frame(X_Image_list)
    Image_Y_DF <- as.data.frame(Y_Image_list)
    Image_XY_DF <<- cbind(Image_X_DF, Image_Y_DF)
    Image_XY_DF$Name <- paste(Image_XY_DF[,1], ",", Image_XY_DF[,2], sep = "")
    Image_XY_DF$ID <- 1:nrow(Image_XY_DF)
  }
  
  #Use nearest neighbors to find closest .im3 file to the coordinates of the boxes created by ROI circling
  setwd("C:/Users/ttu/Documents/ROI_Grid/Dataframe")
  
  list_DF <- list.files("C:/Users/ttu/Documents/ROI_Grid/Dataframe")
  
  for (k in 1:length(list_DF)) {
    query_DF <- read.delim(list_DF[k], header = TRUE)
    query <- query_DF[,1:2]
    data <- Image_XY_DF[,1:2]
    
    nearest_neighbor <<- nn2(data, query = query, k=nrow(data), treetype = "kd", searchtype = "standard")
    write.table(nearest_neighbor, file = paste("C:/Users/ttu/Documents/ROI_Grid/Nearest Neighbors/NN_Images_", 
                                               list_DF[k], sep = ""), sep = "\t", row.names = FALSE)
  }
}

#Extract stage coordinates for the .im3 files from a merged cell_seg_data.txt file
setwd("C:/Users/ttu/Documents/ROI_Grid/Merged")

list_Merged <- list.files("C:/Users/ttu/Documents/ROI_Grid/Merged")

if (length(list_Merged) > 0) {
  
  merged_X_list <- c()
  merged_Y_list <- c()
  
  for (i in 1:length(list_Merged)) {
    merged <<- fread(list_Merged[i], header = TRUE, stringsAsFactors = FALSE)
    
    SN_Unique <- unique(merged$"Sample Name")
    
    #Extract stage coordinates from Sample_Name#
    for (j in 1:length(SN_Unique)) {
      merged_X <- gsub(".*_\\[+([0-9])", "\\1", SN_Unique)
      merged_X <- gsub("\\,+.*", "\\1", merged_X)
      merged_X <- as.numeric(merged_X)
      merged_X_list <- merged_X
      
      merged_Y <- gsub(".*\\,+([0-9])", "\\1", SN_Unique)
      merged_Y <- gsub("]+.*", "\\1", merged_Y)
      merged_Y <- as.numeric(merged_Y)
      merged_Y_list <- merged_Y
    }
    
    merged_X_DF <- as.data.frame(merged_X_list)
    merged_Y_DF <- as.data.frame(merged_Y_list)
    merged_XY_DF <<- cbind(merged_X_DF, merged_Y_DF)
    merged_XY_DF$Name <- paste(merged_XY_DF[,1], ",", merged_XY_DF[,2], sep = "")
    merged_XY_DF$ID <- 1:nrow(merged_XY_DF)
  }
  
  #Use nearest neighbors to find closest .im3 file to the coordinates of the boxes created by ROI circling
  setwd("C:/Users/ttu/Documents/ROI_Grid/Dataframe")
  
  list_DF <- list.files("C:/Users/ttu/Documents/ROI_Grid/Dataframe")
  
  for (k in 1:length(list_DF)) {
    query_DF <- read.delim(list_DF[k], header = TRUE)
    query <- query_DF[,1:2]
    data <- merged_XY_DF[,1:2]
    
    nearest_neighbor <<- nn2(data, query = query, k=nrow(data), treetype = "kd", searchtype = "standard")
    write.table(nearest_neighbor, file = paste("C:/Users/ttu/Documents/ROI_Grid/Nearest Neighbors/NN_Merged_", 
                                               list_DF[k], sep = ""), sep = "\t", row.names = FALSE)
  }
}

#Copy the images closest to the ROI circling created grid boxes
setwd("C:/Users/ttu/Documents/ROI_Grid/Nearest Neighbors")

list_NN <- list.files("C:/Users/ttu/Documents/ROI_Grid/Nearest Neighbors")

ROI_Image_list <<- c()
ROI_Image_list_2 <<- c()
ROI_Merged_list <<- list()

for (i in 1:length(list_NN)) {
  NN <- read.delim(list_NN[i], header = TRUE, stringsAsFactors = FALSE)
  
  if (grepl("Images", list_NN[i]) == TRUE) {
    for (j in 1:length(NN[,1])) {
      Image_ID <- search(Image_XY_DF$ID, NN[j, 1])
      
      Name_Index <- which(grepl("Name", names(Image_XY_DF)))
      IM3_search <- Image_XY_DF[Image_ID, Name_Index]
      ROI_Image_Index <- which(grepl(IM3_search, list_Images))
      
      ROI_Image_list[j] <- list_Images[ROI_Image_Index]
    }
    
    for (k in 1:length(ROI_Image_list)) {
      ROI_Image_list_2[k] <- paste("C:/Users/ttu/Documents/ROI_Grid/Images/", ROI_Image_list[k], sep = "")
    }
    
    file.copy(from = ROI_Image_list_2, "C:/Users/ttu/Documents/ROI_Grid/ROI_Images/" )
  }
  
  else if (grepl("Merged", list_NN[i]) == TRUE) {
    for (j in 1:length(NN[,1])) {
      merged_ID <- search(merged_XY_DF$ID, NN[j, 1])
      
      m_Name_Index <- which(grepl("Name", names(merged_XY_DF)))
      merged_search <- merged_XY_DF[merged_ID, m_Name_Index]
      ROI_merged_Index <- which(grepl(merged_search, merged$"Sample Name"))
      
      ROI_Merged_list[[j]] <- ROI_merged_Index
    }
    for (k in 1:length(ROI_Merged_list)) {
      ROI_merged <- merged[ROI_Merged_list[[k]], ]
      write.table(ROI_merged, file = paste("C:/Users/ttu/Documents/ROI_Grid/ROI_Merged/",k,"_ROI_Merged.txt", sep =""), 
                  sep = "\t", row.names = FALSE)
    }
    
    #setwd("C:/Users/ttu/Documents/ROI_Grid/ROI_Merged")
    #all_files <- list.files("C:/Users/ttu/Documents/ROI_Grid/ROI_Merged")
    
    #dataset <- do.call("rbind",lapply(all_files, FUN=function(files){ 
      #fread(files, header=TRUE, sep="\t")}))
    
    #fwrite(dataset, file = "ROI_Merged.txt", sep = "\t", row.names = FALSE)
  } 
}

#Deduplicate.R#
#Travis Tu; Last Edited: 6/16/2017#
#Folders needed: Folders from ROI_Grid, Deduplicated#

#De-Duplicate ROI_Merged

special <- function(dataframe) {
  Sample_Name_Index <- which(grepl("Sample.Name", names(dataframe)))
  Stage <- gsub("\\.im3.*", "\\1", dataframe[,Sample_Name_Index])
  Stage1 <- gsub("\\[.*", "\\1", Stage)
  Stage2 <- gsub(".*\\[", "\\1", Stage)
  Stage3 <- gsub("\\].*", "\\1", Stage2)
  string <- paste(Stage1, "\\[", Stage3, "\\]", sep = "")
  dataframe$Sample_ID_Extracted2 <- string
  return(dataframe)
}

special_2 <- function(dataframe) {
  Sample_Name_Index <- which(grepl("Sample_ID_Extracted2", names(dataframe)))
  Stage <- gsub("\\.im3.*", "\\1", dataframe[,Sample_Name_Index])
  Stage1 <- gsub("\\(.*", "\\1", Stage)
  Stage2 <- gsub(".*\\(", "\\1", Stage)
  Stage3 <- gsub("\\).*", "\\1", Stage2)
  Stage2 <- gsub(".*\\)", "\\1", Stage2)
  string <- paste(Stage1, "\\(", Stage3, "\\)", Stage2, sep = "")
  dataframe$Sample_ID_Extracted2 <- string
  return(dataframe)
}

#Remove any duplicated files from the merged_cell_seg_data.txt portion
setwd("C:/Users/ttu/Documents/ROI_Grid/ROI_Merged")
all_files <- list.files("C:/Users/ttu/Documents/ROI_Grid/ROI_Merged")

deduplicate <<- c()
deduplicate_search <<- c()

for (i in 1:length(all_files)) {
  duplicate <- read.delim(all_files[i], header = TRUE, stringsAsFactors = FALSE)
  duplicate <- special(duplicate)
  #duplicate <- special_2(duplicate) #Comment out if there is no () in Sample Name
  duplicate_ID <- unique(duplicate$Sample.Name)
  duplicate_ID_search <- unique(duplicate$Sample_ID_Extracted2)
  deduplicate[i] <- duplicate_ID
  deduplicate_search[i] <- duplicate_ID_search
}

unique_deduplicate <- unique(deduplicate_search)
final_list_index <<- c()

for (j in 1:length(unique_deduplicate)) {
  location <- which(grepl(unique_deduplicate[j], deduplicate))
  final_list_index[j] <- location[1]
}

file_copy <- c()
file_copy2 <- c()

for (k in 1:length(final_list_index)) {
  file_copy[k] <- all_files[final_list_index[k]]
}

for (l in 1:length(file_copy)) {
  file_copy2[l] <- paste("C:/Users/ttu/Documents/ROI_Grid/ROI_Merged/", file_copy[l], sep = "")
}

file.copy(from = file_copy2, "C:/Users/ttu/Documents/ROI_Grid/Deduplicated/")

setwd("C:/Users/ttu/Documents/ROI_Grid/Deduplicated")
all_deduplicated_files <- list.files("C:/Users/ttu/Documents/ROI_Grid/Deduplicated")

dataset <- do.call("rbind",lapply(all_deduplicated_files, FUN=function(files){ 
  fread(files, header=TRUE, sep="\t")}))

write.table(dataset, file = "Pt7p_cell_seg_summary_6%.txt", sep = "\t", row.names = FALSE)

