#NN_Cluster_Image.R#
#Travis Tu; Last Edited: 8/8/2018
#Folders: NN, Images, Files, Results

#Install and load required packages
required = c("RANN", "tiff", "stringr", "gtools", "raster", "magick")
for (lib in required)
{
  if (!require(lib, character.only=TRUE))
  {
    install.packages(lib, repos="http://cran.rstudio.com")
    suppressMessages(library(lib, character.only=TRUE, quietly=TRUE))
  }
}

library(RANN)
library(tiff)
library(stringr)
library(gtools)
library(raster)
library(magick)

#Function to match image files and their cell seg files
#Input: String that is the name of the image file, list of cell seg files
#Extract out the coordinates and the root and match them to the cell seg file
#Output: Index of the cell seg file that matches the image file
file_match <- function(image_name, file_list) {
  coor <- str_extract_all(image_name, "\\[([^\\[\\]]+)\\]")[[1]]
  root <- gsub("_\\[.*(.*)", "\\1", image_name)
  
  for (i in 1:length(file_list)) {
    if (grepl(coor, file_list[i]) == TRUE & 
        grepl(root, file_list[i], ignore.case = TRUE) == TRUE) {
      a <- grep(coor, file_list)
      b <- grep(root, file_list, ignore.case = TRUE)
      index <- which(a %in% b)
      return(index)
    }
  }
  return(NA)
}

#Function to identify cells that are within 25 microns of cell of interest
#Input: Dataframe of cell seg data, phenotypes of interest
#Phenotype1 is the phenotype of cell of interest, phenotype2 is phenotype of target
#Uses RANN to get row number of cells within 25 microns of cell of interest
#Output: Matrix with rows corresponding to cells of interest and
#Columns corresponding to nearest neighbors within 25 microns (50 pixels)
NN <- function(dataframe, phenotype1, phenotype2) {
  data <- dataframe[dataframe$Phenotype == phenotype1, ]
  #data <- dataframe[dataframe$Phenotype != phenotype, ]
  query <- dataframe[dataframe$Phenotype == phenotype2, ]
  
  X <- which(grepl("Cell.X.Position", names(dataframe)))
  Y <- which(grepl("Cell.Y.Position", names(dataframe)))
  names <- c(X, Y)
  
  nearest_neighbor <- nn2(query[,names], data[,names], treetype = "kd", 
                          searchtype = "radius", radius = 50)
  return(nearest_neighbor)
}

#Mergesort#
merge <- function(left, right) {
  merged <- c()
  while(length(left) > 0 && length(right) > 0) { #Run as long as BOTH subarrays are not empty
    if (left[1] <= right[1]) { #compare first element of each subarray with each other
      merged <- c(merged, left[1]) #Add the smaller one into the newly created array and delete it from its subarray
      left <- left[-1]
    }
    else {
      merged <- c(merged, right[1])
      right <- right[-1]
    }
  }
  if(length(left) > 0) { #As soon as one  of the subarrays runs out, add on remainder of that subarray to the array
    merged <- c(merged, left)
  }
  if (length(right) > 0) {
    merged <- c(merged, right)
  }
  return(merged)
}

merge_sort_R <- function(array) {
  if (length(array) <= 1) { #base case (array divded into single elements)
    return(array)
  }
  else {
    midpoint <- floor(length(array)/2) #take midpoint of the array to split array into two subarrays
    left_subarray <- array[1:midpoint]
    right_subarray <- array[(midpoint + 1): length(array)]
    
    left_subarray <- merge_sort_R(left_subarray) #repeat split on each subarray until they are only one element long each
    right_subarray <- merge_sort_R(right_subarray)
    
    #if(left_subarray[length(left_subarray)] <= right_subarray[1]) { #Because each subarray has been sorted, if the last element of the left
    #is smaller than the first and smallest element in the right subarray, the two can just be merged
    #return(c(left_subarray, right_subarray))
    #}
    #else {
    return(merge(left_subarray, right_subarray)) #combine subarrays back together as described
    #} 
  }
} 

#Function to identify duplicates in a vector
#Input: vector containing duplicates
#Output: a vector a numbers that were duplicates in the input vector
id_duplicates <- function(vector) {
  sorted_vector <- merge_sort_R(vector)
  duplicate_vector <- c()
  a <- 1
  for (i in 2:length(sorted_vector)) {
    if (sorted_vector[i-1] == sorted_vector[i]) {
      if (length(duplicate_vector) == 0 || 
          duplicate_vector[length(duplicate_vector)] != sorted_vector[i-1]) {
        duplicate_vector[a] <- sorted_vector[i-1]
        a <- a + 1
      }
    }
  }
  return(duplicate_vector)
} 

#Function to find overlapping nearest neighbors and assign it to nearest neighbor
#Input: Output of NN (NN1 = NN[[1]] and NN2 = NN[[2]])
#Output: Revised version of NN[[1]] with overlapping nearest neighbors only
#assigned to its nearest neighbor
no_duplicates <- function(NN1, NN2) {
  dimension <- dim(NN1)
  #Convert the indices matrix into a list#
  index_list_og <- unlist(as.list(NN1))
  distance_list <- unlist(as.list(NN2))
  
  for (i in 1: length(index_list_og)) {
    if(index_list_og[i] == 0) { #If 0 is a unique cell index, replace it with NA
      index_list_og[i] <- NA
    }
  }
  index_list <- Filter(Negate(is.na), index_list_og) #Remove all NAs
  unique_cell <- id_duplicates(index_list)
  
  indices_list <- list()
  for (j in 1:length(unique_cell)) {
    indices <- which(index_list_og %in% unique_cell[j])
    indices_list[[j]] <- indices
  }
  names(indices_list) <- unique_cell
  #print(indices_list)
  
  distances_list <- indices_list
  for (k in 1:length(indices_list)) {
    for (m in 1:length(indices_list[[k]])) {
      distances_list[[k]][m] <- distance_list[indices_list[[k]][m]]
    }
  }
  
  for(k in 1:length(distances_list)) {
    save_min <- which.min(distances_list[[k]])
    for (m in 1:length(indices_list[[k]])) {
      if (m != save_min) {
        index_list_og[indices_list[[k]][m]] <- 0
      }
    }
  }
  
  #print(matrix(index_list_og, nrow = dimension[1], ncol = dimension[2]))
  
  na_vector <- is.na(index_list_og)
  for (p in 1:length(index_list_og)) {
    if (na_vector[p] == TRUE) {
      index_list_og[p] <- 0
    }
  }
  
  return(matrix(index_list_og, nrow = dimension[1], ncol = dimension[2]))
}

#Function to create a dataframe for plotting line segments
#Input: Deduplicated Nearest Neighbor Matrix, dataframes of only Phenotype 1 and 2
#Use row numbers to identify coordinates of the points
#Also takes in an image to get its dimensions so it can avoid flipping at the end
#Output: Dataframe with X and Y coordinates and start and endpoints (4 columns)
line_seg_df <- function(matrix, set1, set2, y) {
  startx <- c()
  starty <- c()
  endx <- c()
  endy <- c()
  
  b <- 1
  for (i in 1:nrow(matrix)) {
    if (sum(matrix[i, ]) > 0) {
      sorted_row <- sort(matrix[i, ], decreasing = TRUE)
      a <- 1
      while(sorted_row[a] > 0 & a <= ncol(matrix)) {
        endx[b] <- set2[sorted_row[a], 9]
        endy[b] <- set2[sorted_row[a], 10]
        startx[b] <- set1[i, 9]
        starty[b] <- set1[i, 10]
        a <- a + 1
        b <- b + 1
      }
    }
  }
  
  startx_df <- as.data.frame(startx)
  #starty_df <- as.data.frame(y - starty)
  starty_df <- as.data.frame(starty)
  endx_df <- as.data.frame(endx)
  #endy_df <- as.data.frame(y - endy)
  endy_df <- as.data.frame(endy)
  
  final_df <- cbind(startx_df, starty_df, endx_df, endy_df)
  return(final_df)
}

setwd("C:/Users/ttu/Documents/NN/Images")
image_list <- list.files("C:/Users/ttu/Documents/NN/Images")
file_list <- list.files("C:/Users/ttu/Documents/NN/Files")

for (i in 1:length(image_list)) {
  index <- file_match(image_list[i], file_list)
  print(file_list[index])
  print(image_list[i])
  setwd("C:/Users/ttu/Documents/NN/Files")
  cell_seg <- read.delim(file_list[index], header = TRUE)
  phenotypes <- sort(as.vector(unique(cell_seg$Phenotype)))
  perm <- permutations(n = length(phenotypes), r = 2, v = phenotypes)
  
  for (j in 1:nrow(perm)) {
    cell_matrix <- NN(cell_seg, perm[j, 1], perm[j, 2])
    #print(cell_matrix)
    set1 <- cell_seg[cell_seg$Phenotype == perm[j, 1], ]
    set2 <- cell_seg[cell_seg$Phenotype == perm[j, 2], ]
    
    deduplicated_matrix <- no_duplicates(cell_matrix[[1]], cell_matrix[[2]])
    
    line_seg_coor <- line_seg_df(deduplicated_matrix, set1, set2)
    
    setwd("C:/Users/ttu/Documents/NN/Images")
    image_var <- image_read(image_list[i])
    image_var <- image_flip(image_var)
    tiff(filename = paste("C:/Users/ttu/Documents/NN/Results_Flipped/", 
                          perm[j, 2], "_within_25_microns_of_", perm[j, 1], "_", image_list[i], ".tif", sep = ""),
         width = 4, height = 4, units = 'in', res = 300)

    plot(image_var)
    #points(set1[, 9], set1[, 10], col = "white")
    #plot(set1[, 9], set1[, 10], col = "blue", ylim = rev(range(c(0:1008))),pch = 20)
    #points(set1[, 9], set1[, 10], col = "white",pch = 20)
    #points(set2[, 9], set2[, 10], col = "pink", pch = 20)
    
    #Add line segments between query cells and its nearest neighbors#
    for (k in 1:nrow(line_seg_coor)) {
      segments(line_seg_coor[k,1], line_seg_coor[k,2], line_seg_coor[k,3], line_seg_coor[k,4], 
               lwd = 1, col = "white")
    }
    #scalebar(200, below = "Pixels", xy = c(100, 100))
    dev.off()
  }
}

setwd("C:/Users/ttu/Documents/NN/Results_Flipped")
to_flip <- list.files("C:/Users/ttu/Documents/NN/Results_Flipped")

for (i in 1:length(to_flip)) {
  image_var <- image_read(to_flip[i])
  final_image <- image_flip(image_var)
  image_write(final_image, path = paste("C:/Users/ttu/Documents/NN/Final_Results/", 
                                        to_flip[i] , sep = ""), format = "tiff")
}



