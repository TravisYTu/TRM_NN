#Phenotype_Distance_TT.R#
#Travis Tu; Last Edited: 10/25/2017#
#Folders: Distance, Files, Results#

#Install and load required packages
required = c("RANN", "data.table")
for (lib in required)
{
  if (!require(lib, character.only=TRUE))
  {
    install.packages(lib, repos="http://cran.rstudio.com")
    suppressMessages(library(lib, character.only=TRUE, quietly=TRUE))
  }
}

library(RANN)
library(data.table)

#Function to create the nearest neighbor matrix from a cell to its nearest neighbor of a certain phenotype
#Input: dataframe, phenotype
#dataframe = cell_seg_data (recoordinated)
#phenotype = string of the phenotype of interest
#Output: A list with two elements. The first is a matrix with every row corresponding to a cell in dataframe.
#The column is its nearest neighbor of the phenotype of interest. The second is a matrix with the same rows
#as the first one but the column is the distance in pixels to its nearest neighbor of the phenotype of interest.
NN <- function(dataframe, phenotype) {
  data <- dataframe
  #data <- dataframe[dataframe$Phenotype != phenotype, ]
  query <- dataframe[dataframe$Phenotype == phenotype, ]
  
  X <- which(grepl("^Cell.X.Position$", names(dataframe)))
  Y <- which(grepl("^Cell.Y.Position$", names(dataframe)))
  print(X)
  print(Y)
  names <- c(X, Y)
  
  nearest_neighbor <- nn2(query[,names], data[,names], k=1, treetype = "kd", searchtype = "standard")
  #nearest_neighbor <- nn2(query[,names], data[,names], k= 1, treetype = "kd", searchtype = "radius", radius = 50)
  
  return(nearest_neighbor)
}

distance <- function(x1, y1, x2, y2) {
  return(sqrt(((x2- x1)^2) + ((y2- y1)^2))) 
}

setwd("C:/Users/ttu/Documents/Distance/Files")

files_to_dist <- list.files("C:/Users/ttu/Documents/Distance/Files")

for (i in 1:length(files_to_dist)) {
  dist <- read.delim(files_to_dist[i], header = TRUE)
  
  final_dist <- c()
  #final_dist_id <- c()
  unique_phenotype <- as.vector(unique(dist$Phenotype))
  
  a <- 1
  for (j in 1:length(unique_phenotype)) {
    NN_matrix <- NN(dist, unique_phenotype[j])
    nearest_neighbor <- as.vector(NN_matrix[2])
    #id <- as.vector(NN_matrix[1])
    final_dist[a] <- nearest_neighbor
    #final_dist_id[a] <- id
    a <- a + 1
  }
  
  final_dist_df <- as.data.frame(final_dist)
  #final_dist_id_df <- as.data.frame(final_dist_id)
  names(final_dist_df) <- unique_phenotype
  
  for (k in 1:length(names(final_dist_df))) {
    names(final_dist_df)[k] <- paste("Distance to ", names(final_dist_df)[k], sep = "")
  }
  
  distanced <- cbind(dist, final_dist_df)
  #distanced <- cbind(dist, final_dist_df, final_dist_id_df)
  write.table(distanced, file = paste("C:/Users/ttu/Documents/Distance/Results/Dist_", files_to_dist[i], sep =""),
         row.names = FALSE, col.names = TRUE, sep = "\t")
}