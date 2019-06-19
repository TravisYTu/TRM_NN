#TRM_Spatial.R#
#Travis Tu; Last Edited: 7/5/2018#
#Folders: TRM, Files, Bins#

#Install and load required packages
required = c("plyr", "tiff")
for (lib in required)
{
  if (!require(lib, character.only=TRUE))
  {
    install.packages(lib, repos="http://cran.rstudio.com")
    suppressMessages(library(lib, character.only=TRUE, quietly=TRUE))
  }
}


library(plyr)
library(tiff)

#Function good_or_bad takes a vector of names (strings) and determines if
#They are good or bad outcome based on the name of the file
#Patients 1-19 are bad, 20-40 are good
#Returns a list length 2 (good and bad file names)
good_or_bad <- function(vector_names) {
  gob <- list()
  g <- c()
  b <- c()
  a <- 1
  d <- 1
  for (i in 1:length(vector_names)) {
    class <- gsub("_.*", "\\1", vector_names[i])
    if (class == "Bad") {
      b[a] <- vector_names[i]
      a <- a + 1
    }
    else {
      g[d] <- vector_names[i]
      d <- d + 1
    }
  }
  gob[[1]] <- g
  gob[[2]] <- b
  return(gob)
}

#Quicksort for R#
quicksort_R <- function(array) {
  if(length(array) <= 1) { #base case
    return(array)
  }
  else {
    pivot <- array[1] #Take first element of the array
    subarray <- array[-1] #Take rest of array
    
    subarray_low <- subarray[subarray < pivot] #new array with all numbers less than pivot
    subarray_high <- subarray[subarray >= pivot] #new array with all numbers greater than pivot
    
    subarray_low <- quicksort_R(subarray_low) #Repeat with subarrays until base case reached
    subarray_high <- quicksort_R(subarray_high)
    
    return(c(subarray_low, pivot, subarray_high)) #stich vector back together
  }
}

#Function bins takes a dataframe and vector of cutoffs and separates the 
#Dataframe into subsets by Distance to CK column cutoffs
#Output a list of length vector of cutoffs with subsetted dataframe
bins <- function(dataframe, cutoffs, split_list, counter) {
  sorted_cutoffs <- quicksort_R(cutoffs)
  if (length(cutoffs) == 1) {
    split_list[[counter]] <- nrow(subset(dataframe, dataframe$Distance.to.CK < sorted_cutoffs[1]))
    split_list[[counter + 1]] <- nrow(subset(dataframe, dataframe$Distance.to.CK >= sorted_cutoffs[1]))
    return(split_list)
  }
  
  else {
    split <- subset(dataframe, dataframe$Distance.to.CK < sorted_cutoffs[1])
    split2 <- subset(dataframe, dataframe$Distance.to.CK >= sorted_cutoffs[1])
    split_list[[counter]] <- nrow(split)
    return(bins(split2, sorted_cutoffs[-1], split_list, (counter+1)))
  }
}

#Function to split the dataframe by phenotype column
#Input: dataframe, Output: list of dataframes by phenotype
phenotype_list <- function(dataframe) {
  p_list <- list()
  phenotypes <- unique(dataframe$Phenotype)
  phenotypes <- phenotypes[order(phenotypes)]
  for (i in 1:length(phenotypes)) {
    p_list[[i]] <- subset(dataframe, dataframe$Phenotype == phenotypes[i])
  }
  names(p_list) <- phenotypes
  return(p_list)
}

#Function to split the dataframe by phenotype column
#Input: dataframe, Output: list of dataframes by phenotype
phenotype_list2 <- function(dataframe) {
  p_list <- list()
  phenotypes <- unique(dataframe$Phenotype)
  phenotypes <- phenotypes[order(phenotypes)]
  for (i in 1:length(phenotypes)) {
    temp <- subset(dataframe, dataframe$Phenotype == phenotypes[i])
    p_list[[i]] <- temp$Distance.to.CK
  }
  names(p_list) <- phenotypes
  return(p_list)
}

#Function to add up the components of a list
#Input: list; Output: sum of the list's components
sum_list <- function(list) {
  for (i in 1:length(list)) {
    if (i == 1) {
      total <- list[[i]]
    }
    else {
      total <- total + list[[i]]
    }
  }
  return (total)
}

setwd("C:/Users/ttu/Documents/TRM/Files")
TRM_files <- list.files("C:/Users/ttu/Documents/TRM/Files")

TRM_files_gob <- good_or_bad(TRM_files)
cutoffs <- c(25, 50, 100)

final_list <- list()

for (i in 1:length(TRM_files_gob)) {
  for (j in 1:length(TRM_files_gob[[i]])) {
    split_list <- list()
    if (is.null(TRM_files_gob[[i]][j]) == TRUE) {
      next
    }
    else {
      file <- read.delim(TRM_files_gob[[i]][j], header = TRUE)
      p_names <- unique(file$Phenotype)
      p_names <- p_names[order(p_names)]
      p_list <- phenotype_list(file)
      
      subset_list <- list()
      for (l in 1:length(p_list)) {
        subset_list[[l]] <- bins(p_list[[l]], cutoffs, split_list, 1)
      }
      
      for (k in 1:length(subset_list)) {
        denom <- sum_list(subset_list[[k]])
        for (m in 1:length(subset_list[[k]])) {
          subset_list[[k]][[m]] <- subset_list[[k]][[m]]/ denom
        }
      }
      
      final <- data.frame(matrix(0, nrow = length(p_list), ncol = (length(cutoffs) + 1)))
      rownames(final) <- p_names
      colnames(final) <- c("0", cutoffs)
      for (a in 1:length(subset_list)) {
        for (b in 1:length(subset_list[[a]])) {
          final[a, b] <- subset_list[[a]][[b]]
        }
      }
    }
    write.table(final, file = paste("C:/Users/ttu/Documents/TRM/Bins/", TRM_files_gob[[i]][j], 
                                    sep = ""), sep = "\t")
  }
}

#########################################################################################################
#Code to create density plots

gob_merged_list <- list()
for (i in 1:length(TRM_files_gob)) {
  dataset <- do.call("rbind.fill",lapply(TRM_files_gob[[i]], FUN=function(files){ 
    read.delim(files, header=TRUE, sep="\t")}))
  gob_merged_list[[i]] <- dataset
}

gob_phenotype_list <- list()
for (i in 1:length(gob_merged_list)) {
  gob_phenotype_list[[i]]  <- phenotype_list2(gob_merged_list[[i]])
}

g_CD8_CD103x <- density(gob_phenotype_list[[1]][[1]])
g_CD8_CD103 <- density(gob_phenotype_list[[1]][[2]])
b_CD8_CD103x <- density(gob_phenotype_list[[2]][[1]])
b_CD8_CD103 <- density(gob_phenotype_list[[2]][[2]])

tiff(filename = "C:/Users/ttu/Documents/TRM/Density_Plot/Density_plot.tif", 
     width = 4, height = 4, units = 'in', res = 300)
plot(g_CD8_CD103x, xlim = c(0, 100) ,ylim = c(0,0.04), main = "TRM Density Plot", 
     xlab = "Distance to CK (pixels)")
lines(g_CD8_CD103, col = "red", lty = 2)
lines(b_CD8_CD103x, col = "blue", lty = 3)
lines(b_CD8_CD103, col = "green", lty = 4)
legend(60, 0.04, legend = c("Good CD103-", "Good CD103+", "Bad CD103-", "Bad CD103+"),
       col = c("black", "red", "blue", "green"), lty = 1:4, cex = 0.5)
dev.off()

