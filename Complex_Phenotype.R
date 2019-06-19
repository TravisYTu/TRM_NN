#Complex_Phenotype.R#
#Travis Tu; Last Edited: 8/15/2018#
#Folders: Complex_Phenotype, Cell_Seg, Results

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

#Recursive function to create multi-threshold complex phenotypes
#Input: Dataframe, column, threshold, split_list, counter
#Dataframe = cell seg file
#column = string denoting the column name of interest
#threshold = vector of thresholds to use to create phenotypes
#split_list = an empty list that will be filled and returned
#counter = a count recursions to change split_list index
multi_complex <- function(dataframe, column, threshold, split_list, counter) {
  sorted_threshold <- merge_sort_R(threshold)
  index <- which(grepl(column, names(dataframe)))
  
  if (length(sorted_threshold) == 1) {
    split1 <- subset(dataframe, dataframe[,index] <= sorted_threshold[1])
    split2 <- subset(dataframe, dataframe[,index] > sorted_threshold[1])
    split_list[[counter]] <- split1$Cell.ID
    split_list[[counter + 1]] <- split2$Cell.ID
    print(split_list)
    return(split_list)
  }
  else {
    split1 <- subset(dataframe, dataframe[,index] <= sorted_threshold[1])
    split2 <- subset(dataframe, dataframe[,index] > sorted_threshold[1])
    split_list[[counter]] <- split1$Cell.ID
    return(multi_complex(split2, column, sorted_threshold[-1], split_list, (counter + 1)))
  }
}

#Binary Search#
search <- function(array, targetValue = NULL) {
  min <- 1
  max <- length(array)
  if (length(array) == 0 || is.null(targetValue) == TRUE) {
    return(NA)
  }
  else {
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
    return(NA)
  }
}

#Function to concatenate original phenotype and name together to make complex phenotype
#Input: dataframe, array, name
#dataframe = cell seg data with sorted cell ids
#array = vector of indices that need to be changed
#name = name to concatenate onto the original phenotype name
#Note: uses binary search function
phenotype_fix <- function(dataframe, array, name) {
  file_Cell_ID <- which(grepl("Cell.ID", names(dataframe)))
  for (i in 1:length(array)) {
    search_result <- search(dataframe[, file_Cell_ID], array[i])
    file_Phenotype <- which(grepl("^Phenotype$", names(dataframe)))
    og_name <- dataframe[search_result, file_Phenotype]
    dataframe[search_result, file_Phenotype] <- paste(og_name, name, sep = "")  
  }
  return(dataframe)
}

#Function to convert phenotype column into complex phenotypes using thresholds
#Input: Dataframe = cell seg file
#column = list of strings denoting column names
#threshold = list of floats denoting thresholds for +/-
#name = list of names to change phenotype with
#Note: threshold, column, and name indices should match
complex <- function(dataframe, column, threshold, name) {
  final_list <- list()
  split_list <- list()
  dataframe$Cell.ID <- 1:nrow(dataframe)
  for (i in 1:length(threshold)) {
    final_list[[i]] <- multi_complex(dataframe, column[[i]], threshold[[i]], split_list, 1)
  }
  print(final_list)
  
  for (m in 1:length(final_list)) {
    for (n in 1:length(final_list[[m]])) {
      if (length(final_list[[m]][[n]]) == 0) {
        n <- n + 1
        next
      }
      else {
        dataframe <- phenotype_fix(dataframe, final_list[[m]][[n]], name[[m]][n])
      }
    }
  }
  return(dataframe)
}

#Function to check if column, threshold, and name have the same dimensions
check <- function(column, threshold, name) {
  if (length(column) != length(threshold)||
      length(column) != length(name)||
      length(threshold) != length(name)) {
    stop("Parameters are not the same length")
  }
  else {
    for (i in 1:length(threshold)) {
      if (length(threshold[[i]]) != (length(name[[i]]) - 1)) {
        stop("Parameters within threshold or name are not same length")
      }
    }
  }
}

setwd("C:/Users/ttu/Documents/Complex_Phenotype/Cell_Seg")
cs_files <- list.files("C:/Users/ttu/Documents/Complex_Phenotype/Cell_Seg")

for (i in 1:length(cs_files)) {
  dataframe_og <- read.delim(cs_files[i], header = TRUE, stringsAsFactors = FALSE)
  #If you only want to create complex phenotypes on a specific phenotype
  #i.e. only want CD8+PD-1+ and not CK+PD-1+
  #dataframe_pos <- dataframe_og[dataframe_og$Phenotype == "Cancer", ]
  #dataframe_neg <- dataframe_og[dataframe_og$Phenotype != "Cancer", ]
  
  column <- list("Membrane.Opal.650.Mean..Normalized.Counts..Total.Weighting.",
                 "Membrane.Opal.570.Mean..Normalized.Counts..Total.Weighting.")
  threshold <- list(0.5, c(0.1, 0.2))
  name <- list(c("Test1-", "Test1+"), c("Test2lo", "Test2int", "Test2high"))
  #Names should be listed in increasing order (- before + for single threshold)
  #(low before intermediate before high in multi-thresholds)
  
  check(column, threshold, name)
  result <- complex(dataframe_og, column, threshold, name)
  
  #If you only want to create complex phenotypes on a specific phenotype
  #i.e. only want CD8+PD-1+ and not CK+PD-1+
  #result <- complex(dataframe_pos, column, threshold, name)
  #result <- rbind(dataframe_neg, result)
  
  write.table(result, file = paste("C:/Users/ttu/Documents/Complex_Phenotype/Results/Complex_",
                                   cs_files[i], sep =""), sep = "\t", row.names = FALSE)
}
