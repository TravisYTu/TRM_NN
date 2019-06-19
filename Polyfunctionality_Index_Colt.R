#Polyfunctionality_Index_Colt.R#
#Travis Tu; Last Edited: 4/4/2017#
#Files: Polyfunctionality_Index, Files#

#Make sure Split folder is empty before running

setwd("C:/Users/ttu/Documents/Polyfunctionality_Index/Files")

Spice_list <- list.files("C:/Users/ttu/Documents/Polyfunctionality_Index/Files")

Polyfunctionality_Index <- function(Fi, n, q){
  PI <- 0
  for (i in 1:n+1){
    PI <- PI + (((i-1)/n)^q)*Fi[[i]]
  }
  return(PI)
} 

Split <- function(dataframe, comparison){
    split <- subset(dataframe, dataframe$Subject == comparison)
    write.table(split, file = paste("C:/Users/ttu/Documents/Polyfunctionality_Index/Split/", 
                                    toString(comparison), ".txt", sep = ""), 
        sep = "\t", row.names = FALSE)
}

convert <- function(dataframe){
  for (i in 1:nrow(dataframe)){
    for (j in 1:ncol(dataframe)){
      if(dataframe[i, j] == "+"){
        dataframe[i,j] <- 1
      }
      else if(dataframe[i,j] == "-"){
        dataframe[i,j] <- 0
      }
    }
  }
  return(dataframe)
}


Positives <- function(dataframe, cytokine1, cytokine2, cytokine3, value){
  c1_index <- which(grepl(cytokine1, names(dataframe)))
  c2_index <- which(grepl(cytokine2, names(dataframe)))
  c3_index <- which(grepl(cytokine3, names(dataframe)))
  value_index <- which(grepl(value, names(dataframe)))
  
  temp1 <- 0
  temp2 <- 0
  
  storage <- list()
  
  for(i in 1:nrow(dataframe)){
    if (sum(dataframe[i,c1_index], dataframe[i, c2_index], dataframe[i, c3_index])== 3){
      storage[[4]] <- dataframe[i, value_index]
    }
    else if (sum(dataframe[i,c1_index], dataframe[i, c2_index], dataframe[i, c3_index])== 1){
      temp1 <- temp1 + dataframe[i, value_index]
      storage[[2]] <- temp1
    }
    else if (sum(dataframe[i,c1_index], dataframe[i, c2_index], dataframe[i, c3_index])== 0)
      storage[[1]] <- dataframe[i, value_index]
    else {
      temp2 <- temp2 + dataframe[i, value_index]
      storage[[3]] <- temp2
    }
  }
  return(storage)
}

for (file in Spice_list){
  Spice <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  converted_Spice <- convert(Spice)
  
  unique_subject <<- unique(converted_Spice$Subject)
  for (j in 1:length(unique_subject)){
    Split(converted_Spice, unique_subject[j])
  } 
}

setwd("C:/Users/ttu/Documents/Polyfunctionality_Index/Split")

Split_list <- list.files("C:/Users/ttu/Documents/Polyfunctionality_Index/Split")

final <- as.data.frame(matrix(0, ncol = 2, nrow = length(Split_list)))
colnames(final) <- c("Subject", "Polyfunctionality Index")

a <- 1
for (file in Split_list) {
  Split <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  
  storage_list <- Positives(Split, cytokine1= "IFNy", cytokine2 = "IL.2", cytokine3 = "TNFa", value = "Value")
  
  PI_value <- Polyfunctionality_Index(storage_list, 3, 1)
  final[a, 2] <- PI_value
  final[a, 1] <- Split[1,1]
  a <- a + 1
}

write.table(final, file = "C:/Users/ttu/Documents/Polyfunctionality_Index/New_melanoma_12_13_14_PI.txt", sep = "\t", row.names = FALSE)
