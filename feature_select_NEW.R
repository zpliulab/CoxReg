library(dplyr)       
library(tidyr)
library(tidyverse)   


myfile = list.files("Independent_data")                
dir = paste("./Independent_data/", myfile, sep="")    
n = length(dir)                                  


Data1 = read.table(file = dir[1],header=T, check.names = FALSE) 
dim(Data1)
coef_gene <- read.csv("univariate_cox_coef.csv", header = T, sep=',')
gene <- as.matrix(coef_gene[,1])


colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))


genedata <- merge(gene, Data2, by = "gene")
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(Data1[c(1,2),],genedata1)   


name <- dir[1]

## Press ./Split, Extract
name1 <- str_split_fixed(name, "./", 2);
name2 <- str_split_fixed(name1[2], "/", 2)
name3 <- str_split_fixed(name2[2], "[.]", 2);
name4 <- substr(name3[1], 1, 9)
name5 <- substr(name3[1], 23, 27)
name6 <- str_c(name4, name5)

path <- paste("./Feature_data/Data_GEO/",paste(name6,".txt", sep=""), sep="")
write.table(genedata2, path, quote = F)


myfile = list.files("Independent_data")              
dir = paste("./Independent_data/", myfile, sep="")   
n = length(dir)                                 


for (i in 2:n){

  # i <-2
  Data1 = read.table(file = dir[i],header=T, check.names = FALSE) 
  dim(Data1)
  
  coef_gene <- read.csv("univariate_cox_coef.csv", header = T, sep=',')
  gene <- as.matrix(coef_gene[,1])

  colnames(gene) <- c('gene')
  Data2 <- cbind(rownames(Data1), Data1)
  colnames(Data2) <- c('gene', colnames(Data1))
  
  genedata <- merge(gene, Data2, by = "gene")
  genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
  genedata2 <- rbind(Data1[c(1,2),],genedata1)  
  
  name <- dir[i]
 
  name1 <- str_split_fixed(name, "./", 2);
  name2 <- str_split_fixed(name1[2], "/", 2)
  name3 <- str_split_fixed(name2[2], "[.]", 2);
  name4 <- substr(name3[1], 1, 9)
  name5 <- substr(name3[1], 23, 27)
  name6 <- str_c(name4, name5)
  
  
  path <- paste("./Feature_data/Data_GEO/",paste(name6,".txt", sep=""), sep="")
  write.table(genedata2, path, quote = F)

  
  print("***i***")
  print(i)

}                                               






