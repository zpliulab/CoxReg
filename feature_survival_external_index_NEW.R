

rm(list = ls())

library(tidyverse)


source('D:\\E\\博士\\R_程序\\BRCA\\R\\myoverlap_separate2GroupsCox.R')


# 导入数据
setwd("D:\\E\\博士\\R_程序\\BRCA\\Feature_data")


myfile = list.files("Data_GEO")                #list.files命令将input文件夹下所有文件名输入a
dir = paste("./Data_GEO/", myfile, sep="")     #用paste命令构建路径变量dir
n = length(dir) 


data = read.table(file = dir[1],header=T, check.names = FALSE) #读入第一个文件内容（可以不用先读一个，但是为了简单，省去定义data.frame的时间，我选择先读入一个文件

datat <- as.data.frame(t(data))
Data <- as.matrix(data)

x <- model.matrix(status ~., datat)[,-c(1,2)]   # 删除了time, statuss~.
x_hat <- data.frame(x)
y <- as.matrix(datat[,c(1,2)])  
colnames(y) <- c("time", "status")
y_hat <- data.frame(y)

## 系数
setwd("D:\\E\\博士\\R_程序\\BRCA\\Data\\TCGA_NEW\\result")
coef_gene <- read.csv("univariate_cox_coef.csv", header = T, sep=',')
coef <- as.matrix(coef_gene[,2])
rownames(coef) <- coef_gene[,1]


# my_overlap_coef ------------------------------------------------------
my_overlap <- function(x, y){
  # x <- coef
  # y <- Data
  
  coefs.v <- x[,1] %>% { .[. != 0]}
  coefs.v %>% {
    data.frame(gene.name   = names(.),
               coefficient = .,
               stringsAsFactors = FALSE)
  } %>%
    arrange(gene.name) %>%
    knitr::kable()
  
  sele <- rownames(as.matrix(coefs.v))
  gene <- rownames(y)[-c(1,2)]
  overlap <- intersect(sele, gene)
  
  lab <- x[,1] %>% { .[. != 0]} %>% as.matrix
  coefs.v <- lab[overlap,]
  
  my <- list(coefs.v, overlap)
  return(my)
} 



# 画图 ----------------------------------------------------------------------

library(ggpubr)
library(magrittr)
library(survminer)
# library(glmSparseNet)
coef_test <- my_overlap(coef, Data)
plotp_Train <- separate2GroupsCox(as.vector(coef_test[[1]]), x_hat[, coef_test[[2]]], 
                                 plot.title = 'GSE1456_smfs', as.data.frame(y), 
                                 legend.outside = T)
plot_train <- plotp_Train$plot
# plot_train


## for Xtile
p_index <- cbind(y,plotp_Train$index)
colnames(p_index) <- c(colnames(y), "riskscore")



name <- dir[1]

## 按./拆分,提取
name1 <- str_split_fixed(name, "./", 2);
name2 <- str_split_fixed(name1[2], "/", 2)
name3 <- str_split_fixed(name2[2], "[.]", 2);
name4 <- name3[1]
name5 <- str_c(name4,"_ex")
# class(name5)

# out <- rbind(p_train, CI_train)
# colnames(out) <- name4
  

setwd("D:\\E\\博士\\R_程序\\BRCA\\Feature_data")
# path <- paste("./P_CI/",paste(name4,".txt", sep=""), sep="")
# write.table(out, path, quote = F)  
# path <- paste("./Sangerbox/",paste(name4,".csv", sep=""), sep="")
# write.csv(sanger, path, row.names = F)   
# path <- paste("./Xtile/",paste(name4,".csv", sep=""), sep="")
# write.csv(p_index, path, quote = F, row.names = F)  # 需要手动转换为tab分隔的txt
path <- paste("./Xtile/",paste(name4,".txt", sep=""), sep="")
write.table(p_index, path, quote = F, row.names = F, sep="\t")
path <- paste("./Figure/", paste(name5,".pdf", sep=""), sep="")
pdf(path, width = 4.5, height = 4.5, onefile = FALSE)  # 设置onefile 参数为FALSE 后，散点图会覆盖前面的空白
plot_train
dev.off() 



###########################################################################

setwd("D:\\E\\博士\\R_程序\\BRCA\\Feature_data")


myfile_GEO = list.files("Data_GEO")                #list.files命令将input文件夹下所有文件名输入a
dir_GEO = paste("./Data_GEO/", myfile_GEO, sep="")     #用paste命令构建路径变量dir
n = length(dir_GEO) 


## 循环的pdf出问题，其他正常
for (i in 1:n){
  
  # i <- 18

  data = read.table(file = dir_GEO[i],header=T, check.names = FALSE) #读入第一个文件内容（可以不用先读一个，但是为了简单，省去定义data.frame的时间，我选择先读入一个文件
  datat <- as.data.frame(t(data))
  Data <- as.matrix(data)
  dim(Data)
  # View(Data[,1:5])
  
  # status <- as.name(colnames(datat)[2])
  
  x <- as.matrix(datat[,-c(1,2)])   # 删除了time, statuss~.
  colnames(x) <- rownames(Data)[-c(1,2)]
  x_hat <- data.frame(x)
  y <- as.matrix(datat[,c(1,2)])  
  colnames(y) <- c("time", "status")
  y_hat <- data.frame(y)
  
  ## 系数
  setwd("D:\\E\\博士\\R_程序\\BRCA\\Data\\TCGA_NEW\\result")
  coef_gene <- read.csv("univariate_cox_coef.csv", header = T, sep=',')
  coef <- as.matrix(coef_gene[,2])
  rownames(coef) <- coef_gene[,1]
  

  coef_test <- my_overlap(coef, Data)
  plotp_Train <- separate2GroupsCox(as.vector(coef_test[[1]]), x_hat[, coef_test[[2]]], 
                                    plot.title = 'GSE1456_smfs', as.data.frame(y), 
                                    legend.outside = T)

  plot_train <- plotp_Train$plot
  plot_train
  
  ## for Xtile
  p_index <- cbind(y,plotp_Train$index)
  colnames(p_index) <- c(colnames(y), "riskscore")

  
  name <- dir[i]
  
  
  ## 按./拆分,提取
  name1 <- str_split_fixed(name, "./", 2);
  name2 <- str_split_fixed(name1[2], "/", 2)
  name3 <- str_split_fixed(name2[2], "[.]", 2);
  name4 <- name3[1]
  name5 <- str_c(name4,"_ex")
  
  # out <- rbind(p_train, CI_train)
  # colnames(out) <- name4
  
  
  setwd("D:\\E\\博士\\R_程序\\BRCA\\Feature_data")
  # path <- paste("./P_CI/",paste(name4,".txt", sep=""), sep="")
  # write.table(out, path, quote = F)
  # path <- paste("./Sangerbox/",paste(name4,".csv", sep=""), sep="")
  # write.csv(sanger, path, row.names = F)
  path <- paste("./Xtile/",paste(name4,".txt", sep=""), sep="")
  write.table(p_index, path, quote = F, row.names = F, sep="\t")
  path <- paste("./Figure/", paste(name5,".pdf", sep=""), sep="")
  pdf(path, width = 4.5, height = 4.5, onefile = FALSE)
  print(plot_train)
  plot_train
  dev.off()

  
  
  print("***第i个文件***")
  print(i)

  
}




