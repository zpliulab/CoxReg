library(BiocGenerics)
library(parallel)
library(Biobase)
library(dplyr)       
library(tidyr)
library(tidyverse)    
library(fdrtool)      
library(data.table)   
library(fdrtool)      
library(data.table)   


exprSet <- read.table("GSE42568_series_matrix.txt",header=T,sep='\t',fill=TRUE,strip.white = T)
dim(exprSet)    # 54675   122


anno <- read.table("GPL570.txt",header=T,sep='\t',fill=TRUE,strip.white = T,quote = "")  


library(dplyr)
library(tidyr)
anno2 <- anno %>%    
  select(ID,ENTREZ_GENE_ID) %>%           
  filter(ENTREZ_GENE_ID != '') 
colnames(anno2) <- c('ID_REF','EntrezID')   
anno2$ID_REF <- as.character(anno2$ID_REF) 


exprset2 <- exprSet %>%                  
  inner_join(anno2,by='ID_REF') %>%       
  select(ID_REF,EntrezID, everything())    


exprset3 <- exprset2
a <- tibble(exprset3[,1:2])


test1 <- apply(a,1, function(x){
  str_split(x[2],'///',simplify=T)      
} )


test2 <- apply(a, 1, function(x){         
  paste(x[1],str_split(x[2],'///', simplify=T), sep = "---")
})


unlist(test2)                             
x <- tibble(unlist(test2))                
colnames(x) <- "lala"                    


x2 <- separate(x,lala,c("id","entrezID"),sep = '---')   
x2[1:10,1]
exprset3[1:10,1]
x3 <- merge(x2,exprset3,by.x = "id",by.y="ID_REF",all=FALSE) 


x4<-x3[,-c(1,3)]                       
x4[1:6,1]
zz <- as.matrix(apply(as.matrix(x4[,1]),1,function(x) as.numeric(x)))


XX <- x4[,-1]
colnames(XX)[1:3]
XX1 <- cbind(zz,XX)
colnames(XX1) <- c("entrezID",colnames(XX))


homo<-read.table("homo.txt",header=T,sep='\t')
homo[1:6,1]
x5 <- merge(homo,XX1,by.x="GeneID",by.y = "entrezID",all=FALSE) 


expset4 <- x5 %>%
  dplyr::select(-GeneID) %>%             
  mutate(rowIQR = apply(.[,-1],1,IQR)) %>% 
  arrange(desc(rowIQR)) %>%              
  distinct(Symbol,.keep_all = T) %>%     
  dplyr::select(-rowIQR) %>%                 
  tibble::column_to_rownames(colnames(.)[1]) 
dim(expset4)   # 21835   121
class(expset4)
expset4[3,3]
# write.table(expset4,"GSE42568_expr.txt",quote=F,sep="\t")  

expset5 <- expset4

# label ---------------------------------------------------------------------

lable = read.csv("GSE42568_all_NT.csv", header = T, sep=',')
status <- as.matrix(lable$Tissue)
colnames(status) <- c("Type")


sample <- as.matrix(lable$Accession)
colnames(sample) <- c("Sample")
y <- cbind(sample,status)


data = rbind(as.matrix(t(y[,2])), as.matrix(expset5))
rownames(data) <- c("Type", rownames(expset5))
View(data[1:10,1:10])
dim(data)    # 21836   121
data[3,3]
expset5[2,3]
# write.table(data,"GSE42568_outcome_scale.txt",quote=F,sep="\t")  










