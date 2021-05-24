library(dplyr)    
library(tidyr)
library(tidyverse) 


Data1 = read.table("TCGA_BRCA_clin_1142_1080_scale.txt", header = T, check.names = FALSE)
coef_gene <- read.csv("univariate_cox_coef.csv", header = T, sep=',')
gene <- as.matrix(coef_gene[,1])

colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))


genedata <- merge(gene, Data2, by = "gene")
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata1[1,1]
genedata2 <- rbind(Data1[c(1,2),],genedata1)    # 55
write.table(genedata2,"GSE_TCGA_3_ori.txt",quote=F,sep="\t")  # 2020.6.29



# Enter the PRS index obtained by TCGA ----------------------------------------------------------------
data = read.table("GSE_TCGA_3_ori.txt", header = T, check.names = FALSE)
datat <- as.data.frame(t(data))
prs = read.table("TCGA_OS_3gene.txt", header = T, check.names = FALSE)

datat$PRS <- prs$riskscore
data1 <- datat[,c(1,2,6,3,4,5)]
write.table(data1, file = "prs_TCGA_for_hiplot.txt", quote=F,sep="\t",row.names = F)  












