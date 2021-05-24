library(dplyr)     
library(tidyr)
library(tidyverse)    


Data1 = read.table("TCGA_pro_outcome_TN_log.txt", header = T, check.names = FALSE)
coef_gene <- read.csv("univariate_cox_coef.csv", header = T, sep=',')
gene <- as.matrix(coef_gene[,1])
dim(Data1)    # 20223   224

colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))


genedata <- merge(gene, Data2, by = "gene")
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata1[1,1]
genedata1 <- apply(genedata1, 2, as.numeric)


## scale
genedata1 <-  t(scale(t(genedata1)))
rownames(genedata1) <- genedata[,1]
genedata1[3,3]


genedata2 <- rbind(Data1[1,],genedata1)
# write.table(genedata2,"GSE_TCGA_3_ori.txt",quote=F,sep="\t")


# PRS ----------------------------------------------------------------
phi_coef <- read.csv("univariate_cox_coef.csv", header = T, sep=',')
genedata3 <- read.table("GSE_TCGA_3_ori.txt", header = T, check.names = FALSE)


colnames(phi_coef) <- c('gene',"coef")
Data3 <- cbind(rownames(genedata3[-1,]), genedata3[-1,])
colnames(Data3) <- c('gene', colnames(genedata3))
genedata <- merge(phi_coef, Data3, by = "gene")
genedata[2,3]

## Calculation
A <- t(as.matrix(genedata[,3:dim(genedata)[2]]))
A <- apply(A, 2, as.numeric)
dim(A)    # 121  67
A[2,1]
b <- as.matrix(genedata[,2])
dim(b)
b[2]
phi <- A %*% b
row.names(phi) <- colnames(genedata3)
phi[1,1]

phi1 <- as.data.frame(cbind(phi, t(genedata3[1,])))
colnames(phi1) <- c("PHI", "Label")
phi1$PHI <- as.numeric(phi1$PHI)
phi1[1,1]
write.table(phi1, file = "phi_GSE_TCGA_3.txt", quote=F,sep="\t")  


# Box plot --------------------------------------------------------------------
pred_log0 = read.table("phi_GSE_TCGA_3.txt", header = T, check.names = FALSE)
pred_log  <- rbind(pred_log0[113:224,], pred_log0[1:112,])
dim(pred_log)


library(ggplot2)
library(ggpubr)
library(Rcpp)

my_comparisons <- list( c("Tumor", "Normal") )
plot <- ggboxplot(pred_log, x = "Label", y = "PHI",
                  color = "Label", palette = "npg"#,  "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty"
) +  # add = "jitter" , fill = c("#90CAF9", "#F48FB1")
  stat_boxplot(geom = "errorbar", width=0.30, size=0.6) +   
  geom_boxplot(size=0.6, fill=c("#099052", "#EF3122")) +  # #00468B", "#ED0000
  stat_compare_means(comparisons=my_comparisons, method = "t.test") + # Add p-value
  theme(legend.position="none",
        axis.text.x=element_text(colour="black",size=14), 
        axis.text.y=element_text(size=14,face="plain"),
        axis.title.y=element_text(size = 14,face="plain"),
        axis.title.x=element_text(size = 14,face="plain"),
        plot.title = element_text(size=15,face="bold",hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  ylab("Risk score")+xlab("") 

pdf(file = "Box_GSE_TCGA_new_color.pdf",width = 3,height = 4)
plot
dev.off()