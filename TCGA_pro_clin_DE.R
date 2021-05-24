library(glmSparseNet)
library(curatedTCGAData)
library(TCGAutils)
library(dplyr)



# Distinguish between tumor samples and normal tissue expression samples -------------------------------------------------------
brca <- curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm", FALSE)  
brca <- TCGAutils::splitAssays(brca, c('01','11'))
xdata.raw <- t(cbind(assay(brca[[1]]), assay(brca[[2]])))
dim(xdata.raw)    # 1205 20501


# Get matches between survival and assay data
class.v        <- TCGAbiospec(rownames(xdata.raw))$sample_definition %>% factor
names(class.v) <- rownames(xdata.raw)


xdata_raw <- xdata.raw %>%
{ (apply(., 2, sd) != 0) } %>%
{ xdata.raw[, .] }

dim(xdata_raw)    # 1205 20222
View(xdata_raw[,1:10])

small_subset <- colnames(xdata.raw)
dim(small.subset)    # 1205 20222

xdata <- xdata_raw[, small_subset[small_subset %in% colnames(xdata_raw)]]
dim(xdata)    # 1205 20222
View(xdata[,1:10])
xdatat <- t(xdata)

ydata <- class.v
View(ydata)
class(ydata)

ydata[which(ydata == "Primary Solid Tumor")] <- c("Tumor")
ydata[which(ydata == "Solid Tissue Normal")] <- c("Normal")
# save(xdata,ydata,file = 'TCGA_pro_outcome.rdata')


# Save classification data of all samples -------------------------------------------------------------
# xydata <- cbind(ydata, xdata)
# dim(xydata)    # 1205 20223
# xydata[which(xydata[,1] == "1"), 1] <- c("Tumor")
# xydata[which(xydata[,1] == "2"), 1] <- c("Normal")
# write.table(xydata, "TCGA_pro_outcome.txt",quote=F,sep="\t")



# 112 Normal samples are matched with Tumor samples of the same individual -----------------------------------------
which(ydata[,1] == "Primary Solid Tumor")


## Normal sample
sample_N <- rownames(xdata)[which(ydata == "Solid Tissue Normal")]
sample_N1 <- sample_N %>% 
  as_tibble() %>% 
  mutate(sample_N = substr(sample_N, 1, 12)) 

## All sample
sample <- rownames(xdata)
sample1 <- sample %>% 
  as_tibble() %>% 
  mutate(sample = substr(sample, 1, 12)) 

## Tumor sample
lab <- which(as.matrix(sample1[,2]) %in% as.matrix(sample_N1[,2]))
xdata_TN <- xdata[lab,]
dim(xdata_TN)    # 224 20222

ydata_TN <- rbind(as.matrix(rep(c("Tumor"), 112)), as.matrix(rep(c("Normal"), 112)))
colnames(ydata_TN) <- c("outcome")

data_TN <- cbind(ydata_TN, xdata_TN)
dim(data_TN)    # 224 20223
# write.table(t(data_TN), "TCGA_pro_outcome_TN.txt",quote=F,sep="\t")




## Differential expression analysis and normalization-------------------------------------
library(DESeq2)
library(limma)
library(pasilla)


data <- read.table("TCGA_pro_outcome_TN.txt",header=T,sep='\t', check.names = F)
dim(data)   # 20223   224
xdatat <- t(as.matrix(apply(as.matrix(data[-1,]),1,function(x) as.numeric(x))))
colnames(xdatat) <- colnames(data)
xdatat[1,2]
xdata[2,1]
ydata <- data[1,]


class(xdatat)
class(xdata)


# DESeq2 package to do differential analysis of RNA-seq data -------------------------------------------------
exprSet <- round(xdatat)
dim(exprSet)    # 20222  1205    20222   224

ydata_TN <- rbind(as.matrix(rep(c("Tumor"), 112)), as.matrix(rep(c("Normal"), 112)))
colnames(ydata_TN) <- c("outcome")

group_list <- as.factor(ydata_TN)
# group_list <- xydata[,1]
colData <- data.frame(row.names=colnames(exprSet), group_list=group_list)
colnames(colData) <- c("outcome")
dim(colData)

dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ outcome)
dds2 <- DESeq(dds)   
resultsNames(dds2)


# Extract the results of the DEGs, we are here to compare the Tumor group to the Normal group
res <-  results(dds2, contrast=c("outcome","Tumor","Normal"))
summary(res) 

res_order <- res[order(res$padj),]
res_order <- as.data.frame(res_order)
# write.csv(res_order,file= "DEG_res_order_TN.csv")


# Determine the threshold and screen for DEGs -----------------------------------------------------------
## FDR Correct
res1 <-  results(dds2, alpha = 0.01)  # Default FDR<0.1
# write.csv(res1,file= "DEG_res.csv")

diff_gene_deseq2 <- subset(res1, padj < 0.01 & abs(log2FoldChange) > 3.321928)
diff_gene_deseq2 <- as.data.frame(diff_gene_deseq2)
dim(diff_gene_deseq2)
# write.csv(diff_gene_deseq2,file= "DEG_Tumor_vs_Normal_10.csv")


# Data normalization -------------------------------------------------------------------
vst_dat <- vst(dds2, blind = TRUE)
dat111 <- assay(vst_dat)
dim(dat111)    # 20222  1205    20222   224

data0 <- rbind(t(ydata_TN), dat111)
dim(data0)    # 20223   224
data0[3,3]
# write.table(data0,"TCGA_pro_outcome_TN_log.txt",quote=F,sep="\t") 




