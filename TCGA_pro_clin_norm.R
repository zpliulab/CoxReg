
library(glmSparseNet)
library(curatedTCGAData)
library(TCGAutils)
library(dplyr)


# 设置路径 --------------------------------------------------------------------


# Extract expression matrix of TNBC subtypes of the BRCA dataset of the TCGA database ------------------------------------------

# Construct survival analysis model, pure tumor sample expression matrix, survival information of corresponding tumor patient

# params <- list(seed = 29221)  
brca <- curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm", FALSE)
head(getSubtypeMap(brca))
head(getClinicalNames("BRCA"))


brca.primary.solid.tumor <- TCGAutils::splitAssays(brca, '01')
xdata.raw <- t(assay(brca.primary.solid.tumor[[1]]))
dim(xdata.raw)    # 1093 20501


xdata.raw <- xdata.raw %>%
{ (apply(., 2, sd) != 0) } %>%
{ xdata.raw[, .] }
dim(xdata.raw)    # 1093 20220


# Get survival information
ydata.raw <- colData(brca.primary.solid.tumor) %>% as.data.frame %>%
  # Keep only data relative to survival or samples
  select(patientID, vital_status,
         Days.to.date.of.Death, Days.to.Date.of.Last.Contact,
         days_to_death,         days_to_last_followup,
         Vital.Status) %>%
  # Convert days to integer 
  mutate(Days.to.date.of.Death = as.integer(Days.to.date.of.Death)) %>%
  mutate(Days.to.Last.Contact  = as.integer(Days.to.Date.of.Last.Contact)) %>%
  # Find max time between all days (ignoring missings)
  rowwise %>%
  mutate(time = max(days_to_last_followup,        Days.to.date.of.Death,
                    Days.to.Last.Contact, days_to_death, na.rm = TRUE)) %>%
  # Keep only survival variables and codes
  select(patientID, status = vital_status, time) %>%
  # Discard individuals with survival time less or equal to 0
  filter(!is.na(time) & time > 0) %>% as.data.frame
dim(ydata.raw)    # 1080    3


# Set index as the patientID
rownames(ydata.raw) <- ydata.raw$patientID

# Get matches between survival and assay data
xdata.raw_1 <- xdata.raw[TCGAbarcode(rownames(xdata.raw)) %in%
                         rownames(ydata.raw),]
dim(xdata.raw_1)    # 1080 20501


# Order ydata the same as assay
ydata.raw    <- ydata.raw[TCGAbarcode(rownames(xdata.raw_1)), ]
xdata <- xdata.raw_1
ydata <- ydata.raw %>% select(status)
# ydata <- ydata.raw %>% select(time, status)

data <- as.matrix(cbind(ydata, xdata))
rownames(data) <- rownames(xdata) 
dim(data)    # 1080 20221
data1 <- t(data)
write.table(data1,"TCGA_pro_outcome.txt",quote=F,sep="\t") 


# Perform DEGs according to survival status and normalized -------------------------------------------
data_outcome <- read.table("TCGA_pro_outcome.txt",header=T,sep='\t', check.names = F)
dim(data_outcome)   # 20221  1080
data_outcome[2,1]

## Express data
xdatat <- data_outcome[-1,]
dim(xdatat)   # 20220  1080

## State of existence
ydata <- t(data_outcome[1,])
sum(ydata[,1])    # 152 dead

## Label grouping
ydata[which(ydata[,1] == "1"), 1] <- c("Dead")
ydata[which(ydata[,1] == "0"), 1] <- c("Alive")

# DESeq2 package to do DEGs of RNA-seq data -------------------------------------------------
exprSet <- round(xdatat)
dim(exprSet)    # 20220  1080
group_list <- as.factor(ydata)
colData <- data.frame(row.names=colnames(exprSet), group_list=group_list)
colnames(colData) <- c("outcome")
dim(colData)



library(DESeq2)
library(limma)
library(pasilla)
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ outcome)
dds2 <- DESeq(dds)  
resultsNames(dds2)


# res <-  results(dds2, contrast=c("outcome","Dead","Alive"))
# summary(res) 
# res_order <- res[order(res$padj),]
# res_order <- as.data.frame(res_order)
# write.csv(res_order,file= "DEG_res_order_DA.csv")


res1 <-  results(dds2, alpha = 0.01)
write.csv(res1,file= "DEG_res.csv")

diff_gene_deseq2 <- subset(res1, padj < 0.01)
# diff_gene_deseq2 <- subset(res1, padj < 0.05 & abs(log2FoldChange) > 0.5849625)
diff_gene_deseq2 <- as.data.frame(diff_gene_deseq2)
dim(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file= "DEG_Dead_vs_Alive_p001.csv")


# Data normalization -------------------------------------------------------------------
vst_dat <- vst(dds2, blind = TRUE)
dat111 <- assay(vst_dat)
dim(dat111)    # 20220  1080

data_norm_clin <- t(dat111)
dim(data_norm_clin)    # 1080 20220


# 带有生存时间的 -----------------------------------------------------------------
ydata_clin <- ydata.raw %>% select(time, status)
data_clin <- as.matrix(cbind(ydata_clin, data_norm_clin))
rownames(data_clin) <- rownames(data_norm_clin) 
dim(data_clin)    # 1080 20222
data_clin1 <- t(data_clin)
write.table(data_clin1,"TCGA_pro_norm_clin.txt",quote=F,sep="\t") 
