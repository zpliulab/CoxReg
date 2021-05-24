
data = read.table("TCGA_BRCA_clin_1142_1080.txt", header = T, check.names = FALSE)
datat <- as.data.frame(t(data))
Data <- as.matrix(data)

x <- model.matrix(status ~., datat)[,-c(1,2)]  
x_hat <- data.frame(x)
y <- as.matrix(datat[,c(1,2)])  
colnames(y) <- c("time", "status")
y_hat <- data.frame(y)

## coef
coef_gene <- read.csv("univariate_cox_coef.csv", header = T, sep=',')
coef <- as.matrix(coef_gene[,2])
rownames(coef) <- coef_gene[,1]

# my_overlap_coef ------------------------------------------------------
my_overlap <- function(x, y){

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
source('D:\\E\\博士\\R_程序\\BRCA\\R\\myoverlap_separate2GroupsCox.R')
library(ggpubr)
library(magrittr)
library(survminer)

coef_test <- my_overlap(coef, Data)
plotp_Train <- separate2GroupsCox(as.vector(coef_test[[1]]), x_hat[, coef_test[[2]]], 
                                  plot.title = 'GSE1456_smfs', as.data.frame(y), 
                                  legend.outside = T)
plot_train <- plotp_Train$plot



## for Xtile
p_index <- cbind(y,plotp_Train$index)
colnames(p_index) <- c(colnames(y), "riskscore")
# write.table(p_index, file = "TCGA_OS_3gene.txt", quote = F, row.names = F, sep="\t")
# write.table(p_index, file = "TCGA_OS_3gene.csv", quote = F, row.names = F)


# Extract features on TCGA --------------------------------------------------------------

library(survminer)
library(survival)
library(ggplot2)

setwd("D:\\E\\博士\\R_程序\\BRCA\\Data\\TCGA_NEW\\result")

data = read.table("TCGA_OS_3gene.txt", header = T, check.names = FALSE)
gene_name <- colnames(data)[3]
exprSet <- as.data.frame(t(data))

## Set cutoff value  
alpha <- 0.50
risk_score  <- t(as.matrix(exprSet[gene_name,]))
cut_off <- rep(as.numeric(quantile(exprSet[gene_name,],alpha)), dim(exprSet)[2])


data$time <- data$time/365
data$riskscore <- ifelse(risk_score > cut_off, 'high','low')
table(data$risk_score)

fit <- survfit(Surv(time, status)~riskscore, data = data)

p <- ggsurvplot(fit, data = data, 
                conf.int = F, 
                # surv.median.line = "hv", 
                risk.table = TRUE, 
                tables.height= 0.25, 
                cumcensor = T,   
                legend = c(0.83,0.95),
                
                # P value
                pval = TRUE, 
                pval.size=6, 
                font.pval= c(14, "bold", "black"),
                pval.coord = c(0.00, 0.05),
                
                # legend
                legend.title = '', # gene_name
                legend.labs=c("High risk", "Low risk"), 
                font.legend= c(14, "plain", "black"), 
                # font.main = c(100, "bold", "black"),
                # xlim = c(0,72), # present narrower X axis, but not affect
                # survival estimates.
                palette=c("red", "blue"),
                font.x = c(14, "plain", "black"),
                font.y = c(14, "plain", "black"), 
                font.tickslab = c(14, "plain", "black"),
                xlab = "Time in years", 
                break.time.by = 6
) # break X axis in time intervals by 500.



# 添加HR和CI

res_cox <- coxph(Surv(time, status) ~riskscore, data=data)
HR <- round(summary(res_cox)$conf.int[1],2)
ci_l <- round(summary(res_cox)$conf.int[3],2)
ci_r <- round(summary(res_cox)$conf.int[4],2)

p1 <- p
p1$plot = p1$plot + 
  ggplot2::annotate("text",x = 3.13, y = 0.3, label = paste("HR : ", HR), size = 5) + 
  ggplot2::annotate("text",x = 6.28, y = 0.2,
                    label = paste("(","95%CI : ", ci_l,"-",ci_r,")", sep = ""), size = 5)

pdf("TCGA_3_os.pdf", width = 4.8, height = 6, onefile = FALSE)  
p1
dev.off() 