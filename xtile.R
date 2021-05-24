library(survminer)
library(survival)
library(ggplot2)


myfile = list.files("Xtile")               
dir = paste("./Xtile/", myfile, sep="")  
n = length(dir) 


# Take cut-off one by one and draw pictures one by one ---------------------------------------------------------
i <- 1

data = read.table(file = dir[i],header=T, check.names = FALSE) 
dir[i]
gene_name <- colnames(data)[3]
exprSet <- data.frame(t(data))

## Set cutoff value
alpha <- 0.500   # Obtained by Xtile

risk_score  <- t(as.matrix(exprSet[gene_name,]))
cut_off <- rep(as.numeric(quantile(exprSet[gene_name,],alpha)), dim(exprSet)[2])


# data$time <- data$time/12
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
                xlab = "Time in years", # customize X axis label. year
                break.time.by = 2
               ) # break X axis in time intervals by 500.


# Ìí¼ÓHRºÍCI 
res_cox <- coxph(Surv(time, status) ~riskscore, data=data)
HR <- round(summary(res_cox)$conf.int[1],2)
ci_l <- round(summary(res_cox)$conf.int[3],2)
ci_r <- round(summary(res_cox)$conf.int[4],2)

p1 <- p
p1$plot = p1$plot + 
          ggplot2::annotate("text",x = 2.13, y = 0.3, label = paste("HR : ", HR), size = 5) + 
          ggplot2::annotate("text",x = 3.88, y = 0.2,
                            label = paste("(","95%CI : ", ci_l,"-",ci_r,")", sep = ""), size = 5)


library(tidyverse)

name <- dir[i]
name1 <- str_split_fixed(name, "./", 2);
name2 <- str_split_fixed(name1[2], "/", 2)
name3 <- str_split_fixed(name2[2], "[.]", 2);
name4 <- name3[1]
name5 <- str_c(name4,"_ex")

path <- paste("./Xtilefigure/", paste(name5,".pdf", sep=""), sep="")
pdf(path, width = 4.8, height = 6, onefile = FALSE)  
p1
dev.off() 


