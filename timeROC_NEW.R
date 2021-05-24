clin_prs = read.table("prs_TCGA_for_hiplot.txt", header = T, check.names = FALSE)
data_prs <- clin_prs[,c(1,2,3)]

######################################################################################
library(timeROC)

survivalROC_helper <- function(t) {  survivalROC(Stime = data_prs$time,
                                                 status = data_prs$status,
                                                 marker = data_prs$PRS,   
                                                 predict.time =t, method="KM")}
library(tidyverse)
library(survivalROC)


survivalROC_data <- data_frame(t = 365*c(1,3,5)) %>%  
  mutate(survivalROC = map(t, survivalROC_helper),  ## Extract scalar AUC  
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),  ## Put cut off dependent values in a data_frame  
         df_survivalROC = map(survivalROC, function(obj){  as_data_frame(obj[c("cut.values","TP","FP")])  })) %>%  
  dplyr::select(-survivalROC) %>%  
  unnest() %>%  
  arrange(t, FP, TP)


survivalROC_data1 <- survivalROC_data %>%   
  mutate(auc =sprintf("%.3f",auc))%>%   
  unite(year, t, auc, sep = " year AUC: ")



# Replace day in the legend with year -------------------------------------------------------
year <- str_split_fixed(as.matrix(survivalROC_data1[,1]), "[ ]", 3) 
year[which(year[,1] == "365"),1] <- c("1")
year[which(year[,1] == "1095"),1] <- c("3")
year[which(year[,1] == "1825"),1] <- c("5")

year1 <- str_c(year[,1], year[,2], year[,3], sep = " ")
survivalROC_data2 <- cbind(year1, survivalROC_data1[,c(2:4)])
survivalROC_data1 <- survivalROC_data2
################################################################################


AUC =factor(survivalROC_data1$year)
plot <- survivalROC_data1 %>%  
  ggplot(mapping = aes(x = FP, y = TP)) +  
  geom_path(aes(color= AUC))+  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+  
  theme_bw() +  
  theme(legend.position = c(0.7,0.2)) +
  theme(axis.text.x = element_text(size = 14,color="black"),
        axis.text.y = element_text(size = 14,color="black"),
        axis.title.x =element_text(size=14), 
        axis.title.y=element_text(size=14),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12))


pdf(file = "timeROC_NEW.pdf",width = 4.5,height = 4.5)
plot
dev.off()


# risk_factor_analysis ----------------------------------------------------
data_risk <- cbind(Data[,c(1,2)], prs, Data[,-c(1,2)])
dim(data_risk)
data_risk[3,3]
# write.csv(data_risk, file = "TCGA_feature_risk_hitplot.csv", row.names = F)

