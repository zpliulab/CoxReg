
library(survival)
library(plyr)


data = read.table("TCGA_feature_data.txt", header = T, check.names = FALSE)
Data <- data.frame(t(data))
str(Data)


Data$status<-factor(Data$status)
summary(Data$status)



###################  Batch single factor regression -----------------------

#1.Build model y
y <- Surv(time=Data$time, event=Data$status==0)  #0为复发

#2.Batch single factor regression model establishment：Uni_cox_model
Uni_cox_model<- function(x){
  FML <- as.formula(paste0 ("y~",x))
  cox <- coxph(FML,data=Data)
  cox1 <-summary(cox)
  coef <- cox1$coefficients[,1]
  HR <- round(cox1$coefficients[,2],2)
  PValue <- round(cox1$coefficients[,5],3)
  CI5 <- round(cox1$conf.int[,3],2)
  CI95 <- round(cox1$conf.int[,4],2)
  Uni_cox_model<- data.frame('Characteristics' =x,
                             names <-rownames(cox1$conf.int),
                             'coef' = coef,
                             'HR' = HR,
                             'CI5' = CI5,
                             'CI95' = CI95,
                             'Uni_P' = PValue
  )
  return(Uni_cox_model)}  

#3.Input the desired single factor regression variables into the model

#3-(1)View the name and serial number of the variable
names(Data)

#3-(2)Enter the variable number
variable_names <- colnames(Data)[c(3:74)]  

#4.Output result
Uni_cox <- lapply(variable_names, Uni_cox_model)
Uni_cox <- ldply(Uni_cox, data.frame)

#5.Optimize the table, here is an example of HR+95% CI+P style
Uni_cox$HR.CI95 <- paste0(Uni_cox$HR," (",Uni_cox$CI5,'-',Uni_cox$CI95,")")
Uni_cox <- Uni_cox[,-4:-6] #HR (95% CI)+P 

# Uni_cox$CI<-paste(Uni_cox$CI5,'-',Uni_cox$CI95)
# Uni_cox<-Uni_cox[,-4:-6] #方法2：HR+95% CI+P 


##################### 二、Single factor regression p<0.05 do multiple factor regression ----------------------------

#1.Extract single factor p<0.05 variable
Uni_cox$Characteristics[Uni_cox$Uni_P<0.05]

#2.Multi-factor model establishment
mul_cox_model<- as.formula(paste0 ("y~",
                                   paste0(Uni_cox$Characteristics[Uni_cox$Uni_P<0.05],
                                          collapse = "+")))
mul_cox<-coxph(mul_cox_model,data=Data)
cox4 <- summary(mul_cox) 
coef <- cox4$coefficients[,c(1,5)]
mul_coef <- coef[which(coef[,2] <= 0.05),1]
# write.csv(mul_coef, file = "univariate_cox_coef.csv")

#3.Extract the information of multi-factor regression
mul_corf <- cox4$coefficients[,1] 
mul_HR <- round(cox4$coefficients[,2],2) 
mul_PValue<- round(cox4$coefficients[,5],4) 
mul_CI1<-round(cox4$conf.int[,3],2)
mul_CI2<-round(cox4$conf.int[,4],2)

#4.Multi-factor results are optimized and combined into a table: mul_cox1
##4-1：HR(95%CI)+P 
mul_HR.CI95 <- paste(mul_HR,"(",mul_CI1,'-',mul_CI2,")")
mul_cox1 <- data.frame("coef"=mul_corf, "mul_HR.CI95"=mul_HR.CI95,"P"=mul_PValue)

##4-2：HR+95%CI+P 
#mul_CI<-paste(mul_CI1,'-',mul_CI2)
#mul_cox1<- data.frame("HR"=mul_HR,"mul_CI"=mul_CI, "P"=mul_PValue)


#####################  三、Single factor and multiple factors are integrated into one table  ------------------------------
#1.Delete the first column of the uni-factor table Uni_cox
Uni_cox <- Uni_cox[,-1]

#2.After deleting the first column, the first column of the uni_cox table Uni_cox is named Characteristics
colnames(Uni_cox)[1] <- 'Characteristics'

#3.The row name of the multi-factor table mul_cox1 is placed in the cell, also named Characteristics
mul_cox1 <- cbind(rownames(mul_cox1), mul_cox1, row.names=NULL); names(mul_cox1 )[1]<-"Characteristics"

#4.The Uni_cox table and the mul_cox1 table are merged with the Characteristic column as the standard
table2 <- merge.data.frame(Uni_cox, mul_cox1, by="Characteristics", all = T, sort = T)

#5.Final form
table3 <- table2[,c(1,2,4,3,5,6,7)]
colnames(table3) <- c("Feature", "Coef_Uni", "HR(95%CI)_Uni", "P_Uni",
                     "Coef_Mul", "HR(95%CI)_Mul", "P_Mul")
# write.csv(table3, file = "univariate_cox.csv",row.names = F)


# latex ----------------------------------------------------------------
library(stargazer)
stargazer(table3, summary=FALSE, rownames=FALSE) 

