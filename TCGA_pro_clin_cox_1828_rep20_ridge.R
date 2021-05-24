library(caret) 
library(tidyverse)
library(survival)
library(ggpubr)
library(magrittr)
library(survminer)
library(glmSparseNet)
library(glmnet)
library(ncvreg)
library(APML0)
library(ncpen)



# Read data --------------------------------------------------------------------
data = read.table("TCGA_BRCA_clin_1142_1080.txt", header = T, check.names = FALSE)
mean(data[-c(1,2),1])


## scale
Data <- data.frame(cbind( t(data)[,c(1,2)], scale(t(data)[,-c(1,2)]) ))
write.table(t(Data), file = "TCGA_BRCA_clin_1142_1080_scale.txt",quote=F,sep="\t")

l <- 50
mean(Data[,l])
var(Data[,l])

## 将时间day转换为year
# Data$time <- Data$time/365 
# View(Data[,1:10])


## Data-training set + test set
set.seed(123)
training.samples <- Data$status %>% createDataPartition(p = 0.7, list = FALSE)
train.data  <- Data[training.samples, ]
test.data <- Data[-training.samples, ]

Data_train <- t(train.data)
x <- model.matrix(status ~., train.data)[,-c(1,2)]   
x_hat <- data.frame(x)
y <- as.matrix(train.data[,c(1,2)])  
colnames(y) <- c("time", "status")
y_hat <- data.frame(y)

Data_test <- t(test.data)
y1 <- t(Data_test[c(1,2),])
colnames(y1) <- c("time", "status")
y1_hat <- data.frame(y1)
x1 <- t(Data_test[-c(1,2),])
x1_hat <- data.frame(x1)


# funtion -----------------------------------------------------------------
# my_CI function ----------------------------------------------------------
my_CI <- function(feature_plus){
  univ_formulas <- sapply(feature_plus, 
                          function(x) as.formula(paste('Surv(y1_hat$time, y1_hat$status)~', x)))
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = x1_hat)})
  CI <- univ_models[[1]]$concordance[6]
  return(CI)
}

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
  gene <- rownames(y[-c(1,2),])
  overlap <- intersect(sele, gene)
  
  lab <- x[,1] %>% { .[. != 0]} %>% as.matrix
  coefs.v <- lab[overlap,]
  
  my <- list(coefs.v, overlap)
  return(my)
} 


######################## glmnet -- ridge ########################
Coef_ridge <- matrix()    
Ci_ridge <- matrix()   
P_ridge <- matrix()  

for(i in 1:20){

  set.seed(200*i)
  cv.ridge = cv.glmnet(x, y, family="cox", alpha = 0, nfolds = 10)
  lambda.min <- cv.ridge$lambda.min
  lambda.1se <- cv.ridge$lambda.1se
  print("***lambda.min、lambda.1se***")
  print(lambda.min)
  print(lambda.1se)
  plot(cv.ridge)

  ridge.model <- glmnet(x, y, family="cox", alpha = 0, lambda = cv.ridge$lambda.min)
  coef_ridge <- as.matrix(coef(ridge.model))
  which(coef_ridge[,1] != 0)
  feature_ridge <- colnames(x)[which(coef_ridge[,1] != 0)]
  # View(feature_ridge)
  ridge_plus <- paste(feature_ridge,collapse="+")


  CI_ridge <- my_CI(ridge_plus)
  CI_ridge

  coef_test <- my_overlap(coef_ridge, Data_test)
  plotp_ridge <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]],
                                    plot.title = 'Ridge', as.data.frame(y1), legend.outside = T)
  p_ridge <- plotp_ridge$pvalue

  tcross <- rep(i, length(coef_ridge))                 
  step_ridge <- data.frame(cbind(coef_ridge, tcross))
  Coef_ridge <- cbind(Coef_ridge, step_ridge)   


  kcross <- rep(i, length(CI_ridge))
  temp_ridge <- data.frame(cbind(CI_ridge, kcross))
  Ci_ridge <- cbind(Ci_ridge, temp_ridge)    


  pcross <- rep(i, length(p_ridge))
  pemp_ridge <- data.frame(cbind(p_ridge, pcross))
  P_ridge <- cbind(P_ridge, pemp_ridge)   

  print(paste(i))

}

write.csv(Coef_ridge, file = "coef_ridge.csv")
write.csv(P_ridge, file = "P_ridge.csv")
write.csv(Ci_ridge, file = "CI_ridge.csv")



# Set seed and reclassify data ------------------------------------------------------------
set.seed(12345)
training.samples <- Data$status %>% createDataPartition(p = 0.7, list = FALSE)
train.data  <- Data[training.samples, ]
test.data <- Data[-training.samples, ]


Data_train <- t(train.data)
x <- model.matrix(status ~., train.data)[,-c(1,2)]   
x_hat <- data.frame(x)
y <- as.matrix(train.data[,c(1,2)])  
colnames(y) <- c("time", "status")
y_hat <- data.frame(y)


Data_test <- t(test.data)
y1 <- t(Data_test[c(1,2),])
colnames(y1) <- c("time", "status")
y1_hat <- data.frame(y1)
x1 <- t(Data_test[-c(1,2),])
x1_hat <- data.frame(x1)



######################### glmnet -- Lasso  #############################
Coef_lasso <- matrix()      
Ci_lasso <- matrix()  
P_lasso <- matrix()

for(i in 1:20){
  
  set.seed(200*i)
  
  cv.lasso = cv.glmnet(x, y, family="cox", alpha = 1, nfolds = 10)
  lambda.min <- cv.lasso$lambda.min
  lambda.1se <- cv.lasso$lambda.1se
  print("***lambda.min、lambda.1se***")
  print(lambda.min)
  print(lambda.1se)
  plot(cv.lasso)
  
  
  ## 拟合
  lasso.model <- glmnet(x, y, family="cox", alpha = 1, lambda = cv.lasso$lambda.min)
  coef_lasso <- as.matrix(coef(lasso.model))
  which(coef_lasso[,1] != 0)
  feature_lasso <- colnames(x)[which(coef_lasso[,1] != 0)]
  lasso_plus <- paste(feature_lasso,collapse="+")
  
  
  CI_lasso <- my_CI(lasso_plus)
  CI_lasso
  
  coef_test <- my_overlap(coef_lasso, Data_test)
  plotp_lasso <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]], 
                                    plot.title = 'Lasso', as.data.frame(y1), legend.outside = T)
  p_lasso <- plotp_lasso$pvalue
  
  tcross <- rep(i, length(coef_lasso))               
  step_lasso <- data.frame(cbind(coef_lasso, tcross))
  Coef_lasso <- cbind(Coef_lasso, step_lasso)  
  
  
  kcross <- rep(i, length(CI_lasso)) 
  temp_lasso <- data.frame(cbind(CI_lasso, kcross))
  Ci_lasso <- cbind(Ci_lasso, temp_lasso)  
  
  
  pcross <- rep(i, length(p_lasso)) 
  pemp_lasso <- data.frame(cbind(p_lasso, pcross))
  P_lasso <- cbind(P_lasso, pemp_lasso)  
  
  print(paste(i)) 

}
write.csv(Coef_lasso, file = "coef_lasso.csv")
write.csv(P_lasso, file = "P_lasso.csv")
write.csv(Ci_lasso, file = "CI_lasso.csv")



#################################### elastic net -------------------------------------------------------------

Coef_elastic <- matrix()       
Ci_elastic <- matrix()   
P_elastic <- matrix() 


for(i in 1:20){
  
  set.seed(200*i)
  
  cv.elastic = cv.glmnet(x, y, family="cox", alpha = 0.5, nfolds = 10)
  lambda.min <- cv.elastic$lambda.min
  lambda.1se <- cv.elastic$lambda.1se
  print("***lambda.min、lambda.1se***")
  print(lambda.min)
  print(lambda.1se)
  plot(cv.elastic)
  
  
  ## 拟合
  elastic.model <- glmnet(x, y, family="cox", alpha = 0.5, lambda = cv.elastic$lambda.min)
  coef_elastic <- as.matrix(coef(elastic.model))
  which(coef_elastic[,1] != 0)
  feature_elastic <- colnames(x)[which(coef_elastic[,1] != 0)]
  elastic_plus <- paste(feature_elastic,collapse="+")
  
  
  CI_elastic <- my_CI(elastic_plus)
  CI_elastic
  
  coef_test <- my_overlap(coef_elastic, Data_test)
  plotp_elastic <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]], 
                                    plot.title = 'Elastic', as.data.frame(y1), legend.outside = T)
  p_elastic <- plotp_elastic$pvalue
  
  tcross <- rep(i, length(coef_elastic))               
  step_elastic <- data.frame(cbind(coef_elastic, tcross))
  Coef_elastic <- cbind(Coef_elastic, step_elastic)   
  
  
  kcross <- rep(i, length(CI_elastic)) 
  temp_elastic <- data.frame(cbind(CI_elastic, kcross))
  Ci_elastic <- cbind(Ci_elastic, temp_elastic)   

  pcross <- rep(i, length(p_elastic)) 
  pemp_elastic <- data.frame(cbind(p_elastic, pcross))
  P_elastic <- cbind(P_elastic, pemp_elastic)  
  
  print(paste(i)) 
  
}

write.csv(Coef_elastic, file = "coef_elastic.csv")
write.csv(P_elastic, file = "P_elastic.csv")
write.csv(Ci_elastic, file = "CI_elastic.csv")



######################### ncvreg -- SCAD 惩罚 #############################

X <- x
y <- y


Coef_SCAD <- matrix()      
Ci_SCAD <- matrix()   
P_SCAD <- matrix()


for(i in 1:20){
  
  set.seed(200*i)
  
  cv.SCAD <- cv.ncvsurv(X, y, family ="cox", penalty="SCAD")  
  lambda.min <- cv.SCAD$lambda.min 
  plot(cv.SCAD)
  print("*** lambda.min ***")
  print(lambda.min)
  
  
  SCAD.model <- ncvsurv(X, y, family ="cox", lambda = cv.SCAD$lambda.min, penalty="SCAD") 
  coef_SCAD <- as.matrix(coef(SCAD.model))
  which(coef_SCAD[,1] != 0)
  feature_SCAD <- colnames(x)[which(coef_SCAD[,1] != 0)]
  SCAD_plus <- paste(feature_SCAD,collapse="+")
  
  
  CI_SCAD <- my_CI(SCAD_plus)
  CI_SCAD
  
  coef_test <- my_overlap(coef_SCAD, Data_test)
  plotp_SCAD <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]], 
                                    plot.title = 'SCAD', as.data.frame(y1), legend.outside = T)
  p_SCAD <- plotp_SCAD$pvalue 

  tcross <- rep(i, length(coef_SCAD))               
  step_SCAD <- data.frame(cbind(coef_SCAD, tcross))
  Coef_SCAD <- cbind(Coef_SCAD, step_SCAD)   
  
  
  kcross <- rep(i, length(CI_SCAD)) 
  temp_SCAD <- data.frame(cbind(CI_SCAD, kcross))
  Ci_SCAD <- cbind(Ci_SCAD, temp_SCAD)  
  
  
  pcross <- rep(i, length(p_SCAD)) 
  pemp_SCAD <- data.frame(cbind(p_SCAD, pcross))
  P_SCAD <- cbind(P_SCAD, pemp_SCAD)   
  
  print(paste(i)) 


}  

write.csv(Coef_SCAD, file = "coef_SCAD.csv")
write.csv(P_SCAD, file = "P_SCAD.csv")
write.csv(Ci_SCAD, file = "CI_SCAD.csv")



######################### ncvreg -- MCP #############################

X <- x
y <- y


Coef_MCP <- matrix()     
Ci_MCP <- matrix()   
P_MCP <- matrix()

# set.seed(123456)

for(i in 1:20){
  
  set.seed(31*i)
  
  cv.MCP <- cv.ncvsurv(X, y, family ="cox", penalty="MCP")  
  lambda.min <- cv.MCP$lambda.min 
  print("*** lambda.min ***")
  print(lambda.min)
  plot(cv.MCP)
  
  MCP.model <- ncvsurv(X, y, family ="cox", lambda = cv.MCP$lambda.min, penalty="MCP") 
  coef_MCP <- as.matrix(coef(MCP.model))
  which(coef_MCP[,1] != 0)
  feature_MCP <- colnames(x)[which(coef_MCP[,1] != 0)]
  MCP_plus <- paste(feature_MCP,collapse="+")
  if(MCP_plus == ""){
    set.seed(31)
    cv.MCP <- cv.ncvsurv(X, y, family ="cox", penalty="MCP")  
    lambda.min <- cv.MCP$lambda.min 
    print("*** lambda.min ***")
    print(lambda.min)
    

    MCP.model <- ncvsurv(X, y, family ="cox", lambda = cv.MCP$lambda.min, penalty="MCP") 
    coef_MCP <- as.matrix(coef(MCP.model))
    which(coef_MCP[,1] != 0)
    feature_MCP <- colnames(x)[which(coef_MCP[,1] != 0)]
    MCP_plus <- paste(feature_MCP,collapse="+")
  }

    
  CI_MCP <- my_CI(MCP_plus)
  CI_MCP
  
  coef_test <- my_overlap(coef_MCP, Data_test)
  plotp_MCP <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]], 
                                    plot.title = 'Ridge', as.data.frame(y1), legend.outside = T)
  p_MCP <- plotp_MCP$pvalue
  p_MCP
  
  tcross <- rep(i, length(coef_MCP))             
  step_MCP <- data.frame(cbind(coef_MCP, tcross))
  Coef_MCP <- cbind(Coef_MCP, step_MCP)   
  
  
  kcross <- rep(i, length(CI_MCP)) 
  temp_MCP <- data.frame(cbind(CI_MCP, kcross))
  Ci_MCP <- cbind(Ci_MCP, temp_MCP)  

  pcross <- rep(i, length(p_MCP)) 
  pemp_MCP <- data.frame(cbind(p_MCP, pcross))
  P_MCP <- cbind(P_MCP, pemp_MCP)   
  
  
  print(paste(i)) 
  
}

write.csv(Coef_MCP, file = "coef_MCP.csv")
write.csv(P_MCP, file = "P_MCP.csv")
write.csv(Ci_MCP, file = "CI_MCP.csv")



############################ APML0 L0 ###################################


y <- Data[,c(1,2)]
colnames(y) <- c("time", "status")
x <- as.matrix(Data[,-c(1,2)])
dim(x)    # 1080  687

library(APML0)
set.seed(123456)
l0.model = APML0(x, y, family="cox", penalty="Lasso", nfolds=10) # Lasso
print(l0.model)

# lambda.min <- l0.model$lambda.min
# print("*** lambda.min ***")
# print(lambda.min)

# coef_l0 <- as.matrix(l0.model$Beta)
# rownames(coef_l0) <- rownames(coef_lasso)

lambda.opt <- l0.model$lambda.opt
print("*** lambda.opt ***")
print(lambda.opt)

coef_l0 <- as.matrix(l0.model$Beta0)
rownames(coef_l0) <- rownames(data)[-c(1,2)]

which(coef_l0[,1] != 0)
feature_l0 <- colnames(x)[which(coef_l0[,1] != 0)]
l0_plus <- paste(feature_l0,collapse="+")


CI_l0 <- my_CI(l0_plus)
CI_l0


coef_test <- my_overlap(coef_l0, Data_test)
plotp_l0 <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]], 
                               plot.title = 'L0', as.data.frame(y1), legend.outside = T)

P_l0 <- plotp_l0$pvalue
P_l0


write.csv(coef_l0, file = "coef_l0.csv")
write.csv(feature_l0, file = "feature_l0.csv", row.names = F)
write.csv(CI_l0, file = "CI_l0.csv")
write.csv(P_l0, file = "P_l0.csv")
pdf(file = "l0.pdf",width = 4,height = 3)
plotp_l0$plot
dev.off()



######################### ncpen  Bridge 1/2 #############################


x <- model.matrix(status ~., train.data)[,-c(1,2)]   # 删除了time, statuss~.
x_hat <- data.frame(x)
y <- as.matrix(train.data[,c(1,2)])  
colnames(y) <- c("time", "status")
y_hat <- data.frame(y)


Data_test <- t(test.data)
y1 <- t(Data_test[c(1,2),])
colnames(y1) <- c("time", "status")
y1_hat <- data.frame(y1)
x1 <- t(Data_test[-c(1,2),])
x1_hat <- data.frame(x1)


x_1 <- as.matrix(cbind(x,y[,2]))
y_1 <- as.matrix(y[,1])


library(ncpen)

set.seed(12345)
Bridge.model <- ncpen(y.vec = y_1, x.mat = x_1, family="cox", penalty="mbridge")
opt.lambda <- gic.ncpen(Bridge.model, pch="*", type="b")$opt.lambda
print("*** optional lambda ***")
print(opt.lambda)

n <- which(Bridge.model$lambda == opt.lambda)


coef <- coef(Bridge.model)
coef_Bridge <- as.matrix(coef[, n])   

which(coef_Bridge[,1] != 0)
feature_Bridge <- colnames(x)[which(coef_Bridge[,1] != 0)]
Bridge_plus <- paste(feature_Bridge,collapse="+")


CI_Bridge <- my_CI(Bridge_plus)
CI_Bridge


coef_test <- my_overlap(coef_Bridge, Data_test)
plotp_Bridge <- separate2GroupsCox(as.vector(coef_test[[1]]), x1_hat[, coef_test[[2]]], 
                                   plot.title = 'Bridge', as.data.frame(y1), legend.outside = T)
P_Bridge <- plotp_Bridge$pvalue
P_Bridge

                          
write.csv(coef, file = "coef_Bridge.csv")
write.csv(coef_Bridge, file = "coef_Bridge.csv")
write.csv(feature_Bridge, file = "feature_Bridge.csv", row.names = F)
write.csv(CI_Bridge, file = "CI_Bridge.csv")
write.csv(P_Bridge, file = "P_Bridge.csv")
pdf(file = "Bridge.pdf",width = 4,height = 3)
plotp_Bridge$plot
dev.off()



# rm(list = ls())

######################################## ###########################################################

ridge <- read.table("coef_ridge.csv", header=TRUE, sep = ',')
lasso <- read.table("coef_lasso.csv", header=TRUE, sep = ',')
Elastic <- read.table("coef_elastic.csv", header=TRUE, sep = ',')
Bridge <- read.table("coef_Bridge.csv", header=TRUE, sep = ',')
l0 <- read.table("coef_l0.csv", header=TRUE, sep = ',')
SCAD <- read.table("coef_SCAD.csv", header=TRUE, sep = ',')
MCP <- read.table("coef_MCP.csv", header=TRUE, sep = ',')



# my_cbind 20coef -----------------------------------------------------
my_cbind <- function(x){
  x1 <- matrix()
  x1 <- as.character(x[,1])
  for (i in 0:19){
    x1 <- cbind(x1, x[, 3+2*i])
  }
  return(x1)
}

coef_ridge1 <- my_cbind(ridge)
coef_lasso1 <- my_cbind(lasso)
coef_Elastic1 <- my_cbind(Elastic) 
# coef_Bridge1 <- my_cbind(Coef_Bridge)
# coef_l01 <- my_cbind(Coef_l0)  
coef_SCAD1 <- my_cbind(SCAD) 
coef_MCP1 <- my_cbind(MCP)  



# my_union_coef 20 non-0 --------------------------------------------------
my_union_coef <- function(x){
  x1 <- matrix(data=NA)
  j <- 1
  for (i in 1:dim(x)[1]){
    # i <- 2
    if (x[i,2] !=0 | x[i,3] !=0 | x[i,4] !=0 | x[i,5] !=0
        | x[i,6] !=0 | x[i,7] !=0 | x[i,8] !=0 | x[i,9] !=0
        | x[i,10] !=0 | x[i,11] !=0 | x[i,12] !=0 | x[i,13] !=0
        | x[i,14] !=0 | x[i,15] !=0 | x[i,16] !=0 | x[i,17] !=0
        | x[i,18] !=0 | x[i,19] !=0 | x[i,20] !=0 | x[i,21] !=0) {
      x1[j] <- as.character(x[i,1])
      j <- j+1
    }else{
      print("Wrong!")
    }
  }
  return(as.matrix(x1))
}


# Selected gene -----------------------------------------------------------
coef_ridge2 <- my_union_coef(coef_ridge1)
coef_lasso2 <- my_union_coef(coef_lasso1)
coef_Elastic2 <- my_union_coef(coef_Elastic1)
coef_Bridge2 <- as.matrix(Bridge[which(Bridge[,2] != 0),1]) 
coef_l02 <- as.matrix(l0[which(l0[,2] != 0),1])
coef_SCAD2 <- my_union_coef(coef_SCAD1)
coef_MCP2 <- my_union_coef(coef_MCP1)



## 读出 gene 数据
write.csv(coef_ridge2, "feature_ridge.csv",row.names = F)
write.csv(coef_lasso2, "feature_lasso.csv",row.names = F)
write.csv(coef_Elastic2, "feature_Elastic.csv",row.names = F)
# write.csv(coef_Elastic2, "feature_Elastic62.csv")  
write.csv(coef_l02, "feature_l0.csv",row.names = F)
write.csv(coef_Bridge2, "feature_Bridge.csv",row.names = F)
write.csv(coef_SCAD2, "feature_SCAD.csv",row.names = F)
write.csv(coef_MCP2, "feature_MCP.csv",row.names = F)


# 计算个数 --------------------------------------------------------------------
sum(coef_ridge2[,1] != 0)
sum(coef_lasso2[,1] != 0)
sum(coef_Elastic2[,1] != 0)
sum(coef_l02[,1] != 0)
sum(coef_Bridge2[,1] != 0)
sum(coef_SCAD2[,1] != 0)
sum(coef_MCP2[,1] != 0)


# 计算CI --------------------------------------------------------------------
Ci_ridge <- read.table("CI_ridge.csv", header=TRUE, sep = ',')
Ci_lasso <- read.table("CI_lasso.csv", header=TRUE, sep = ',')
Ci_Elastic <- read.table("CI_elastic.csv", header=TRUE, sep = ',')
Ci_Bridge <- read.table("CI_Bridge.csv", header=TRUE, sep = ',')
Ci_l0 <- read.table("CI_l0.csv", header=TRUE, sep = ',')
Ci_SCAD <- read.table("CI_SCAD.csv", header=TRUE, sep = ',')
Ci_MCP <- read.table("CI_MCP.csv", header=TRUE, sep = ',')

# my_ci 提取20次的ci -----------------------------------------------------

my_ci <- function(x){
  x1 <- matrix()
  for (i in 1:20){
    x1 <- cbind(x1, x[, 2*i+1])
  }
  return(x1)
}

Ci_ridge1 <- as.matrix(my_ci(Ci_ridge))[,-1]
Ci_lasso1 <- as.matrix(my_ci(Ci_lasso))[,-1]
Ci_Elastic1 <- as.matrix(my_ci(Ci_Elastic))[,-1]
Ci_l01 <- Ci_l0[,-1]
Ci_Bridge1 <- Ci_Bridge[,-1]
Ci_SCAD1 <- as.matrix(my_ci(Ci_SCAD))[,-1]
Ci_MCP1 <- as.matrix(my_ci(Ci_MCP))[,-1]


mean(Ci_ridge1)
mean(Ci_lasso1)
mean(Ci_Elastic1)
Ci_l01
Ci_Bridge1
mean(Ci_SCAD1)
mean(Ci_MCP1)

# 计算 P ---------------------------------------------------------------------
P_ridge <- read.table("P_ridge.csv", header=TRUE, sep = ',')
P_lasso <- read.table("P_lasso.csv", header=TRUE, sep = ',')
P_Elastic <- read.table("P_elastic.csv", header=TRUE, sep = ',')
P_Bridge <- read.table("P_Bridge.csv", header=TRUE, sep = ',')
P_l0 <- read.table("P_l0.csv", header=TRUE, sep = ',')
P_SCAD <- read.table("P_SCAD.csv", header=TRUE, sep = ',')
P_MCP <- read.table("P_MCP.csv", header=TRUE, sep = ',')



P_ridge1 <- as.matrix(my_ci(P_ridge))[,-1]
P_lasso1 <- as.matrix(my_ci(P_lasso))[,-1]
P_Elastic1 <- as.matrix(my_ci(P_Elastic))[,-1]
P_l01 <- P_l0[,-1]
P_Bridge1 <- P_Bridge[,-1]
P_SCAD1 <- as.matrix(my_ci(P_SCAD))[,-1]
P_MCP1 <- as.matrix(my_ci(P_MCP))[,-1]


mean(P_ridge1)
mean(P_lasso1)
mean(P_Elastic1)
P_l01
P_Bridge1
mean(P_SCAD1)
mean(P_MCP1)



# Integrate 7 features ---------------------------------------------------------------

ridge <- as.matrix(read.table("feature_ridge.csv", header=TRUE, sep = ','))
lasso <- as.matrix(read.table("feature_lasso.csv", header=TRUE, sep = ','))
Elastic <- as.matrix(read.table("feature_Elastic.csv", header=TRUE, sep = ','))
Bridge <- as.matrix(read.table("feature_Bridge.csv", header=TRUE, sep = ','))
l0 <- as.matrix(read.table("feature_l0.csv", header=TRUE, sep = ','))
SCAD <- as.matrix(read.table("feature_SCAD.csv", header=TRUE, sep = ','))
MCP <- as.matrix(read.table("feature_MCP.csv", header=TRUE, sep = ','))


## CI>0.6  lasso elastic L0 SCAD
uniondata <- union(lasso, union(Elastic, union(l0, SCAD)))
# write.csv(uniondata, file = "union_gene.csv", row.names = F)

a1 <- intersect(lasso, Elastic)
a2 <- intersect(lasso, l0)
a3 <- intersect(lasso, SCAD)
a4 <- intersect(Elastic, l0)
a5 <- intersect(Elastic, SCAD)
a6 <- intersect(l0, SCAD)

interdata <- union(a1, union(a2, union(a3, union(a4, union(a5, a6)))))

interdata <- intersect(lasso, intersect(Elastic, intersect(l0, SCAD)))
# uniondata <- union(Elastic, union(l0, SCAD))
# write.csv(interdata, file = "intersect_gene.csv", row.names = F)


