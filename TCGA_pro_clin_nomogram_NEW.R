library(glmSparseNet)
library(curatedTCGAData)
library(TCGAutils)
library(dplyr)


# params <- list(seed = 29221)  
brca <- curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm", FALSE)
dim(brca)

brca.primary.solid.tumor <- TCGAutils::splitAssays(brca, '01')
xdata.raw <- t(assay(brca.primary.solid.tumor[[1]]))
dim(xdata.raw)    # 1093 20501

# Get survival information
ydata.raw <- colData(brca.primary.solid.tumor) %>% as.data.frame %>%
  # Keep only data relative to survival or samples
  select(patientID, vital_status,
         Days.to.date.of.Death, Days.to.Date.of.Last.Contact,
         days_to_death,         days_to_last_followup,
         Vital.Status,    
         years_to_birth, pathologic_stage, pathology_T_stage, pathology_N_stage, pathology_M_stage) %>%
  # Convert days to integer  
  mutate(Days.to.date.of.Death = as.integer(Days.to.date.of.Death)) %>%
  mutate(Days.to.Last.Contact  = as.integer(Days.to.Date.of.Last.Contact)) %>%
  # Find max time between all days (ignoring missings)
  rowwise %>%
  mutate(time = max(days_to_last_followup,        Days.to.date.of.Death,
                    Days.to.Last.Contact, days_to_death, na.rm = TRUE)) %>%
  # Keep only survival variables and codes
  select(patientID, status = vital_status, time,     
         years_to_birth, pathologic_stage, pathology_T_stage, pathology_N_stage, pathology_M_stage) %>%
  # Discard individuals with survival time less or equal to 0
  filter(!is.na(time) & time > 0) %>% as.data.frame


# Set index as the patientID
rownames(ydata.raw) <- ydata.raw$patientID

# Get matches between survival and assay data
xdata.raw_1 <- xdata.raw[TCGAbarcode(rownames(xdata.raw)) %in%
                           rownames(ydata.raw),]
dim(xdata.raw_1)    # 1080 20501

# Order ydata the same as assay
ydata.raw    <- ydata.raw[TCGAbarcode(rownames(xdata.raw_1)), ]
xdata <- xdata.raw_1
ydata <- ydata.raw 


data <- as.matrix(cbind(ydata, xdata))
rownames(data) <- rownames(xdata) 
dim(data)    # 1080 20221


ydata.raw_2 <- data[,1:8]
ydata.raw_2[,1] <- rownames(ydata.raw_2)
write.csv(ydata.raw_2, file = "TCGA_clin_information.csv", row.names = F)


## Be sure to ensure that the samples of PRS correspond to the samples of clin one-to-one
clin <- read.csv("TCGA_clin_information.csv", header = T)


# pathologic_stage      8 NA
clin[which(clin[,5] == "stage i"), 5] <- c(0)   
clin[which(clin[,5] == "stage ia"), 5] <- c(0)  
clin[which(clin[,5] == "stage ib"), 5] <- c(0)
clin[which(clin[,5] == "stage ii"), 5] <- c(0)  
clin[which(clin[,5] == "stage iia"), 5] <- c(0)
clin[which(clin[,5] == "stage iib"), 5] <- c(0)  
clin[which(clin[,5] == "stage iii"), 5] <- c(1)
clin[which(clin[,5] == "stage iiia"), 5] <- c(1)
clin[which(clin[,5] == "stage iiib"), 5] <- c(1)
clin[which(clin[,5] == "stage iiic"), 5] <- c(1)
clin[which(clin[,5] == "stage iv"), 5] <- c(1)
clin[which(clin[,5] == "stage x"), 5] <- c(1)


# pathology_T_stage
clin[which(clin[,6] == "t1"), 6] <- c(0)   
clin[which(clin[,6] == "t1a"), 6] <- c(0)  
clin[which(clin[,6] == "t1b"), 6] <- c(0)
clin[which(clin[,6] == "t1c"), 6] <- c(0)  
clin[which(clin[,6] == "t2"), 6] <- c(0)
clin[which(clin[,6] == "t2a"), 6] <- c(0)
clin[which(clin[,6] == "t2b"), 6] <- c(0)
clin[which(clin[,6] == "t3"), 6] <- c(1)
clin[which(clin[,6] == "t3a"), 6] <- c(1)
clin[which(clin[,6] == "t4"), 6] <- c(1)
clin[which(clin[,6] == "t4b"), 6] <- c(1) 
clin[which(clin[,6] == "t4d"), 6] <- c(1)
clin[which(clin[,6] == "tx"), 6] <- c(1) 


# pathology_N_stage
clin[!duplicated(clin[,7]),7]
clin[which(clin[,7] == "n0"), 7] <- c(0)   
clin[which(clin[,7] == "n0 (i-)"), 7] <- c(0) 
clin[which(clin[,7] == "n0 (i+)"), 7] <- c(0)
clin[which(clin[,7] == "n0 (mol+)"), 7] <- c(0)
clin[which(clin[,7] == "n1"), 7] <- c(0) 
clin[which(clin[,7] == "n1a"), 7] <- c(0)
clin[which(clin[,7] == "n1b"), 7] <- c(0) 
clin[which(clin[,7] == "n1c"), 7] <- c(0) 
clin[which(clin[,7] == "n1mi"), 7] <- c(0)
clin[which(clin[,7] == "n2"), 7] <- c(1)
clin[which(clin[,7] == "n2a"), 7] <- c(1)
clin[which(clin[,7] == "n3"), 7] <- c(1) 
clin[which(clin[,7] == "n3a"), 7] <- c(1)
clin[which(clin[,7] == "n3b"), 7] <- c(1) 
clin[which(clin[,7] == "n3c"), 7] <- c(1)
clin[which(clin[,7] == "nx"), 7] <- c(1)


# pathology_M_stage
clin[!duplicated(clin[,8]),8]
clin[which(clin[,8] == "m0"), 8] <- c(0)
clin[which(clin[,8] == "cm0 (i+)"), 8] <- c(0)
clin[which(clin[,8] == "m1"), 8] <- c(1) 
clin[which(clin[,8] == "mx"), 8] <- c(1)
class(clin)


clinn <- as.matrix(clin[,-1])
prs <- read.table("TCGA_OS_3gene.txt",header = T, check.names = FALSE)
clin1 <- cbind(clinn[,1], clinn[,2], prs[,3], clinn[,3], clinn[,4], clinn[,5], 
               clinn[,6], clinn[,7])


colnames(clin1) <- c("status", "time", "PHI", "years_to_birth", "pathologic_stage", 
                     "pathology_T_stage", "pathology_N_stage", "pathology_M_stage")
clin <- as.data.frame(clin1)
clin[3,3]
class(clin)


# Change the number in the variable of data.frame type to num numeric type
clin_n <- as.data.frame(lapply(clin[,c(1,2,3,4,5,6,7,8)],as.numeric))
clin_n[3,3]
# write.csv(clin_n, file = "TCGA_clin_information_hitplot_NEW.csv", row.names = F)


library(survival)
## Error in summary.rms(f) : adjustment values not defined here or with datadist for years_to_birth
dd <- datadist(clin_n)
options(datadist="dd")


# Establish COX regression equation
f <- cph(Surv(time, status) ~ PHI+years_to_birth+pathologic_stage+pathology_T_stage+pathology_N_stage+pathology_M_stage,
         x=T, y=T, surv=T, data=clin_n, time.inc=365)


# Analyze the COX regression equation
summary(f)


# Draw a nomogram
surv <- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(365, x), function(x) surv(1095, x),   # 36-3year, 60-5year
                            function(x) surv(1825, x)), 
                lp=F, funlabel=c("1-year survival", "3-year survival", "5-year survival"),
                maxscale=100)
pdf(file = "nomogram_prs.pdf",width = 21,height = 8)
plot(nom)
dev.off()


# Draw calibration curve -------------------------------------------------------------------
library(rms) 
s <- Surv(clin_n$time, clin_n$status, type="right")

## 365 day
f <- cph(s~PHI+years_to_birth+pathologic_stage+
           pathology_T_stage+pathology_N_stage+pathology_M_stage, 
         x=TRUE, y=TRUE, surv = TRUE, time.inc=365, data=clin_n)
cal <- calibrate(f,u=365,cmethod='KM',m=180)   

pdf(file = "calibration_1_NEW.pdf",width = 4.5,height = 5.0)
plot(cal, xlim = c(0.88,1), ylim= c(0.88,1),   
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),
     col=c(rgb(255,0,0,maxColorValue=255)))
abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))
dev.off()


## 365*3   1095
f <- cph(s~PHI+years_to_birth+pathologic_stage+
           pathology_T_stage+pathology_N_stage+pathology_M_stage, 
         x=TRUE, y=TRUE, surv = TRUE, time.inc=1095, data=clin_n)
cal <- calibrate(f,u=1095,cmethod='KM',m=180)  

pdf(file = "calibration_3_NEW.pdf",width = 4.5,height = 5.0)
plot(cal, xlim = c(0.68,1), ylim= c(0.68,1),  
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),
     col=c(rgb(255,0,0,maxColorValue=255)))
abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255)))
dev.off()


## 365*3   1825
f <- cph(s~PHI+years_to_birth+pathologic_stage+
           pathology_T_stage+pathology_N_stage+pathology_M_stage, 
         x=TRUE, y=TRUE, surv = TRUE, time.inc=1825, data=clin_n)
cal <- calibrate(f,u=1825,cmethod='KM',m=180)  


pdf(file = "calibration_5_NEW.pdf",width = 4.5,height = 5.0)
plot(cal, xlim = c(0.4,1), ylim= c(0.4,1),  
     errbar.col=c(rgb(0,0,0,maxColorValue=255)),
     col=c(rgb(255,0,0,maxColorValue=255)))
abline(0,1,lty=3,lwd=2,col=c(rgb(0,0,255,maxColorValue= 255))) 
dev.off()
