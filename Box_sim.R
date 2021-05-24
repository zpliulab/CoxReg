library(GOSim)
library(GOSemSim)


setOntology(ont = "BP", loadIC=TRUE, DIR=NULL)
H <-getAncestors()
HH <-getOffsprings()
HHH <-unlist(HH)
d <- godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)


b <-as.matrix(read.csv("cluster_GO_NEW.csv",header = T)[1:60,1]) 
a <- b
## selected 23 go terms
a <- b[c(1,2,3,12,13,14,15,17,19,20,22,23,24,25,26,28,30,34,35,37,38,46,48),]

aa <-as.matrix(read.csv("bc shishi.csv"))[,2]
amlsim <- as.matrix(mgoSim(a,aa,semData=d,measure = "Wang",combine=NULL))
# write.csv(amlsim,"consimi.csv", row.names = F)


# my go terms --------------------------------------------------------------------
library(corrplot)
corrplot(amlsim)  

dasimr <- as.matrix(apply(amlsim,1,max))
dasimc <- as.matrix(apply(amlsim,2,max))
congoterm <-as.matrix(rbind(dasimr, dasimc))

l <- length(dasimr)
n <- length(dasimr) + length(dasimc)

type1 <- as.matrix(rep("Selected",n))
congoterm1 <- data.frame(cbind(congoterm,type1))
colnames(congoterm1) <- c("sim","type")
congoterm1$sim <-as.numeric(as.vector(congoterm1$sim))
congoterm1hou <-as.matrix(congoterm1)

dasimrmean <-mean(dasimr)
dasimcmean <-mean(dasimc)


# random go terms --------------------------------------------------------------
Boxdata  <- matrix() 
 
for(i in 1:10){
  
  set.seed(4689*i)
  rand <- as.matrix(sample(HHH, size=l, replace = FALSE))
  row.names(rand) <-NULL
  randsim <-as.matrix(mgoSim(rand,aa,semData=d,measure = "Wang",combine=NULL))
  
  randsimr <- as.matrix(apply(randsim,1,max))
  randsimc <- as.matrix(apply(randsim,2,max))
  randomterm <- as.matrix(rbind(randsimr,randsimc))
  type2 <- as.matrix(rep("Random",n))
  randomterm1 <- data.frame(cbind(randomterm,type2))
  randomterm1$X1 <- as.numeric(as.vector(randomterm1$X1))
  randomterm1hou <- as.matrix(randomterm1)
  
  
  # Consistent go term and random integration --------------------------------------------------------
  boxdata <- data.frame(rbind(congoterm1hou,randomterm1hou))
  boxdata$sim <- as.numeric(as.vector(boxdata$sim))
  
  randsimrmean <- mean(randsimr)
  randsimcmean <- mean(randsimc)
  
  
  tcross <- rep(i, length(boxdata))               
  step <- data.frame(cbind(boxdata, tcross))
  Boxdata <- cbind(Boxdata, step)                  

  print(paste(i)) 
  
}


# my_cbind 20 coef -----------------------------------------------------
my_cbind <- function(x){
  x <- Boxdata
  x1 <- matrix()
  x1 <- as.character(x[,3])
  for (i in 0:9){
    x1 <- as.matrix(cbind(x1, x[, 2+3*i]))
  }
  return(x1)
}

coef <- my_cbind(Boxdata)
coef2 <- apply(coef[,-1], 2, as.numeric)
class(boxdata)
## Averaging each row
boxdata1 <- as.data.frame(cbind(rowMeans(coef2), boxdata[,2]))
colnames(boxdata1) <- c("SS-measure", "Type")
rownames(boxdata1) <- rownames(boxdata)
boxdata1$`SS-measure` <- as.numeric(as.vector(boxdata1$`SS-measure`))
# write.csv(boxdata1, "boxdata_cluster.csv", row.names = F)

class(boxdata[2,1])
class(boxdata1[2,1])

# Box plot --------------------------------------------------------------------
library(ggpubr)
my_comparisons <- list( c("Selected", "Random") )
p <- ggboxplot(boxdata1, x = "Type", y = "SS-measure",
          color = "Type", palette = "lancet", # "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty"
          add = "jitter", fill = c("#90CAF9", "#F48FB1"))+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "#C5E1A5")# Add global p-value

pdf(file = "TCGA_NEW\\result\\box_sim_cluster_NEW.pdf",width = 4.5,height = 4.5)
p
dev.off()

# select cluster GO ---------------------------------------------------------------
pr <- 0.7
which(dasimr[,1] >= pr)
bb0 <- as.matrix(read.csv("cluster_GO_NEW.csv",header = T)[,c(1,2)]) 
aa0 <- as.matrix(bb0[c(1,2,3,12,13,14,15,17,19,20,22,23,24,25,26,28,30,34,35,37,38,46,48),])


goselect <- aa0[which(dasimr[,1] >= pr),]
goselect[,2]
# write.csv(goselect, "go_6_select_NEW.csv", row.names = F)

