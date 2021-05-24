

# cor --------------------------------------------------------------
data  <- read.table("TCGA_BRCA_clin_1142_1080.txt", header = T, check.names = FALSE)
data1 <- as.matrix(data[-c(1,2),])
exp <- cbind(row.names(data1),data1)
dim(exp)    # 1142 1081

ex_v <- apply(as.matrix(exp[,-1]),2,
              function(x)as.numeric(x))
colnames(ex_v) <- colnames(exp)[-1]
rownames(ex_v) <- exp[,1]
matPCC <- cor(t(ex_v),method='pearson')
colnames(matPCC) <- exp[,1]
rownames(matPCC) <- exp[,1]
diag(matPCC) <- 0
used_matPCC <- matPCC
used_matPCC[upper.tri(matPCC)] <- 0

# set threshold  
thre <- 0.7
se_net <- rbind() 
for(i in 1:dim(used_matPCC)[1])
{
  for(j in 1:dim(used_matPCC)[2])
  {
    if(abs(used_matPCC[i,j]) > thre)
    {
      se_net <- rbind(se_net,c(colnames(used_matPCC)[i],
                               colnames(used_matPCC)[j],used_matPCC[i,j] ) ) 
    }
  }
}
dim(se_net)    # 520

#### cor network
library(igraph)
nnet <- graph.data.frame(se_net[,c(1,2)],directed = F) 
# plot(nnet, vertex.color="pink",vertex.size=8,
     # label.font=2,label.cex=2,label.color='black',main='candidate_net')

#### network add
net7 <- rbind(net1, se_net[,c(1,2)])
dim(net7)    # 6402    2
net8 <- net7[!duplicated(net7),]
dim(net8)    # 6402    2
# write.csv(net8, file = "Regnetwork_cor.csv",row.names = F)


# MI --------------------------------------------------------------------
refN <- as.matrix(read.csv("net_in_inter_genes.csv",header = T))
library(igraph)
nnet <- graph.data.frame(refN,directed = F) 
sim_nnet <- simplify(nnet,remove.multiple = T,remove.loops = T)
ed <- get.edgelist(sim_nnet)

data  <- read.table("TCGA_BRCA_clin_1142_1080.txt", header = T, check.names = FALSE)
data1 <- as.matrix(data[-c(1,2),])
exp <- cbind(row.names(data1),data1)
dim(exp)    # 1142 1081


ex_v <- apply(as.matrix(exp[,-1]),2,
              function(x)as.numeric(x))
colnames(ex_v) <- colnames(exp)[-1]
rownames(ex_v) <- exp[,1]
re <- ex_v


library(infotheo)
nbins <- sqrt(NROW(t(re)))
dis_data <- discretize(t(re), disc="equalfreq", nbins) 
colnames(dis_data) <- rownames(ex_v)
ex_v[1,1]
rownames(ex_v)
###### 0-MI
MI_0 <- rbind()
for(i in 1:dim(ed)[1])
{
  #i <- 1
  loc1 <- which(rownames(ex_v) == ed[i,1])
  loc2 <- which(rownames(ex_v) == ed[i,2])
  MI <- mutinformation(dis_data[,loc1],dis_data[,loc2],method="emp")
  MI_0 <- rbind(MI_0,MI)
}
median(MI_0)
mean(MI_0)

# hist(MI_0,main='MI_0_distribution')

MI_0_thre <- 0.55
MI_0_left <- rbind()
for(i in 1:length(MI_0))
{
  if(MI_0[i] >= MI_0_thre)
  {
    MI_0_left <- rbind(MI_0_left, c(ed[i,],MI_0[i]))
  }
}

dim(MI_0_left)    # 404   3
MI_0_net <- graph.data.frame(MI_0_left[,c(1,2)],directed = F)

#### network add
net9 <- rbind(net8, MI_0_left[,c(1,2)])
dim(net9)    # 7922    2
net10 <- net9[!duplicated(net9),]
dim(net10)    # 6402    2
# write.csv(net10, file = "Regnetwork_cor_mi.csv",row.names = F)



# 72gene is in the network ---------------------------------------------------------
node <- as.matrix(read.csv("union_gene.csv",header = T))  
node[1]
node1 <- as.matrix(c("TP53", "E2F1", "E2F2","JUN", "SP1", "MYC",
                     "EP300", "SPI1", "NKX6-1", "GLI1", "ESR1",
                     "HSF1", "STAT1", "XBP1", "SOX9" ))
node2 <- rbind(node, node1)
node_used <- node2


net_used <- as.matrix(read.csv('Regnetwork_cor_mi.csv', header=T))
k1 <- which(net_used[,1] %in% node_used)   # 2221
k2 <- which(net_used[,2] %in% node_used)   # 389
length(intersect(k1,k2))    # 190
used <- net_used[intersect(k1,k2),]
library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP, remove.loops = T, remove.multiple = T)  
ed <- as_edgelist(p1, names = TRUE)


## save
node0 <- cbind(node1, as.matrix(rep(1, 15)))
colnames(node0) <- c("add", "lab")
write.csv(node0,"corcmi_72_node.csv",row.names = F,quote = F)
write.csv(ed,"corcmi_72_net.csv",row.names = F,quote = F)



