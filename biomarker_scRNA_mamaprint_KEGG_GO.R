
Data = read.table("TCGA_pro_norm_clin.txt", header = T, check.names = FALSE)
dim(Data)    # 20222  1080
gene <- as.vector(rownames(Data)[-c(1,2)])


DE_outcome <- read.csv("DEG_Tumor_vs_Normal_10.csv", header = T, sep=',')    # 489
DE_clin <- read.csv("DEG_Dead_vs_Alive_p001.csv", header = T, sep=',')    # 501
Degene <- as.vector(union(DE_outcome[,1], DE_clin[,1]))    # 956


intersect(gene, Degene)
sub <- setdiff(gene, Degene)    # From the overall data, delete DEgene and


# biomarker ---------------------------------------------------------------
mark <- as.matrix(read.csv("128 biomarkers.csv", header = T, sep=','))
scmark <- as.matrix(c("KRT15","UBE2C","TOP2A","KRT6B","MKI67","HMGB2","ASPM","CDC20","KIF20A","CDK"))
mama_70 <- as.matrix(read.csv("70_mama.csv", header = T, sep=',') )
KEGG_147 <- as.matrix(read.csv("KEGG147.csv", header = T, sep=',') )


brcaterm <- as.matrix(read.csv('bc7.csv',header = T))
goterm <- as.matrix(read.csv('GO_annotations-9606-inferred-allev.csv',header = T))
term <- merge(goterm, brcaterm, by.x="go_id",by.y = "go_id_brca",all=FALSE) 
gene_symbol <- cbind(as.matrix(term$go_id), as.matrix(term$gene_symbols))
library(tidyverse)
gene <- str_split_fixed(gene_symbol[,2], "[|]", 453)
genet <- t(gene)
GO <- as.matrix(Reduce(union, genet[,1: dim(brcaterm)[1] ]))   #  取并 1128  取交无
GO_3104 <- as.matrix(GO[-which(GO[,1] == ""),])
dim(GO_3104)    # 519
colnames(GO_3104) <- c("gene_go")    


## intersect GO_3104 
intersect(mark, scmark)    
intersect(mark, mama_70)      
intersect(mark, KEGG_147)       
intersect(mark, GO_3104)       
intersect(scmark, mama_70)     
intersect(scmark, KEGG_147)    
intersect(scmark, GO_3104)     
intersect(mama_70, KEGG_147)    
intersect(mama_70, GO_3104)     
intersect(KEGG_147, GO_3104)    


# gene Integration ------------------------------------------------------------
add_gene <- as.matrix( union(union(union(union(intersect(sub, mark), 
                                               intersect(sub, scmark)), 
                                         intersect(sub, mama_70)), 
                                   intersect(sub, KEGG_147)), 
                             intersect(sub, GO_3104)) )    # 677
DEgene1 <- as.matrix(Degene)    # 956
allgene <- as.matrix(rbind(add_gene, DEgene1))    #1633

rownames(allgene) <- allgene[,1] 
data <- Data[rownames(allgene),]
dim(data)    # 1633 1080

all_data <- rbind(Data[c(1,2),], data)
dim(all_data)    # 2568 1080

## 1089 gene
gene1089 <- as.matrix(rownames(all_data)[-c(1,2)]) 
colnames(gene1089) <- c("gene")


# Extract gene ------------------------------------------------------------------
inter_gene <- row.names(all_data)[-c(1,2)]
# write.csv(inter_gene, file = "inter_gene.csv", row.names = F)


# network match -----------------------------------------------------------
net <- as.matrix(read.csv("renewed_Regnetwork_10_5.csv",header = T))
net[1,]

# node <- as.matrix(read.csv("inter_gene.csv",header = T))#[-1,]
node <- as.matrix(inter_gene)#[-1,]
node[1]

node_used <- node
net_used <- net[,c(2,4)]
k1 <- which(net_used[,1] %in% node_used)   # 48593
k2 <- which(net_used[,2] %in% node_used)   # 18360
length(intersect(k1,k2))    # 5229
# length(union(k1,k2))
used <- net_used[intersect(k1,k2),]
# used <- net_used[union(k1,k2),]


library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP, remove.loops = T, remove.multiple = T)  # 最终的数对
ed <- as_edgelist(p1, names = TRUE)    # 4734
# write.csv(ed,"net_in_inter_genes.csv",row.names = F,quote = F)


# Calculate the largest (weak or strong) connected component of the graph ---------------------------------------------------------
clu <- components(p1)
groups(clu)
a <- groups(clu)
a$`1`
a$`2`

g <- p1
pdf(file = "net_in_genes_adjp.pdf",width = 15,height = 15)
plot(g, layout=layout.fruchterman.reingold,  
     vertex.size=4,  
     vertex.label = V(g)$name,  
     vertex.label.cex=0.7,  
     vertex.label.dist=0.4, 
     vertex.label.color = "black"  
)
dev.off()


# The number of nodes in the network -----------------------------------------------------------
kk1 <- which( node_used %in% used[,1])   # 217
kk2 <- which( node_used %in% used[,2])   # 1406
length(union(kk1,kk2))    # 1143
# feature_new <- as.matrix(node_used[union(kk1,kk2),])   
feature_new <- as.matrix(a$`1`)
write.csv(feature_new,"gene_in_net.csv",row.names = F,quote = F)


# Extract the data of gene in the network ---------------------------------------------------------
net_gene <- read.csv("gene_in_net.csv", header = T, row.names=1)  
net_data <- Data[rownames(net_gene),]
dim(net_data)    # 1142 1080
net_clin_data <- rbind(Data[c(1,2),], net_data)
dim(net_clin_data)    # 1283 1080
write.table(net_clin_data, file = "TCGA_BRCA_clin_1142_1080.txt",quote = F, sep = "\t")


