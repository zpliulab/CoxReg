
gene <- as.matrix(read.csv("union_gene.csv",header = T))
colnames(gene) <- c('gene')


library(org.Hs.eg.db)
library(clusterProfiler)


# genes symbol ----------------------------------------------------------
genelist <- as.character(gene)
eg <- bitr(genelist, fromType="SYMBOL", toType=c("ENTREZID","GENENAME"), OrgDb="org.Hs.eg.db"); 
head(eg)


# go ----------------------------------------------------------------------
geneList <- eg$ENTREZID
go_BP <- enrichGO(gene = geneList,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
head(go_BP)[,1:6]
class(go_BP)
write.csv(go_BP@result, file = "cluster_GO_NEW.csv", row.names = F)

# Perform simple visualization ----------------------------------------------------------------
barplot(go_BP, showCategory=10,drop=T)    # showCategory=10,title="Enrichment GO"
dotplot(go_BP, showCategory=10,title="Enrichment GO Top10") #ÅİÅİÍ¼
# plotGOgraph(go_BP) 	
library(ggnewscale)
cnetplot(go_BP)

aa <- go_BP@result[1:30,]
bb <- aa[-c(9,17,20),]

cnetplot(go_BP,
        showCategory = 5,
        foldChange = NULL,
        layout = "kk",
        colorEdge = T,
        circular = T,
        node_label = "all"
)
  
  
enrichGO = DOSE::setReadable(go_BP, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
enrichGO
library(ggnewscale)
cnetplot(enrichGO)


library(enrichplot)
library(DOSE)
# data(geneList)
# de <- names(geneList)[1:100]
# x <- enrichDO(de)
# x2 <- pairwise_termsim(x)

x2 <- go_BP
cnetplot(x2)
# use `layout` to change the layout of map
cnetplot(x2, layout = "star")
# use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
cnetplot(x2, showCategory = 10)
categorys <- c("cell fate commitment", "canonical Wnt signaling pathway",
               "cardiac chamber development", "mitral valve morphogenesis",
               "Notch signaling pathway", "regulation of Wnt signaling pathway")

plot <-  cnetplot(go_BP,
         showCategory = categorys,
         foldChange = NULL,
         layout = "kk",    # kk, gem
         colorEdge = T,
         circular = F,
         node_label = "gene",
         cex_category = 1.0,
         cex_gene = 1.0,
         cex_label_category = 1.0,
         cex_label_gene = 0.6,
         shadowtext = "none")

pdf(file = "cnetplot_cluster_NEW.pdf", width = 6.5, height = 4.5)
plot
dev.off()
## ¡¯star¡¯, ¡¯circle¡¯, ¡¯gem¡¯, ¡¯dh¡¯, ¡¯graphopt¡¯, ¡¯grid¡¯, ¡¯mds¡¯, ¡¯randomly¡¯, ¡¯fr¡¯, ¡¯kk¡¯, ¡¯drl¡¯, ¡¯lgl¡¯


