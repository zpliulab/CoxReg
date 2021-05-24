elastic <- as.matrix(read.table("union_gene.csv", header=TRUE, sep = ','))
pvalue <- as.matrix(read.csv("DEG_res_order_TN.csv", header=TRUE, sep = ','))

gene_pvalue <- as.matrix(pvalue[which(pvalue[,1] %in% elastic[,1]),c(1,7)])
colnames(gene_pvalue) <- c("72gene", "padj")

write.csv(gene_pvalue, file = "union_gene_pvalue.csv",row.names = F, quote = F)


