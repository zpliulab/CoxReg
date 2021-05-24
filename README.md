# Cox_Reg
Cox_Reg: A computational method to detect prognostic biomarkers of breast cancer based on gene expression data and regularized Cox proportional hazards models. 


1. TCGA_pro_clin_DE.R  --  First process the data to get the data of all samples. Then select 112 Tumor + 112 Normal samples to get DEGs.
2. TCGA_pro_clin_norm.R -- First process the data to get the data of all samples 1080. Then follow 152 dead + 928 alive to get DEGs.
3. biomarker_scRNA_mamaprint_KEGG_GO.R -- First integrate prior information (Biomarker 128, scRNA 10, Mamaprint 70, KEGG 147, GO terms). Then combine with RegNatwork to obtain genes with connected structures in the network.
4. TCGA_pro_clin_cox_1828_rep20_ridge.R -- Feature selection results of 7 methods, and finally union, 72.
5. network_match.R -- Observe the network structure of the union gene.
6. net_cor_mi.R - Enter net_in_inter_genes, integrate cor, increase MI. Then add 15 genes to form a network with 72 genes.
