# Cox_Reg
Cox_Reg: A computational method to detect prognostic biomarkers of breast cancer based on gene expression data and regularized Cox proportional hazards models. 


1. TCGA_pro_clin_DE.R  --  First process the data to get the data of all samples. Then select 112 Tumor + 112 Normal samples to get DEGs.
2. TCGA_pro_clin_norm.R -- First process the data to get the data of all samples 1080. Then follow 152 dead + 928 alive to get DEGs.
3. biomarker_scRNA_mamaprint_KEGG_GO.R -- First integrate prior information (Biomarker 128, scRNA 10, Mamaprint 70, KEGG 147, GO terms). Then combine with RegNatwork to obtain genes with connected structures in the network.
4. TCGA_pro_clin_cox_1828_rep20_ridge.R -- Feature selection results of 7 methods, and finally union, 72.
5. network_match.R -- Observe the network structure of the union gene.
6. net_cor_mi.R -- Enter net_in_inter_genes, integrate cor, increase MI. Then add 15 genes to form a network with 72 genes.
7. gene_id_pvalue.R -- Use the P-value of DEG_res_order_TN.csv, two genes are NA.
8. feature_select_use.R -- Use the data TCGA_BRCA_clin_1142_1080_scale.txt to extract 72 feature genes of all samples.
9. Univariate_cox_for_72.R - Perform single-factor and multi-factor Cox regression for 72 features, save univariate_cox.csv. According to the multi-factor Cox, obtain the gene and coef of Risk Score.
10. TCGA_pro_clin_nomogram.R - Used to count the clinical information of TCGA, that is, Table 1.
11. feature_select_TCGA.R - Extract coefficients of 3 genes from TCGA_pro_outcome_TN_log.txt, compare between normal and tumor, get box plot.
