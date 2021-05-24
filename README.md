# Cox_Reg
Cox_Reg: A computational method to detect prognostic biomarkers of breast cancer based on gene expression data and regularized Cox proportional hazards models. 


1. TCGA_pro_clin_DE.R  --  First process the data to get the data of all samples. Then select 112 Tumor + 112 Normal samples to get DEGs.
2. TCGA_pro_clin_norm.R -- First process the data to get the data of all samples 1080. Then follow 152 dead + 928 alive to get DEGs.
3. biomarker_scRNA_mamaprint_KEGG_GO.R -- First integrate prior information (Biomarker 128, scRNA 10, Mamaprint 70, KEGG 147, GO terms). Then combine with RegNatwork to obtain genes with connected structures in the network.
