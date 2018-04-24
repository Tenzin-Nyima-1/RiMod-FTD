# Date: 26 May, 2017;
# Project: RIMOD-FTD;
# Aim: to perform hypergeometric test to identify over-represented GO-BP-MF-CC in the DE(hits) from the MAPT vs controls test on frontal samples;
# NOTE: the annotation of the count-table was carried using RIKEN's CAT gtf file as the gene model;
# NOTE 2: the same script could be adapted for the count table from: all 7 brain regions; and, temporal brain; 
# Workflow (main steps);
# Step 1: create gene-set collections;
# Step 2: load the workspace that performed the DEG-test;
# Step 3: prepare the universe gene list;
# Step 4: prepare the hit gene list;
# Step 5: perform the hypergeometric test;
# Step 6: extract the significant processes;
# Step 4: save this workspace;


# load the relevant libraries;
options(stringsAsFactors=FALSE)
library(edgeR)
library(statmod)
library(gdata)
library(stringr)
library(gdata)
library(GO.db)
library(AnnotationDbi)
library(Biobase)
library(org.Hs.eg.db)
library(reshape2)
library(KEGG.db)
library(HTSanalyzeR)
library(snow)
library(muStat)



# ########################################;
# I. create a list of gene-set collections;
# ########################################;
GO_BP <- GOGeneSets(species="Hs", ontologies=c("BP"))
GO_MF <- GOGeneSets(species="Hs", ontologies=c("MF"))
GO_CC <- GOGeneSets(species="Hs", ontologies=c("CC"))
PW_KEGG <- KeggGeneSets(species="Hs")
ListGSC <- list(GO_BP=GO_BP, GO_MF=GO_MF, GO_CC=GO_CC, PW_KEGG=PW_KEGG)


# #########################################################################################################################################################################################;
# II. load the workspace that performed DGE analysis between MAPT vs. control samples #####################################################################################################;
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc_fro/S7_FronSamps_MAPTvsCont_RikAnn_SVA_DEG.RData")



# ###################################;
# III. Prepare the universe gene list;
# ###################################;

# 1. Extract all the expressed genes;
Fron_MAPTvsCont_SVA_toptags_All <- topTags(Fron_MAPTvsCont_SVA_glmQlFtest, n=Inf, adjust.method = "BH", sort.by = "PValue", p.value = 1)
dim(Fron_MAPTvsCont_SVA_toptags_All)
# [1] 23422    28;

# 2. select only the protein coding ones;
Fron_MAPTvsCont_SVA_toptags_All_mRNA <- Fron_MAPTvsCont_SVA_toptags_All$table[grep("coding_mRNA", Fron_MAPTvsCont_SVA_toptags_All$table$CAT_geneClass, ignore.case = T),]
dim(Fron_MAPTvsCont_SVA_toptags_All_mRNA)
# [1] 20630    28;
identical(grep("coding_mRNA", Fron_MAPTvsCont_SVA_toptags_All$table$CAT_geneClass, ignore.case = T), which(Fron_MAPTvsCont_SVA_toptags_All$table$CAT_geneClass == "coding_mRNA"))
# [1] TRUE;

# 3. now, select only the promoters;
Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt <- Fron_MAPTvsCont_SVA_toptags_All_mRNA[grep("Promoter", Fron_MAPTvsCont_SVA_toptags_All_mRNA$annotation, ignore.case = T),]
dim(Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt)
# [1] 20385    28;

# @Just a check;
unique(Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt$CAT_geneClass)
# [1] "coding_mRNA"
unique(Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt$annotation)
# [1] "Promoter (<=1kb)" "Promoter (2-3kb)" "Promoter (1-2kb)"

# @housekeeping;
which.na(Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt$entrez_ID)
# integer(0);
unique(Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt$entrez_ID[grep("na", Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt$entrez_ID)])
# [1] "__na";
length(Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt$entrez_ID[grep("na", Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt$entrez_ID)])
# [1] 326;

# 4. Cleaning: since there are na entrez, have to remove them;
Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt_cln <- Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt[-(which(Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt$entrez_ID == "__na")),]
dim(Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt_cln)
# [1] 20059    28;
identical(which(Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt$entrez_ID == "__na"), grep("na", Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt$entrez_ID))
# [1] TRUE;

# 5. Finally, select the universe genelist with gene-ids as entrez and values as logfc;
Fron_MAPTvsCont_Univ <- Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt_cln$logFC
names(Fron_MAPTvsCont_Univ) <- Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt_cln$entrez_ID
length(Fron_MAPTvsCont_Univ)
# [1] 20059;

# @crosschecks;
identical(as.numeric(Fron_MAPTvsCont_Univ), as.numeric(Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt_cln$logFC))
# [1] TRUE;
identical(as.character(names(Fron_MAPTvsCont_Univ)), as.character(Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt_cln$entrez_ID))
# [1] TRUE;
grep("na", names(Fron_MAPTvsCont_Univ), ignore.case = T)
# integer(0);



# ##############################;
# IV. Prepare the hits gene list;
# ##############################;

# 1. extract the hits i.e, < 0.05 FDR adjusted pvalue;
Fron_MAPTvsCont_Hits_df <- Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt_cln[which(Fron_MAPTvsCont_SVA_toptags_All_mRNA_promt_cln$FDR < 0.05),]
dim(Fron_MAPTvsCont_Hits_df)
# [1] 4710   28;

# @Crosscheck: if the above hit is similar to the DEG list from the main analysis;
identical(Fron_MAPTvsCont_Hits_df, Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod_promt[-(grep("__na", Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod_promt$entrez_ID, ignore.case = T)),])
# [1] TRUE; 
# NOTE: after excluding the na entrez, both tables are similar;

# 2. Now, save the hit entrez into a vector;
Fron_MAPTvsCont_Hits <- Fron_MAPTvsCont_Hits_df$entrez_ID
length(Fron_MAPTvsCont_Hits)
# [1] 4710;

# @just a check;
grep("na", Fron_MAPTvsCont_Hits, ignore.case = T)
# integer(0);



# ######################;
# V. Hypergeometric test;
# ######################;

# 1. Initialize and pre-process;
Fron_MAPTvsCont_gsca <- new("GSCA", listOfGeneSetCollections=ListGSC, geneList=Fron_MAPTvsCont_Univ, hits=Fron_MAPTvsCont_Hits)
# NOTE: tried generating gsca object with universe as just the entrez-ids, however error report: "'geneList' should be a named numeric or integer vector with length > 0!";
# -moreover this value could be used for GSEA;
Fron_MAPTvsCont_gsca_pp <- preprocess(Fron_MAPTvsCont_gsca, species="Hs", initialIDs="Entrez.gene", keepMultipleMappings=TRUE, duplicateRemoverMethod="max", orderAbsValue=FALSE)

# 2. now, perform the analysis;
Fron_MAPTvsCont_gsca_pp_ana <- analyze(Fron_MAPTvsCont_gsca_pp, para=list(pValueCutoff=0.05, pAdjustMethod="BH", minGeneSetSize=15), doGSOA=TRUE, doGSEA=FALSE)
Fron_MAPTvsCont_gsca_pp_ana_annot <- appendGSTerms(Fron_MAPTvsCont_gsca_pp_ana, keggGSCs="PW_KEGG", goGSCs=c("GO_BP", "GO_MF", "GO_CC"))

# 3. report the results in a html;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Frontal_RikAnno/HTSanalyzeR_HyperGeom_Results/")
report(Fron_MAPTvsCont_gsca_pp_ana_annot, experimentName="Fron_MAPTvsCont", species="Hs", allSig=TRUE, keggGSCs="PW_KEGG", goGSCs=c("GO_BP", "GO_MF", "GO_CC"), reportDir="Fron_MAPTvsCont_HypGeom")



# #########################;
# VI. subset the GO results;
# #########################;
summarize(Fron_MAPTvsCont_gsca_pp_ana_annot)

# A) GOBP results; 
Fron_MAPTvsCont_GOBP <- attributes(Fron_MAPTvsCont_gsca_pp_ana_annot)$result$HyperGeo.results$GO_BP
dim(Fron_MAPTvsCont_GOBP)
# [1] 1218    8;

# 1. hits;
Fron_MAPTvsCont_GOBP_Sig <- Fron_MAPTvsCont_GOBP[which(Fron_MAPTvsCont_GOBP$Adjusted.Pvalue < 0.05),]
dim(Fron_MAPTvsCont_GOBP_Sig)
# [1] 82  8;

# B) GOCC results; 
Fron_MAPTvsCont_GOCC <- attributes(Fron_MAPTvsCont_gsca_pp_ana_annot)$result$HyperGeo.results$GO_CC
dim(Fron_MAPTvsCont_GOCC)
# [1] 303   8;

# 1. hits;
Fron_MAPTvsCont_GOCC_Sig <- Fron_MAPTvsCont_GOCC[which(Fron_MAPTvsCont_GOCC$Adjusted.Pvalue < 0.05),]
dim(Fron_MAPTvsCont_GOCC_Sig)
# [1] 46  8;

# C) GOMF results; 
Fron_MAPTvsCont_GOMF <- attributes(Fron_MAPTvsCont_gsca_pp_ana_annot)$result$HyperGeo.results$GO_MF
dim(Fron_MAPTvsCont_GOMF)
# [1] 343   8;

# 1. hits;
Fron_MAPTvsCont_GOMF_Sig <- Fron_MAPTvsCont_GOMF[which(Fron_MAPTvsCont_GOMF$Adjusted.Pvalue < 0.05),]
dim(Fron_MAPTvsCont_GOMF_Sig)
# [1] 29  8;

# D) KEGG results; 
Fron_MAPTvsCont_KEGG <- attributes(Fron_MAPTvsCont_gsca_pp_ana_annot)$result$HyperGeo.results$PW_KEGG
dim(Fron_MAPTvsCont_KEGG)
# [1] 188   8;

# 1. hits;
Fron_MAPTvsCont_KEGG_Sig <- Fron_MAPTvsCont_KEGG[which(Fron_MAPTvsCont_KEGG$Adjusted.Pvalue < 0.05),]
dim(Fron_MAPTvsCont_KEGG_Sig)
# [1] 19  8;



# ########################;
# VII. Save this workspace;
# ########################;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc_fro/")
save.image("Fron_MAPTvsCont_DEG_HTSanalyzeR.RData")
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc_fro/Fron_MAPTvsCont_DEG_HTSanalyzeR.RData")
# END _______________________________________________________________________________________________________________________________________________________________________________;
