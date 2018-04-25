# Date: 26-27 March 2018;
# Project: RIMOD-FTD;
# Aim: Weighted Gene Coexpression Network Analysis(WGCNA) of frontal samples: all cases and control samples;
# NOTE: the annotation of the count-table was carried using RIKEN's CAT gtf file as the gene model;
# NOTE 2: the same script could be adapted for the count table from: all 7 brain regions; and, temporal brain; 
# Workflow (main steps);
# Step 2: step-by-step network construction and module detection;


options(stringsAsFactors=FALSE)
library(AnnotationDbi)
library(edgeR)
library(stringi)
library(stringr)
library(muStat)
library(WGCNA)
library(edgeR)
library(impute)
library(GO.db)
library(preprocessCore)
enableWGCNAThreads()


# #################################################################;
# I. Upload the workspace that prepared the data for WGCNA analysis;
# _________________________________________________________________;
#load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc/WGCNA_S1_FronSamps_AllCasesNConts_CAGEseq.RData")

# 1. recall the final expression data;
dim(reqFro_alCasesVsCont_normd_logcpm_Trs)
# [1]    48 21051;
Fro_alCasesVsCont_nGenes
# [1] 21051;
Fro_alCasesVsCont_nSamples
# [1] 48;

# -@sample-id prepared from the above table;
# reqFro_AlCaseVSCont_SampIDs <- stri_split_fixed(stri_split_fixed(rownames(reqFro_alCasesVsCont_normd_logcpm_Trs), "sample_", simplify = T)[,2], "_no", simplify = T)[,1]
length(reqFro_AlCaseVSCont_SampIDs)
# [1] 48;

# 2. trait info ordered as per the order of samples in the expres set;
dim(reqFroAlCaseVsCont_TraitInfo)
# [1] 48 14; NOTE: contains cols with char class;

# -@all col classes converted into numeric classes;
dim(finFroAlCaseVsCont_TraitInfo)
# [1] 48 13;

# --@Crosschecks if the transformation was done correctly;
TemChkCols <- c("Number_of_input_reads", "Number_of_Uniquely_mapped_reads", "Ratio_of_uniquely_mapped_reads", "RIN", "LINKERS", "AGE", "GENDER", "PMD.MIN.", "PH")
unique(unlist(
  lapply(seq(TemChkCols), function(i)
    if(is.numeric(reqFroAlCaseVsCont_TraitInfo[, TemChkCols[i]])) {
      identical(finFroAlCaseVsCont_TraitInfo[, TemChkCols[i]], as.numeric(reqFroAlCaseVsCont_TraitInfo[, TemChkCols[i]]))
    } else {
      identical(finFroAlCaseVsCont_TraitInfo[, TemChkCols[i]], as.numeric(factor(reqFroAlCaseVsCont_TraitInfo[, TemChkCols[i]])))
    })
))
# [1] TRUE;  

# ---@keep only the required objects;
length(ls())
# [1] 32;
rm(list = ls()[-(match(c("reqFro_alCasesVsCont_normd_logcpm_Trs", "Fro_alCasesVsCont_nGenes", "Fro_alCasesVsCont_nSamples", "reqFro_AlCaseVSCont_SampIDs", "reqFroAlCaseVsCont_TraitInfo", "finFroAlCaseVsCont_TraitInfo"), ls()))])
length(ls())
# [1] 6;

collectGarbage();



# ######################################;
# II. Unsigned network topology analysis;
# ______________________________________;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Fron_WGCNA/Unsigned/All_CasesVSCont/Step2")

# 1. Construct weighted gene network: choose soft thresholding power beta to be raised to coexpression similarity to calculate adjacency;
powers <- 1:20
cex1 <- 0.9;

# -@Recommends: Choose the smallest power for which R^2>0.8 or if a saturation curve results, choose the power at the kind of the saturation curve;
sft_UnSign <- pickSoftThreshold(reqFro_alCasesVsCont_normd_logcpm_Trs, powerVector = powers, networkType = "unsigned", verbose = 5)

# a.1) plot the result;
h1 <- 0.8781268
pdf("Fro_alCasesVsCont_SFT_UnSign", width = 14, height = 8)
par(mfrow = c(1,2));
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft_UnSign$fitIndices[,1], -sign(sft_UnSign$fitIndices[,3])*sft_UnSign$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
text(sft_UnSign$fitIndices[,1], -sign(sft_UnSign$fitIndices[,3])*sft_UnSign$fitIndices[,2], labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=h1,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft_UnSign$fitIndices[,1], sft_UnSign$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft_UnSign$fitIndices[,1], sft_UnSign$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# -@get the Rsq value for the chosen beta value;
sft_UnSignplotDF <- data.frame(pwr = sft_UnSign$fitIndices[,1], RSq = -sign(sft_UnSign$fitIndices[,3])*sft_UnSign$fitIndices[,2])
dim(sft_UnSignplotDF)
# [1] 20  2;
sft_UnSignplotDF[7:8,]
#   pwr       RSq
#   7   7 0.8781268
#   8   8 0.8869769

# a.1.1) inspect scalefree plot for both: beta = 7; and, beta = 8;
Fro_alCasesVsCont_k7 <- softConnectivity(reqFro_alCasesVsCont_normd_logcpm_Trs, power = 7, type = "unsigned")
Fro_alCasesVsCont_k8 <- softConnectivity(reqFro_alCasesVsCont_normd_logcpm_Trs, power = 8, type = "unsigned")

pdf("Fro_alCasesVsCont_UnSign_ScaleFreePlot", width = 14, height = 8)
par(mfrow = c(1,2))
scaleFreePlot(Fro_alCasesVsCont_k7, main = "Scale free plot(beta=7) - unsigned Fro(all cases and controls)\n", truncated = T)
scaleFreePlot(Fro_alCasesVsCont_k8, main = "Scale free plot(beta=8) - unsigned Fro(all cases and controls)\n", truncated = T)
dev.off()

# a.2) calculate coexpression similarity and adjacency using the chosen soft thresholding power;
sft_UnSign_Pwr <- 7;
# Note: chose this power coz-least power @saturation-mean connectivity >power8-Rsq almost 0.9;
Fro_alCasesVsCont_Adj <- adjacency(reqFro_alCasesVsCont_normd_logcpm_Trs, type = "unsigned", power = sft_UnSign_Pwr)
dim(Fro_alCasesVsCont_Adj)
# [1] 21051 21051;

# a.3) convert adjacency into TOM in order to avoid noisy association and calculate the dissimilarities;
Fro_alCasesVsCont_Adj2TOM <- TOMsimilarity(Fro_alCasesVsCont_Adj, TOMType = "unsigned")
dim(Fro_alCasesVsCont_Adj2TOM)
# [1] 21051 21051;

Fro_alCasesVsCont_Adj2TOM_disim <- 1 - Fro_alCasesVsCont_Adj2TOM
dim(Fro_alCasesVsCont_Adj2TOM_disim)
# [1] 21051 21051;

# a.4) clustering using TOM;
Fro_alCasesVsCont_Adj2TOM_disim_geneTre <- hclust(as.dist(Fro_alCasesVsCont_Adj2TOM_disim), method = "average")

# a.4.1 plot the results;
pdf("Fro_alCasesVsCont_TOMdisim_geneTre_UnSign", width = 10, height = 6)
plot(Fro_alCasesVsCont_Adj2TOM_disim_geneTre, xlab="", sub="", main = "Frontal: gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

# a.5) module identification using dynamic tree cut;
Fro_alCasesVsCont_dynamicMods <- cutreeDynamic(dendro = Fro_alCasesVsCont_Adj2TOM_disim_geneTre, distM = Fro_alCasesVsCont_Adj2TOM_disim, deepSplit = 3, pamRespectsDendro = FALSE, minClusterSize = 50)
sort(table(Fro_alCasesVsCont_dynamicMods))

# a.5.1 convert numeric labels into colors;
Fro_alCasesVsCont_dynamicMods_Cols <- labels2colors(Fro_alCasesVsCont_dynamicMods)
table(Fro_alCasesVsCont_dynamicMods_Cols)

# a.5.2 plot the dendogram and colors underneath;
pdf("Fro_alCasesVsCont_TOMdisim_geneTre_min50ModMap", width = 15, height = 10)
plotDendroAndColors(Fro_alCasesVsCont_Adj2TOM_disim_geneTre, Fro_alCasesVsCont_dynamicMods_Cols, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Frontal-brain: gene dendrogram and module colors")
dev.off()

# a.6) Merge modules with similar expression profiles; 

# a.6.1 Calculate eigengenes;
reqFro_alCasesVsCont_normd_logcpm_Trs_MELs <- moduleEigengenes(reqFro_alCasesVsCont_normd_logcpm_Trs, colors = Fro_alCasesVsCont_dynamicMods_Cols, softPower = sft_UnSign_Pwr)
reqFro_alCasesVsCont_normd_logcpm_Trs_MEs <- reqFro_alCasesVsCont_normd_logcpm_Trs_MELs$eigengenes
dim(reqFro_alCasesVsCont_normd_logcpm_Trs_MEs)
# [1] 48 43;

# a.6.2 Calculate dissimilarity of module eigengenes
reqFro_alCasesVsCont_normd_logcpm_Trs_MEDiss <- 1-cor(reqFro_alCasesVsCont_normd_logcpm_Trs_MEs)
dim(reqFro_alCasesVsCont_normd_logcpm_Trs_MEDiss)
# [1] 43 43;

# a.6.3 Cluster module eigengenes;
reqFro_alCasesVsCont_normd_logcpm_Trs_disMETree <- hclust(as.dist(reqFro_alCasesVsCont_normd_logcpm_Trs_MEDiss), method = "average");

# @Plot the module eigengenes cluster;
pdf("unSgnFro_alCasesVsCont_MEsclusters_min50", width = 15, height = 10)
plot(reqFro_alCasesVsCont_normd_logcpm_Trs_disMETree, main = "Frontal: clustering of module eigengenes", xlab = "", sub = "")
MEDissThres = 0.1
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# a.6.4 Call an automatic merging function;
Fro_alCasesVsCont_normlogcpmTrs_meMrg <- mergeCloseModules(reqFro_alCasesVsCont_normd_logcpm_Trs, Fro_alCasesVsCont_dynamicMods_Cols, cutHeight = MEDissThres, verbose = 3)
dim(Fro_alCasesVsCont_normlogcpmTrs_meMrg$newMEs)
# [1] 48 39;

# @the merged module colors
Fro_alCasesVsCont_normlogcpmTrs_meMrgCols <- Fro_alCasesVsCont_normlogcpmTrs_meMrg$colors;

# @Eigengenes of the new merged modules:
Fro_alCasesVsCont_normlogcpmTrs_meMrg_MEs <- Fro_alCasesVsCont_normlogcpmTrs_meMrg$newMEs;
dim(Fro_alCasesVsCont_normlogcpmTrs_meMrg_MEs)
# [1] 48 39;


# 5. visualize the clusters along with the profile before and after merging the modules;
pdf("Fro_alCasesVsCont_TOMdisimgeneTre_b4afModMrg_min50Cpt1", width = 15, height = 10)
plotDendroAndColors(Fro_alCasesVsCont_Adj2TOM_disim_geneTre, cbind(Fro_alCasesVsCont_dynamicMods_Cols, Fro_alCasesVsCont_normlogcpmTrs_meMrgCols), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Frontal: cluster dendrogram")
dev.off()



# #######################;
# III. Save the workspace;
# _______________________;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc")
save.image("WGCNA_S2_FronSamps_AllCasesNConts_CAGEseq.RData")
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc/WGCNA_S2_FronSamps_AllCasesNConts_CAGEseq.RData")
# END ____________________________________________________________________________________________________________________________________________________________________________________________;