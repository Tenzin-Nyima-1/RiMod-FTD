# Date: 17-18 April 2018;
# Project: RIMOD-FTD;
# Aim: Weighted Gene Coexpression Network Analysis(WGCNA) of frontal samples: all cases and control samples;
# NOTE: the annotation of the count-table was carried using RIKEN's CAT gtf file as the gene model;
# NOTE 2: the same script could be adapted for the count table from: all 7 brain regions; and, temporal brain; 
# Workflow (main steps);
# Step 4: Unsigned - checking intramodular connectivity;


options(stringsAsFactors=FALSE)
library(WGCNA)
library(psych)
library(stringr)
# *
enableWGCNAThreads()


# ##############################################################################################################;
# I. Load the workspace from S3: GS vs MM;
# ______________________________________________________________________________________________________________;
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc/WGCNA_S3_FronSamps_AllCasesNConts_CAGEseq.RData")



# ##############################################################################################################;
# II. Intramodular connectivity analysis;
# ______________________________________________________________________________________________________________;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Fron_WGCNA/Unsigned/All_CasesVSCont/Step3_GSvsIntModConn")

# 1. intramodular connectivity;
Fro_alCasesVsCont_intrModConn <- intramodularConnectivity(adjMat = Fro_alCasesVsCont_Adj, colors = Fro_alCasesVsCont_normlogcpmTrs_meMrgCols)
dim(Fro_alCasesVsCont_intrModConn)
# [1] 21051     4;


# [A] Modules significantly correlated to the biological traits: all cases vs cont and each group categorized separately ###;

# 1. Recall: only one module found significantly correlated to the each gp;
length(Fro_bothStats_FDRsigMEnm)
# [1] 2;

# @save the plots;
for(i in seq(Fro_bothStats_FDRsigMEnm)){
  pdfNms <- unlist(lapply(seq(Fro_bothStats_FDRsigMEnm), function(x) paste("unsgnFro_eachCasesNCont_", Fro_bothStats_FDRsigMEnm[x], "_intKvsGS_min50Cpt1.pdf", sep = "")))
  pdf(pdfNms[i], width = 12, height = 7)
  par(mfrow = c(1,2))
  #P1: all cases combined and controls;
  verboseScatterplot(x = Fro_alCasesVsCont_intrModConn[Fro_CorbothTrt_FDRSig_Mod_ClusIDs[[Fro_bothStats_FDRsigMEnm[i]]], "kWithin"],
                     y = abs(Fro_FDRSigCor_Mod_GSolClusNCont[[Fro_bothStats_FDRsigMEnm[i]]]),
                     main = paste("Frontal brain: module - ", Fro_bothStats_FDRsigMEnm[i], " (n=", length(Fro_FDRSigCor_Mod_GSeachClusNCont[[Fro_bothStats_FDRsigMEnm[i]]]), ")", "\n", sep = ""),
                     xlab = "Intramodular connectivity",
                     ylab= "Gene significance in all cases combined and controls",
                     col = newFro_bothStats_FDRsigMEnm[i], cex.main = 1, cex.lab = 1, cex.axis = 1, cex = 0.5)
  text(x = Fro_alCasesVsCont_intrModConn[Fro_CorbothTrt_FDRSig_Mod_ClusIDs[[Fro_bothStats_FDRsigMEnm[i]]], "kWithin"],
       y = abs(Fro_FDRSigCor_Mod_GSolClusNCont[[Fro_bothStats_FDRsigMEnm[i]]]),
       labels = Fro_alCasesVsCont_WGCNA_GeneInf3[Fro_CorbothTrt_FDRSig_Mod_ClusIDs[[Fro_bothStats_FDRsigMEnm[i]]], "geneName"], 
       cex = 0.5, pos = 3, font = 2, col = "mediumvioletred")
  #P2: each cases grouped individually and controls;
  verboseScatterplot(x = Fro_alCasesVsCont_intrModConn[Fro_CorbothTrt_FDRSig_Mod_ClusIDs[[Fro_bothStats_FDRsigMEnm[i]]], "kWithin"],
                     y = abs(Fro_FDRSigCor_Mod_GSeachClusNCont[[Fro_bothStats_FDRsigMEnm[i]]]),
                     main = paste("Frontal brain: module - ", Fro_bothStats_FDRsigMEnm[i], " (n=", length(Fro_FDRSigCor_Mod_GSeachClusNCont[[Fro_bothStats_FDRsigMEnm[i]]]), ")", "\n", sep = ""),
                     xlab = "Intramodular connectivity",
                     ylab= "Gene significance in each cases grouped separately and controls",
                     col = newFro_bothStats_FDRsigMEnm[i], cex.main = 1, cex.lab = 1, cex.axis = 1, cex = 0.5)
  text(x = Fro_alCasesVsCont_intrModConn[Fro_CorbothTrt_FDRSig_Mod_ClusIDs[[Fro_bothStats_FDRsigMEnm[i]]], "kWithin"],
       y = abs(Fro_FDRSigCor_Mod_GSeachClusNCont[[Fro_bothStats_FDRsigMEnm[i]]]),
       labels = Fro_alCasesVsCont_WGCNA_GeneInf3[Fro_CorbothTrt_FDRSig_Mod_ClusIDs[[Fro_bothStats_FDRsigMEnm[i]]], "geneName"], 
       cex = 0.5, pos = 3, font = 2, col = "mediumvioletred")
  dev.off()
}


# [B] Modules containing key FTD genes #########################################################;

# 1. Generate a new anno df with the k-info, but without the MM cols;
Fro_alCasesVsCont_GeneInf_intK_df <- cbind(Fro_alCasesVsCont_WGCNA_GeneInf3[,-(grep("MM", colnames(Fro_alCasesVsCont_WGCNA_GeneInf3)))], 
                                           Fro_alCasesVsCont_intrModConn[rownames(Fro_alCasesVsCont_WGCNA_GeneInf3),])
dim(Fro_alCasesVsCont_GeneInf_intK_df)
# [1] 21051    27;

# @crosschecks;
identical(rownames(Fro_alCasesVsCont_WGCNA_GeneInf3), as.character(Fro_alCasesVsCont_WGCNA_GeneInf3$Clus_ID))
# [1] TRUE;
identical(rownames(Fro_alCasesVsCont_WGCNA_GeneInf3), rownames(Fro_alCasesVsCont_intrModConn[rownames(Fro_alCasesVsCont_WGCNA_GeneInf3),]))
# [1] TRUE;

# 2. Recall important variables;
length(Fro_KeyFTDGenes_asocModNms)
# [1] 13;

# 3. call all the functions I designed to get values specific for each cluster IDs repestively for each module;
source("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_Scripts/TenzinNyima_functions4WGCNA_analysis.R")

# *use the function (Fro_KeyFTDGenesasocMod_AnnoIntkDF);
Fro_alCasesVsCont_GeneInfintK_lsPMod <- Fro_KeyFTDGenesasocMod_AnnoIntkDF(ModNames = Fro_KeyFTDGenes_asocModNms, ModAnnoDF = Fro_alCasesVsCont_GeneInf_intK_df)
length(Fro_alCasesVsCont_GeneInfintK_lsPMod)
# [1] 13;

unlist(lapply(Fro_alCasesVsCont_GeneInfintK_lsPMod, nrow))
# blue  lightgreen greenyellow      purple        pink       green      orange         red   turquoise       brown       black darkmagenta   royalblue; 
# 5475         187         316         326         363         651         798         566        6393         918         480         101         178;

# @crosscheck;
identical(names(Fro_alCasesVsCont_GeneInfintK_lsPMod), Fro_KeyFTDGenes_asocModNms)
# [1] TRUE;
unique(unlist(lapply(seq(Fro_KeyFTDGenes_asocModNms), function(x) 
  identical(rownames(Fro_KeyFTDGenes_AnnoDF_perModls[[x]]), rownames(Fro_alCasesVsCont_GeneInfintK_lsPMod[[x]])) )))
# [1] TRUE;

# 4. execute the function (NonHubGens_convBlnk) to label only genes with the top 10% intramodular connectiviting values;
Fro_alCasesVsCont_top20percK_genLab_ls <- NonHubGens_convBlnk(KeyGens_AnnDFls = Fro_alCasesVsCont_GeneInfintK_lsPMod, keyGeNms = KeyFTDGenes)
length(Fro_alCasesVsCont_top20percK_genLab_ls)
# [1] 13;

# @crosscheck;
identical(names(Fro_alCasesVsCont_top20percK_genLab_ls), Fro_KeyFTDGenes_asocModNms)
# [1] TRUE;

# @save the plot;
for(i in seq(Fro_KeyFTDGenes_asocModNms)){
  pdfNms <- unlist(lapply(seq(Fro_KeyFTDGenes_asocModNms), function(x) paste("unsgnFro_KeyFTDGenes_", Fro_KeyFTDGenes_asocModNms[x], "_intKvsGS_min50Cpt1.pdf", sep = "")))
  pdf(pdfNms[i], width = 12, height = 7)
  par(mfrow = c(1, 2))
  #P1: all cases combined and controls;
  verboseScatterplot(x = Fro_alCasesVsCont_intrModConn[Fro_KeyFTDGenes_asocMod_ClusIDs[[Fro_KeyFTDGenes_asocModNms[i]]], "kWithin"],
                     y = abs(Fro_KeyFTDGenes_asocMod_GSolClusNCont[[Fro_KeyFTDGenes_asocModNms[i]]]),
                     xlab = "Intramodular connectivity",
                     ylab = paste("Gene significance in all cases combined and controls"),
                     main = paste("Frontal brain: module - ", Fro_KeyFTDGenes_asocModNms[i], " (n=", length(Fro_KeyFTDGenes_asocMod_GSolClusNCont[[Fro_KeyFTDGenes_asocModNms[i]]]), ")", "\n", sep = ""),
                     col = Fro_KeyFTDGenes_asocModNms[i], cex.main = 1, cex.lab = 1, cex.axis = 1, cex = 0.5)
  text(x = Fro_alCasesVsCont_intrModConn[Fro_KeyFTDGenes_asocMod_ClusIDs[[Fro_KeyFTDGenes_asocModNms[i]]], "kWithin"],
       y = abs(Fro_KeyFTDGenes_asocMod_GSolClusNCont[[Fro_KeyFTDGenes_asocModNms[i]]]),
       labels = Fro_alCasesVsCont_top20percK_genLab_ls[[Fro_KeyFTDGenes_asocModNms[i]]][, "geneName"], 
       col = Fro_alCasesVsCont_top20percK_genLab_ls[[Fro_KeyFTDGenes_asocModNms[i]]][,"LablCol"],
       cex = 0.5, pos = 3, font = 2)
  #P2: each cases grouped individually and controls;
  verboseScatterplot(x = Fro_alCasesVsCont_intrModConn[Fro_KeyFTDGenes_asocMod_ClusIDs[[Fro_KeyFTDGenes_asocModNms[i]]], "kWithin"],
                     y = abs(Fro_KeyFTDGenes_asocMod_GSeachClusNCont[[Fro_KeyFTDGenes_asocModNms[i]]]),
                     xlab = "Intramodular connectivity",
                     ylab = paste("Gene significance in each individual cases grouped separately and controls"),
                     main = paste("Frontal brain: module - ", Fro_KeyFTDGenes_asocModNms[i], " (n=", length(Fro_KeyFTDGenes_asocMod_GSeachClusNCont[[Fro_KeyFTDGenes_asocModNms[i]]]), ")", "\n", sep = ""),
                     col = Fro_KeyFTDGenes_asocModNms[i], cex.main = 1, cex.lab = 1, cex.axis = 1, cex = 0.5)
  text(x = Fro_alCasesVsCont_intrModConn[Fro_KeyFTDGenes_asocMod_ClusIDs[[Fro_KeyFTDGenes_asocModNms[i]]], "kWithin"],
       y = abs(Fro_KeyFTDGenes_asocMod_GSeachClusNCont[[Fro_KeyFTDGenes_asocModNms[i]]]),
       labels = Fro_alCasesVsCont_top20percK_genLab_ls[[Fro_KeyFTDGenes_asocModNms[i]]][, "geneName"], 
       col = Fro_alCasesVsCont_top20percK_genLab_ls[[Fro_KeyFTDGenes_asocModNms[i]]][,"LablCol"],
       cex = 0.5, pos = 3, font = 2)
  dev.off()
}



# ##############################################################################################################;
# III. Save this workspace;
# ______________________________________________________________________________________________________________;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc")
save.image("WGCNA_S4_FronSamps_AllCasesNConts_CAGEseq.RData")
#load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc/WGCNA_S4_FronSamps_AllCasesNConts_CAGEseq.RData")
# End ___________________________________________________________________________________________________________________________________________________________________________________________;
