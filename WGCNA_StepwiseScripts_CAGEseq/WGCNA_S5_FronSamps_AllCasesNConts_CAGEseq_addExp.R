# Date: 23-24 April, 2018;
# Project: RIMOD-FTD;
# Aim: Weighted Gene Coexpression Network Analysis(WGCNA) of frontal samples: all cases and control samples;
# NOTE: the annotation of the count-table was carried using RIKEN's CAT gtf file as the gene model;
# NOTE 2: the same script could be adapted for the count table from: all 7 brain regions; and, temporal brain; 
# Workflow (main steps);
# Step 4: Unsigned - additional exploration of the data;

options(stringsAsFactors=FALSE)
library(WGCNA)
library(psych)
library(stringr)
# *
enableWGCNAThreads()



# ##############################################################################################################;
# I. Load the previous workspace;
# ______________________________________________________________________________________________________________;
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc/WGCNA_S4_FronSamps_AllCasesNConts_CAGEseq.RData")

# @set the working dir;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Fron_WGCNA/Unsigned/All_CasesVSCont/Step5_addExplor")



# ##############################################################################################################;
# II. Generate the MDS plot;
# ______________________________________________________________________________________________________________;

# 1. generate a 2-D multi-dimensional scaling plot;
Fro_unsgn_AllCasesNCont <- cmdscale(d = as.dist(Fro_alCasesVsCont_Adj2TOM_disim), k = 2)
dim(Fro_unsgn_AllCasesNCont)
# [1] 21051     2;

# 2. save the plot;
pdf("Fro_AllCasesNCont_MDS2d_min50Cpt1.pdf")
plot(Fro_unsgn_AllCasesNCont, cex = 0.5, col = as.character(Fro_alCasesVsCont_normlogcpmTrs_meMrgCols), main = "Frontal brain: MDS plot", xlab = "Scaling dimension 1", ylab = "Scaling dimension 2")
dev.off()



# ##############################################################################################################;
# III. Checking the signed KME outputs as well;
# ______________________________________________________________________________________________________________;
# NOTE: intramodular connectivity for all genes on the array as opposed to intramodular connectivity that is limited to within a module only;
# NOTE2: to find genes with high gene significance and high intramodular connectivity;

# 1. get genes with high GS and high intramodular connectivity;
Fro_alCasesVsCont_ink <- Fro_alCasesVsCont_GeneInf_intK_df[unlist(lapply(seq(Fro_bothStats_FDRsigMEnm), 
                                                                         function(x) which(Fro_alCasesVsCont_GeneInf_intK_df$MrgdMECols == Fro_bothStats_FDRsigMEnm[x]))),]
dim(Fro_alCasesVsCont_ink)
# [1] 372  27;

# Fro_alCasesVsCont_ink hub using: GS vs intramodular connectiity;
(Fro_alCasesVsCont_inkgens <- Fro_alCasesVsCont_ink$geneName[with(Fro_alCasesVsCont_ink, which(kWithin > 10 & abs(GSStat_MutCat) > 0.4))])
# [1] "S100A10"  "VAMP5"    "ID3"      "FXYD5"    "TGFBR3"   "NUPR1"    "VIM"      "C1R"      "PLTP"     "LTBR"     "RIT1"     "RAB13"    "TMBIM1"   "SERPING1" "MYL12A"   "PAWR"     "TFPI"     "AKAP13"   "CLIC1"   
# [20] "PLP2"     "UBE2L6"   "ITPKB"    "LMOD1"    "TENC1"    "TIMP1"    "MARCH3"   "PLIN2"    "NAGA"     "HSPB1"    "HSD3B7"   "C1RL"     "CLIC1"    "TTC23"    "MSN"      "TMED10"   "IFI35"    "DTNA"     "AKR1C3"  
# [39] "MAOB"     "PELO"     "CLEC14A"  "NALCN"    "SAP30BP"  "CD109"    "RAB34"    "MARK2"    "STOM"     "IL13RA1"  "CMTM6"    "CD63"     "SLC9A9";

# 2. Generalizing intramodular connectivity for all genes on the array (when simulated manual was followed);
# @Signed eigengene-based connectivity or module membership;
Fro_alCasesVsCont_kme <- signedKME(datExpr = reqFro_alCasesVsCont_normd_logcpm_Trs, Fro_alCasesVsCont_mrgcolMEs, outputColumnName = "MM")
dim(Fro_alCasesVsCont_kme)
# [1] 21051    39;

# @find the interesting ones;
Filter_skyblue3Genes <- which(abs(Fro_eachCasesNCont_geneTrtSig) > .3 & abs(Fro_alCasesVsCont_kme[match("MMskyblue3", colnames(Fro_alCasesVsCont_kme))]) > .8)
length(Filter_skyblue3Genes)
# [1] 13;

# @crosscheck;
identical(str_sub(rownames(Fro_eachCasesNCont_geneTrtSig), 1, -2), str_sub(rownames(Fro_alCasesVsCont_kme), 1, -2))
# [1] TRUE;

# extract the hub genes;
Filter_skyblue3Genes_hubDF <- Fro_alCasesVsCont_GeneInf_intK_df[rownames(Fro_eachCasesNCont_geneTrtSig),][Filter_skyblue3Genes,]
# *it was written in the part 7 of the simulated data manual that;
# Note: that the definition does not require that gene i is a member of the particular module(here-mediumpurple);

# 3. the above(signedKME) seem to use the same explanation, so checking if previously calculated MM value matches the above;
fro_MMcopy <- Fro_alCasesVsCont_geneModMem
rownames(fro_MMcopy) <- str_sub(rownames(fro_MMcopy), 1, -2)
fro_kmecopy <- Fro_alCasesVsCont_kme
rownames(fro_kmecopy) <- str_sub(rownames(fro_kmecopy), 1, -2)

# @crosscheck;
identical(fro_MMcopy, fro_kmecopy)
# [1] TRUE
# NOTE: So, both have the same values;
# NOTE: either way outputs the same values, hence, the previous method is followed;
# ================================================================================;



# ##############################################################################################################;
# IV. Correlation of MEs with biological traits;
# ______________________________________________________________________________________________________________;

# 1. recall the ordered ME values;
# Fro_alCasesVsCont_mrgcolMEs <- orderMEs(Fro_alCasesVsCont_mrgcolMEs0)
dim(Fro_alCasesVsCont_mrgcolMEs)
# [1] 48 39;

# 2. MEs and the traits into one ordered DF;
Fro_allCasesNCont_ordMEnT <- orderMEs(cbind(Fro_alCasesVsCont_mrgcolMEs, Fro_alCasesVsCont_stat, Fro_alMutCat))
dim(Fro_allCasesNCont_ordMEnT)
# [1] 48 41;

# 3. plot the dendogram and the heatmap together;
pdf("unsgnFro_AllCasesnCont_MET_cor_min50Cpt1.pdf", width = 12, height = 15)
plotEigengeneNetworks(multiME = Fro_allCasesNCont_ordMEnT, 
                      marDendro = c(0,2,1,1),
                      marHeatmap = c(4,4,0,0),
                      setLabels = "",
                      plotDendrograms = T,
                      plotHeatmaps = T,
                      cex.lab = 0.5)
dev.off()



# ##############################################################################################################;
# V. Measure of module significance as a measure of gene significance;
# ______________________________________________________________________________________________________________;

# 1. combine all the modules that are interesting;
Fro_intMod_Nms <- c(Fro_bothStats_FDRsigMEnm, Fro_KeyFTDGenes_asocModNms)  
length(Fro_intMod_Nms)  
# [1] 15;

# @extract the df for the interested mods only;
Fro_intMod_alCasesNCont_GS <- Fro_alCasesVsCont_GeneInf_intK_df[unlist(lapply(seq(Fro_intMod_Nms), function(x) which(Fro_alCasesVsCont_GeneInf_intK_df$MrgdMECols == Fro_intMod_Nms[x]))), ]
dim(Fro_intMod_alCasesNCont_GS)
# [1] 17124    27;

# 2. module significance defined as a measure of GS;

# @GS-all cases combined and controls;
Fro_intMod_alCasesNCont_modsigasGS  <- with(Fro_intMod_alCasesNCont_GS, tapply(abs(GSStat_CaseCont), factor(MrgdMECols, levels = unique(MrgdMECols)), mean))
length(Fro_intMod_alCasesNCont_modsigasGS)
# [1] 15;

# @GS-each cases and controls;
Fro_intMod_eachCasesNCont_modsigasGS  <- with(Fro_intMod_alCasesNCont_GS, tapply(abs(GSStat_MutCat), factor(MrgdMECols, levels = unique(MrgdMECols)), mean))
length(Fro_intMod_eachCasesNCont_modsigasGS)
# [1] 15;

# @visualize;
pdf("Fro_intMod_alCasesNCont_modsigasGS_min50Cpt1.pdf", width = 30, height = 12)
par(mfrow = c(1,2))
with(Fro_intMod_alCasesNCont_GS, plotModuleSignificance(geneSignificance = abs(GSStat_CaseCont),  
                                                        colors = factor(MrgdMECols, levels = unique(MrgdMECols)),
                                                        ylim = c(0, 0.40),
                                                        main = "Frontal brain: gene significance across modules, ",
                                                        ylab = "Gene significance in all cases combined and controls", cex.names = 0.55))
abline(v = 2.5, col = "darkorchid4", lty = 4)
with(Fro_intMod_alCasesNCont_GS, plotModuleSignificance(geneSignificance = abs(GSStat_MutCat),  
                                                        colors = factor(MrgdMECols, levels = unique(MrgdMECols)),
                                                        ylim = c(0, 0.40),
                                                        main = "Frontal brain: gene significance across modules, ",
                                                        ylab = "Gene significance in each cases separately and controls", cex.names = 0.55))
abline(v = 2.5, col = "darkorchid4", lty = 4)
dev.off()



# ##############################################################################################################;
# VI. Save this workspace;
# ______________________________________________________________________________________________________________;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc")
save.image("WGCNA_S5_FronSamps_AllCasesNConts_CAGEseq_addExp.RData")
#load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc/WGCNA_S5_FronSamps_AllCasesNConts_CAGEseq_addExp.RData")
# END __________________________________________________________________________________________________________________________________________________________________________________________________;

