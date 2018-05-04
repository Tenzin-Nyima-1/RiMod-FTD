# Date: 17-18 April 2018;
# Project: RIMOD-FTD;
# Aim: Weighted Gene Coexpression Network Analysis(WGCNA) of frontal samples: all cases and control samples;
# NOTE: the annotation of the count-table was carried using RIKEN's CAT gtf file as the gene model;
# NOTE 2: the same script could be adapted for the count table from: all 7 brain regions; and, temporal brain; 
# Workflow (main steps);
# Step 3: Unsigned - relating modules to the external information and to identify hub genes;


options(stringsAsFactors=FALSE)
library(WGCNA)
library(psych)
# *
enableWGCNAThreads()



# ##############################################################################################################;
# I. Load the workspace that constructed the network;
# ______________________________________________________________________________________________________________;
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc/WGCNA_S2_FronSamps_AllCasesNConts_CAGEseq.RData")

# 1. Recall some vital objects;

# @the final transposed expression data(deleted one of the outlier sample);
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



# ##############################################################################################################;
# II. Unsigned Network: relate the modules to the external trait informations;
# ______________________________________________________________________________________________________________;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Fron_WGCNA/Unsigned/All_CasesVSCont/Step3_GSvsMM")

# 1. Recalculate MEs with the merged module colors;
# Fro_alCasesVsCont_normlogcpmTrs_meMrgCols <- Fro_alCasesVsCont_normlogcpmTrs_meMrg$colors;
Fro_alCasesVsCont_mrgcolMEs0 <- moduleEigengenes(reqFro_alCasesVsCont_normd_logcpm_Trs, colors = Fro_alCasesVsCont_normlogcpmTrs_meMrgCols, softPower = sft_UnSign_Pwr)$eigengenes
dim(Fro_alCasesVsCont_mrgcolMEs0)
# [1] 48 39;

# @order similar MEs together;
Fro_alCasesVsCont_mrgcolMEs <- orderMEs(Fro_alCasesVsCont_mrgcolMEs0)
dim(Fro_alCasesVsCont_mrgcolMEs)
# [1] 48 39;

# @crosscheck;
setdiff(colnames(Fro_alCasesVsCont_mrgcolMEs0), colnames(Fro_alCasesVsCont_mrgcolMEs))
# character(0);
setdiff(colnames(Fro_alCasesVsCont_mrgcolMEs), colnames(Fro_alCasesVsCont_mrgcolMEs0))
# character(0);
length(intersect(colnames(Fro_alCasesVsCont_mrgcolMEs0), colnames(Fro_alCasesVsCont_mrgcolMEs)))
# [1] 39;
identical(Fro_alCasesVsCont_mrgcolMEs, Fro_alCasesVsCont_normlogcpmTrs_meMrg_MEs)
# [1] TRUE #NOTE: the merged MEs in step 2 and the above(ordered) are exactly the same!;


# [A] Quantify module trait associations #########################################################################; 

# 1. Calculate the correlation between the MEs and the data-trait;
Fro_alCasesVsCont_MEsTrait_Cor <- cor(Fro_alCasesVsCont_mrgcolMEs, finFroAlCaseVsCont_TraitInfo, use = "p")
dim(Fro_alCasesVsCont_MEsTrait_Cor)
# [1] 39 13;

# # 1.1 a script to figure out how the pearson correlation was calculated: using one ME and the variable containing case-control status;
# MEdorg <-  Fro_alCasesVsCont_mrgcolMEs[,1] #variable-x;
# CaseCont_stat <- finFroAlCaseVsCont_TraitInfo[,1] #variable-x;
# ManuCorCalc <- data.frame(MEdorg = MEdorg, CaseCont_stat = CaseCont_stat, ME_Stat = MEdorg * CaseCont_stat, MEsq = MEdorg^2, Statsq = CaseCont_stat^2)
# ManuCorCalc[nrow(ManuCorCalc)+1,] <- colSums(ManuCorCalc)
# nSp <- nrow(ManuCorCalc)-1
# R_numer <- (nSp * ManuCorCalc$ME_Stat[nSp+1]) - (ManuCorCalc$MEdorg[nSp+1] * ManuCorCalc$CaseCont_stat[nSp+1])
# R_denom <- sqrt( (((nSp * ManuCorCalc$MEsq[nSp+1]) - (ManuCorCalc$MEdorg[nSp+1])^2) * ((nSp * ManuCorCalc$Statsq[nSp+1]) - (ManuCorCalc$CaseCont_stat[nSp+1])^2)))
# all.equal(R_numer/R_denom, Fro_alCasesVsCont_MEsTrait_Cor[1,1])
# # [1] TRUE;
# # 1.2 Also, check if the univ. addition of a value to the case-control numeric status changes the cor;
# all.equal(cor(Fro_alCasesVsCont_mrgcolMEs[,1], finFroAlCaseVsCont_TraitInfo[,1], use = "p"), cor(Fro_alCasesVsCont_mrgcolMEs[,1], finFroAlCaseVsCont_TraitInfo[,1]+10, use = "p"))
# # [1] TRUE; NOTE: seem to stay roughly the same: meaning only the trend is considered;

# 2. check the multiple test corrected significance level of the obtained correlation values;
Fro_alCasesVsCont_MEsTrait_CorPVal_fdr0 <- corr.p(Fro_alCasesVsCont_MEsTrait_Cor, n = Fro_alCasesVsCont_nSamples, adjust = "BH")
Fro_alCasesVsCont_MEsTrait_CorPVal_fdr <- Fro_alCasesVsCont_MEsTrait_CorPVal_fdr0$p
dim(Fro_alCasesVsCont_MEsTrait_CorPVal_fdr)
# [1] 39 13;

# @crosschecks: if ran without the pval adjustment renders the same result but using the above funtion?;
all.equal(corr.p(Fro_alCasesVsCont_MEsTrait_Cor, n = Fro_alCasesVsCont_nSamples, adjust = "none")$p, corPvalueStudent(Fro_alCasesVsCont_MEsTrait_Cor, Fro_alCasesVsCont_nSamples))
# [1] TRUE; NOTE: meaning that the function calculates the same raw pvalue just as the WGCNA tutorial;
identical(Fro_alCasesVsCont_MEsTrait_CorPVal_fdr0$r, Fro_alCasesVsCont_MEsTrait_Cor)
# [1] TRUE;

# 3. Visualize the corrlation value and the adjusted pvalue for each ME and the corresponding biological trait;
Fro_alCasesVsCont_MEsTrait_CorPVal_fdr_txtmat <-  paste(signif(Fro_alCasesVsCont_MEsTrait_Cor, 2), 
                                                        " (", 
                                                        signif(Fro_alCasesVsCont_MEsTrait_CorPVal_fdr, 1), 
                                                        ")", sep = "");
dim(Fro_alCasesVsCont_MEsTrait_CorPVal_fdr_txtmat) <- dim(Fro_alCasesVsCont_MEsTrait_Cor)
dim(Fro_alCasesVsCont_MEsTrait_CorPVal_fdr_txtmat)
# [1] 39 13;

# @visualize;
pdf("UnSgnFDR_Fro_alCasesVsCont_MEsTrait_CoradjPVal_min50Cpt1", width = 25, height = 20)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = Fro_alCasesVsCont_MEsTrait_Cor, xLabels = colnames(Fro_alCasesVsCont_MEsTrait_Cor), yLabels = rownames(Fro_alCasesVsCont_MEsTrait_Cor), ySymbols = rownames(Fro_alCasesVsCont_MEsTrait_Cor), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = Fro_alCasesVsCont_MEsTrait_CorPVal_fdr_txtmat, setStdMargins = FALSE, cex.text = 1, zlim = c(-1,1), main = paste("Frontal brain: module-trait relationships"), cex.lab.x = 0.7, cex.lab.y = 0.8)
dev.off()


# [B] Gene significance and Module Membership #######################################################################; 

# 1. Universal variables for the downstream analysis;

# - the case-control status;
Fro_alCasesVsCont_stat <- as.data.frame(finFroAlCaseVsCont_TraitInfo$CASE.CONTROL)
rownames(Fro_alCasesVsCont_stat) <- rownames(finFroAlCaseVsCont_TraitInfo)
colnames(Fro_alCasesVsCont_stat) <- "Stat_CaseCont"
dim(Fro_alCasesVsCont_stat)
# [1] 48  1;

# - all mutation category;
Fro_alMutCat <- as.data.frame(finFroAlCaseVsCont_TraitInfo$MUTATED.GENE)
rownames(Fro_alMutCat) <- rownames(finFroAlCaseVsCont_TraitInfo)
colnames(Fro_alMutCat) <- "Stat_MutCat"
dim(Fro_alMutCat)
# [1] 48  1;

# - just the colors of the merged and ordered ME;
Fro_alCasesVsCont_mrgMECols <- substring(names(Fro_alCasesVsCont_mrgcolMEs), 3)
length(Fro_alCasesVsCont_mrgMECols)
# [1] 39;

# 2. Module membership (Kme): cor bw gene expression profile & ME;

# *cor;
Fro_alCasesVsCont_geneModMem <- as.data.frame(cor(reqFro_alCasesVsCont_normd_logcpm_Trs, Fro_alCasesVsCont_mrgcolMEs, use = "p"))
colnames(Fro_alCasesVsCont_geneModMem) <- paste("MM", Fro_alCasesVsCont_mrgMECols, sep = "")
dim(Fro_alCasesVsCont_geneModMem)
# [1] 21051    39;

# # *significance test;
# Fro_alCasesVsCont_geneModMemPVal0 <- corr.p(as.matrix(Fro_alCasesVsCont_geneModMem), Fro_alCasesVsCont_nSamples, adjust = "BH")
# # note: too slow since >20k genes, hence following the method suggested by WGCNA package;
Fro_alCasesVsCont_geneModMemPVal <- as.data.frame(corPvalueStudent(as.matrix(Fro_alCasesVsCont_geneModMem), Fro_alCasesVsCont_nSamples))
colnames(Fro_alCasesVsCont_geneModMemPVal) <- paste("p.MM", Fro_alCasesVsCont_mrgMECols, sep = "")
dim(Fro_alCasesVsCont_geneModMemPVal)
# [1] 21051    39;    

# @crosscheck;
identical(substring(colnames(as.data.frame(cor(reqFro_alCasesVsCont_normd_logcpm_Trs, Fro_alCasesVsCont_mrgcolMEs, use = "p"))), 3), Fro_alCasesVsCont_mrgMECols)
# [1] TRUE;
identical(substring(names(Fro_alCasesVsCont_geneModMem), 3), Fro_alCasesVsCont_mrgMECols)
# [1] TRUE;
identical(substring(colnames(as.data.frame(corPvalueStudent(as.matrix(Fro_alCasesVsCont_geneModMem), Fro_alCasesVsCont_nSamples))), 3), Fro_alCasesVsCont_mrgMECols)
# [1] TRUE;
identical(substring(names(Fro_alCasesVsCont_geneModMemPVal), 5), Fro_alCasesVsCont_mrgMECols)
# [1] TRUE;


# 3. Gene significance: abs cor between the gene and the trait;

# 3.1 Gene significance with the CASE-Control status;

# *cor;
Fro_alCasesVsCont_geneTrtSig <- as.data.frame(cor(reqFro_alCasesVsCont_normd_logcpm_Trs, Fro_alCasesVsCont_stat, use = "p"))
colnames(Fro_alCasesVsCont_geneTrtSig) <- paste("GS", names(Fro_alCasesVsCont_stat), sep = "")
dim(Fro_alCasesVsCont_geneTrtSig)  
# [1] 21051     1;

# *significance test;
Fro_alCasesVsCont_geneTrtSigPVal <- as.data.frame(corPvalueStudent(as.matrix(Fro_alCasesVsCont_geneTrtSig), Fro_alCasesVsCont_nSamples))
colnames(Fro_alCasesVsCont_geneTrtSigPVal) <- paste("p.GS", names(Fro_alCasesVsCont_stat), sep = "")
dim(Fro_alCasesVsCont_geneTrtSigPVal)
# [1] 21051     1;

# 3.2 Gene significance with the mutated gene-category;

# *cor;
Fro_eachCasesNCont_geneTrtSig <- as.data.frame(cor(reqFro_alCasesVsCont_normd_logcpm_Trs, Fro_alMutCat, use = "p"))
colnames(Fro_eachCasesNCont_geneTrtSig) <- paste("GS", names(Fro_alMutCat), sep = "")
dim(Fro_eachCasesNCont_geneTrtSig)  
# [1] 21051     1;

# *significance test;
Fro_eachCasesNCont_geneTrtSigPVal <- as.data.frame(corPvalueStudent(as.matrix(Fro_eachCasesNCont_geneTrtSig), Fro_alCasesVsCont_nSamples))
colnames(Fro_eachCasesNCont_geneTrtSigPVal) <- paste("p.GS", names(Fro_alMutCat), sep = "")
dim(Fro_eachCasesNCont_geneTrtSigPVal)
# [1] 21051     1;



# [C] Summarize the outputs of the network analysis #############################################################################################################################; 

# 1. recall the normalised dgelist object;
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc/WGCNA_S1_FronSamps_AllCasesNConts_CAGEseq.RData")
dim(Fron_AlCasesVsCont_DGls_1in16_nrmd)
# [1] 21051    49;

# -crosscheck;
identical(rownames(Fron_AlCasesVsCont_DGls_1in16_nrmd), as.character(Fron_AlCasesVsCont_DGls_1in16_nrmd$genes$Clus_ID))
# [1] TRUE;

# 2. save the annotations for the tags;
# NOTE: it doesn't matter that one of the outlier sample is present in the above object, only focussing on the geneinfo list;
Fro_alCasesVsCont_WGCNA_GeneInf <- Fron_AlCasesVsCont_DGls_1in16_nrmd$genes[colnames(reqFro_alCasesVsCont_normd_logcpm_Trs), 
                                                                            c("Clus_ID", "width", "distanceToTSS", "annotation", "geneId", "transcriptId", "geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand", "geneName", "geneType", "CAT_geneClass", "CAT_DHS_type", "HGNC_name", "HGNC_symbol", "entrez_ID")]
dim(Fro_alCasesVsCont_WGCNA_GeneInf)
# [1] 21051    18;

# @crosscheck;
identical(rownames(Fro_alCasesVsCont_WGCNA_GeneInf), as.character(Fro_alCasesVsCont_WGCNA_GeneInf$Clus_ID))
# [1] TRUE;
identical(rownames(Fro_alCasesVsCont_WGCNA_GeneInf), colnames(reqFro_alCasesVsCont_normd_logcpm_Trs))
# [1] TRUE;
identical(rownames(Fro_alCasesVsCont_WGCNA_GeneInf), rownames(Fron_AlCasesVsCont_DGls_1in16_nrmd))
# [1] TRUE;
identical(rownames(Fro_alCasesVsCont_WGCNA_GeneInf), as.character(Fron_AlCasesVsCont_DGls_1in16_nrmd$genes$Clus_ID))
# [1] TRUE;

# 3. attach the gene-trait correlation and the associated significance values to the above DF;

# @a crosscheck before attaching the file; 
identical(rownames(Fro_alCasesVsCont_WGCNA_GeneInf), rownames(Fro_alCasesVsCont_geneTrtSig))
# [1] TRUE;
identical(rownames(Fro_alCasesVsCont_WGCNA_GeneInf), rownames(Fro_alCasesVsCont_geneTrtSigPVal))
# [1] TRUE;
identical(rownames(Fro_alCasesVsCont_WGCNA_GeneInf), rownames(Fro_eachCasesNCont_geneTrtSig))
# [1] TRUE;
identical(rownames(Fro_alCasesVsCont_WGCNA_GeneInf), rownames(Fro_eachCasesNCont_geneTrtSigPVal))
# [1] TRUE;

# *attach;
Fro_alCasesVsCont_WGCNA_GeneInf1 <- data.frame(Fro_alCasesVsCont_WGCNA_GeneInf, 
                                               MrgdMECols = Fro_alCasesVsCont_normlogcpmTrs_meMrgCols,
                                               Fro_alCasesVsCont_geneTrtSig[rownames(Fro_alCasesVsCont_WGCNA_GeneInf),], 
                                               Fro_alCasesVsCont_geneTrtSigPVal[rownames(Fro_alCasesVsCont_WGCNA_GeneInf),],
                                               Fro_eachCasesNCont_geneTrtSig[rownames(Fro_alCasesVsCont_WGCNA_GeneInf),], 
                                               Fro_eachCasesNCont_geneTrtSigPVal[rownames(Fro_alCasesVsCont_WGCNA_GeneInf),])
colnames(Fro_alCasesVsCont_WGCNA_GeneInf1) <- c(colnames(Fro_alCasesVsCont_WGCNA_GeneInf), 
                                                "MrgdMECols",
                                                colnames(Fro_alCasesVsCont_geneTrtSig), 
                                                colnames(Fro_alCasesVsCont_geneTrtSigPVal),
                                                colnames(Fro_eachCasesNCont_geneTrtSig), 
                                                colnames(Fro_eachCasesNCont_geneTrtSigPVal))
dim(Fro_alCasesVsCont_WGCNA_GeneInf1)
# [1] 21051    23;

# 4. Include an ordered MM and the corresponding pvalue to the aboce table;  

# 4.1) order(default, decreasing=F) MEs by their significance for the status: each cases vs cont separately;
Fro_alCasesVsCont_MEsMutGen_CorOrd <- order(-abs(as.data.frame(Fro_alCasesVsCont_MEsTrait_Cor)$MUTATED.GENE))
length(Fro_alCasesVsCont_MEsMutGen_CorOrd)
# [1] 39;

# 4.2) avoiding overwriting, save the above DF into a new variable;
Fro_alCasesVsCont_WGCNA_GeneInf2 <- Fro_alCasesVsCont_WGCNA_GeneInf1
dim(Fro_alCasesVsCont_WGCNA_GeneInf2)
# [1] 21051    23;

# 4.3) then in the identified order bind the MM and MMpvalue to the new variable one by one;
for(i in seq(ncol(Fro_alCasesVsCont_geneModMem))){
  Fro_alCasesVsCont_WGCNA_GeneInf2  <- data.frame(Fro_alCasesVsCont_WGCNA_GeneInf2, 
                                                  Fro_alCasesVsCont_geneModMem[, Fro_alCasesVsCont_MEsMutGen_CorOrd[i]],
                                                  Fro_alCasesVsCont_geneModMemPVal[, Fro_alCasesVsCont_MEsMutGen_CorOrd[i]])
}
dim(Fro_alCasesVsCont_WGCNA_GeneInf2)
# [1] 21051   101;

# *get the column names for the merged columns;
Fro_alCasesVsCont_WGCNA_ordMMpVAl_Nms <- unlist(lapply(seq(ncol(Fro_alCasesVsCont_geneModMem)), function(x) c(colnames(Fro_alCasesVsCont_geneModMem)[Fro_alCasesVsCont_MEsMutGen_CorOrd][x],
                                                                                                              colnames(Fro_alCasesVsCont_geneModMemPVal)[Fro_alCasesVsCont_MEsMutGen_CorOrd][x])))
length(Fro_alCasesVsCont_WGCNA_ordMMpVAl_Nms)
# [1] 78;

# *assign column names to the above main data frame;
colnames(Fro_alCasesVsCont_WGCNA_GeneInf2) <- c(colnames(Fro_alCasesVsCont_WGCNA_GeneInf1), Fro_alCasesVsCont_WGCNA_ordMMpVAl_Nms)
dim(Fro_alCasesVsCont_WGCNA_GeneInf2)
# [1] 21051   101;

# @crosschecks;
identical(rownames(Fro_alCasesVsCont_WGCNA_GeneInf2), as.character(Fro_alCasesVsCont_WGCNA_GeneInf2$Clus_ID))
# [1] TRUE;
identical(Fro_alCasesVsCont_WGCNA_GeneInf2[,colnames(Fro_alCasesVsCont_WGCNA_GeneInf1)], Fro_alCasesVsCont_WGCNA_GeneInf1)
# [1] TRUE;
identical(Fro_alCasesVsCont_WGCNA_GeneInf2[,colnames(Fro_alCasesVsCont_geneModMem)], Fro_alCasesVsCont_geneModMem)
# [1] TRUE;
identical(Fro_alCasesVsCont_WGCNA_GeneInf2[,colnames(Fro_alCasesVsCont_geneModMemPVal)], Fro_alCasesVsCont_geneModMemPVal)
# [1] TRUE;

# 4.4) order the genes by the module colour and then by the mutation category significance;
Fro_alCasesVsCont_WGCNA_GeneInf3 <- Fro_alCasesVsCont_WGCNA_GeneInf2[with(Fro_alCasesVsCont_WGCNA_GeneInf2, order(MrgdMECols, -abs(GSStat_MutCat))),]
dim(Fro_alCasesVsCont_WGCNA_GeneInf3)
# [1] 21051   101;

# @crosscheck;
identical(Fro_alCasesVsCont_WGCNA_GeneInf3[rownames(Fro_alCasesVsCont_WGCNA_GeneInf2),], Fro_alCasesVsCont_WGCNA_GeneInf2)
# [1] TRUE;
identical(rownames(Fro_alCasesVsCont_WGCNA_GeneInf3), as.character(Fro_alCasesVsCont_WGCNA_GeneInf3$Clus_ID))
# [1] TRUE;


# [C] Intramodular analysis: identifying genes with high GS and MM #######################################################################;

# 1. test code;
dem <- head(Fro_alCasesVsCont_MEsTrait_CorPVal_fdr)
dem_sig <- unlist(lapply(seq(nrow(dem)), function(x) paste(length(which(dem[x,1] < 0.05)),
                                                           length(which(dem[x,2] < 0.05)),
                                                           length(which(dem[x,3:ncol(dem)] < 0.05)), sep = "_") ))
grep("1_0_0", dem_sig)
# integer(0)
grep("1_1_0", dem_sig)
# [1] 3;
grep("0_1_0", dem_sig)
# integer(0)

# 2. identify MEs whose adjusted p-value < 0.05 (for the obtained correlation ME and traits);
dim(Fro_alCasesVsCont_MEsTrait_CorPVal_fdr)
# [1] 39 13;

# @generate 2 universe variables;
(nMods <- nrow(Fro_alCasesVsCont_MEsTrait_CorPVal_fdr)) #[1] 39;
(nTrts <- ncol(Fro_alCasesVsCont_MEsTrait_CorPVal_fdr)) #[1] 13;

# @identify significant modules;
Fro_alCasesVsCont_MECasCont_FDRSig <- unlist(lapply(seq(nMods), function(x) paste(length(which(Fro_alCasesVsCont_MEsTrait_CorPVal_fdr[x,1] < 0.05)),
                                                                                  length(which(Fro_alCasesVsCont_MEsTrait_CorPVal_fdr[x,2] < 0.05)),
                                                                                  length(which(Fro_alCasesVsCont_MEsTrait_CorPVal_fdr[x,3:nTrts] < 0.05)), sep = "_") ))
length(Fro_alCasesVsCont_MECasCont_FDRSig)
# [1] 39;

# 3. identify modules which is significantly correlated with either case-cont/mut-genCat status or both; 

# *only all case combined and cont status;
grep("1_0_0", Fro_alCasesVsCont_MECasCont_FDRSig, ignore.case = T)
# integer(0);

# *only mutation category status;
grep("0_1_0", Fro_alCasesVsCont_MECasCont_FDRSig, ignore.case = T)
# integer(0);

# *both statuses;
grep("1_1_0", Fro_alCasesVsCont_MECasCont_FDRSig, ignore.case = T)
# [1]  3 22;
(Fro_bothStats_FDRsigMEnm <- substring(rownames(Fro_alCasesVsCont_MEsTrait_CorPVal_fdr)[grep("1_1_0", Fro_alCasesVsCont_MECasCont_FDRSig, ignore.case = T)], 3))
# [1] "skyblue3"    "lightyellow";
length(Fro_bothStats_FDRsigMEnm)
# [1] 2;


# [D] Intramodular analysis: data prep for visualizations - GS vs MM #######################################################################;
# @call all the functions I designed to get values specific for each cluster IDs repestively for each module;
source("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_Scripts/TenzinNyima_functions4WGCNA_analysis.R")


# 1) Extract cluster IDs relevant to each modules ________________________________________________________________________________________; 

# .1 use the function (Extract_ClusIDs) to extract the cluster-ids of each individual module separately;

# .2 extract clus-ids using the above function;
Fro_CorbothTrt_FDRSig_Mod_ClusIDs <- Extract_ClusIDs(ModNames = Fro_bothStats_FDRsigMEnm, ModAnnoDF = Fro_alCasesVsCont_WGCNA_GeneInf3)
length(Fro_CorbothTrt_FDRSig_Mod_ClusIDs)  
# [1] 2;

# @crosscheck;
identical(names(Fro_CorbothTrt_FDRSig_Mod_ClusIDs), Fro_bothStats_FDRsigMEnm)
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_CorbothTrt_FDRSig_Mod_ClusIDs)), function(x)
  identical(Fro_CorbothTrt_FDRSig_Mod_ClusIDs[[x]], as.character(with(Fro_alCasesVsCont_WGCNA_GeneInf3, Clus_ID[which(MrgdMECols == names(Fro_CorbothTrt_FDRSig_Mod_ClusIDs)[x])]) )) )))
# [1] TRUE;

# 2) Extract MMvalues for the cluster-IDs in to each modules ________________________________________________________________________________; 

# .1 use the function (Extract_MM) to extract the cluster-ids of each individual module separately;

# .2 extract MM for clus-ids using the above function;
Fro_CorbothTrt_FDRSig_Modclus_MM <- Extract_MM(ModNames = names(Fro_CorbothTrt_FDRSig_Mod_ClusIDs), ModAnnoDF = Fro_alCasesVsCont_WGCNA_GeneInf3)
length(Fro_CorbothTrt_FDRSig_Modclus_MM)
# [1] 2;

# @crosscheck;
identical(names(Fro_CorbothTrt_FDRSig_Modclus_MM), names(Fro_CorbothTrt_FDRSig_Mod_ClusIDs))
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_CorbothTrt_FDRSig_Modclus_MM)), function(x) 
  identical(Fro_CorbothTrt_FDRSig_Modclus_MM[[x]], Fro_alCasesVsCont_WGCNA_GeneInf3[which(Fro_alCasesVsCont_WGCNA_GeneInf3$MrgdMECols == names(Fro_CorbothTrt_FDRSig_Modclus_MM)[x]),
                                                                                    match(paste("MM", names(Fro_CorbothTrt_FDRSig_Modclus_MM)[x], sep = ""), colnames(Fro_alCasesVsCont_WGCNA_GeneInf3))] ) )))
# [1] TRUE;

# @crosscheck using clusID names;
identical(names(Fro_CorbothTrt_FDRSig_Modclus_MM), names(Fro_CorbothTrt_FDRSig_Mod_ClusIDs))
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_CorbothTrt_FDRSig_Modclus_MM)), function(x)
  identical(Fro_CorbothTrt_FDRSig_Modclus_MM[[x]], Fro_alCasesVsCont_WGCNA_GeneInf3[Fro_CorbothTrt_FDRSig_Mod_ClusIDs[[x]],
                                                                                    match(paste("MM", names(Fro_CorbothTrt_FDRSig_Modclus_MM)[x], sep = ""), colnames(Fro_alCasesVsCont_WGCNA_GeneInf3))] ) )))
# [1] TRUE;

# 3) Extract GSvalues for significant modules _______________________________________________________________________________________________; 

# .1 use the function (Extract_GS_both) to extract GS-values for both traits: all cases combined and cont; and, each case n cont;


# 3.2 Modules significant only for the trait: each case and cont - extract their GS values;
Fro_FDRSigCor_Mod_GSolClusNCont <- Extract_GS_both(ModNames = Fro_bothStats_FDRsigMEnm, ModAnnoDF = Fro_alCasesVsCont_WGCNA_GeneInf3, TraitColNm = "GSStat_CaseCont")
length(Fro_FDRSigCor_Mod_GSolClusNCont)
# [1] 2;  

# @crosscheck;
identical(names(Fro_FDRSigCor_Mod_GSolClusNCont), Fro_bothStats_FDRsigMEnm)
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_FDRSigCor_Mod_GSolClusNCont)), function(x)
  identical(Fro_FDRSigCor_Mod_GSolClusNCont[[x]], Fro_alCasesVsCont_WGCNA_GeneInf3[which(Fro_alCasesVsCont_WGCNA_GeneInf3$MrgdMECols == names(Fro_FDRSigCor_Mod_GSolClusNCont)[x]), "GSStat_CaseCont"]) )))
# [1] TRUE;

# @crosscheck using the identified clus-ids;
identical(names(Fro_FDRSigCor_Mod_GSolClusNCont), names(Fro_CorbothTrt_FDRSig_Mod_ClusIDs))
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_FDRSigCor_Mod_GSolClusNCont)), function(x) 
  identical(Fro_FDRSigCor_Mod_GSolClusNCont[[x]], Fro_alCasesVsCont_WGCNA_GeneInf3[Fro_CorbothTrt_FDRSig_Mod_ClusIDs[[x]], "GSStat_CaseCont"]) )))
# [1] TRUE;


# 3.3 Modules(2) significant only for the trait: each case and cont - extract their GS values;
Fro_FDRSigCor_Mod_GSeachClusNCont <- Extract_GS_both(ModNames = Fro_bothStats_FDRsigMEnm, ModAnnoDF = Fro_alCasesVsCont_WGCNA_GeneInf3, TraitColNm = "GSStat_MutCat")
length(Fro_FDRSigCor_Mod_GSeachClusNCont)
# [1] 2;  

# @crosscheck;
identical(names(Fro_FDRSigCor_Mod_GSeachClusNCont), Fro_bothStats_FDRsigMEnm)
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_FDRSigCor_Mod_GSeachClusNCont)), function(x)
  identical(Fro_FDRSigCor_Mod_GSeachClusNCont[[x]], Fro_alCasesVsCont_WGCNA_GeneInf3[which(Fro_alCasesVsCont_WGCNA_GeneInf3$MrgdMECols == names(Fro_FDRSigCor_Mod_GSeachClusNCont)[x]), "GSStat_MutCat"]) )))
# [1] TRUE;

# @crosscheck using the identified clus-ids;
identical(names(Fro_FDRSigCor_Mod_GSeachClusNCont), names(Fro_CorbothTrt_FDRSig_Mod_ClusIDs))
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_FDRSigCor_Mod_GSeachClusNCont)), function(x) 
  identical(Fro_FDRSigCor_Mod_GSeachClusNCont[[x]], Fro_alCasesVsCont_WGCNA_GeneInf3[Fro_CorbothTrt_FDRSig_Mod_ClusIDs[[x]], "GSStat_MutCat"]) )))
# [1] TRUE;


# [E] Intramodular analysis: scatter plots(GS vs MM) for mods sig cor to mutation cateogry #######################################################################;

# 1. save names for the pdf files;
Fro_bothStats_FDRsigCor_MMvsGSnms <- unlist(lapply(seq(Fro_bothStats_FDRsigMEnm), 
                                                   function(x) paste("UnSgn_Fro_bothStats_",  Fro_bothStats_FDRsigMEnm[x], "_MMvsGS_min50Cpt1.pdf", sep = "")))

# *save;
newFro_bothStats_FDRsigMEnm <- ifelse(Fro_bothStats_FDRsigMEnm == "lightyellow", "orange", "blue")
for(i in seq(Fro_bothStats_FDRsigMEnm)){
  pdf(Fro_bothStats_FDRsigCor_MMvsGSnms[i], width = 12, height = 7)
  par(mfrow = c(1,2))
  #P1: all cases combined and controls;
  verboseScatterplot(x = abs(Fro_CorbothTrt_FDRSig_Modclus_MM[[Fro_bothStats_FDRsigMEnm[i]]]),
                     y = abs(Fro_FDRSigCor_Mod_GSolClusNCont[[Fro_bothStats_FDRsigMEnm[i]]]),
                     xlab = paste("Module Membership in", Fro_bothStats_FDRsigMEnm[i], "module"),
                     ylab = paste("Gene significance in all cases combined and controls (Frontal brain)"),
                     main = paste("Module membership vs. gene significance \n"),
                     cex.main = 1, cex.lab = 0.5, cex.axis = 1, pch = 20, cex = 0.5, 
                     col = newFro_bothStats_FDRsigMEnm[i], abline = T, abline.color = newFro_bothStats_FDRsigMEnm[i])
  text(x = abs(Fro_CorbothTrt_FDRSig_Modclus_MM[[Fro_bothStats_FDRsigMEnm[i]]]),
       y = abs(Fro_FDRSigCor_Mod_GSolClusNCont[[Fro_bothStats_FDRsigMEnm[i]]]),
       labels = Fro_alCasesVsCont_WGCNA_GeneInf3[Fro_CorbothTrt_FDRSig_Mod_ClusIDs[[Fro_bothStats_FDRsigMEnm[i]]], "geneName"], 
       cex = 0.5, pos = 3, font = 2, col = "mediumvioletred")
  abline(h = 0.6, v = 0.8, col = "chartreuse3", lty = 2)
  #P2: each cases grouped individually and controls;
  verboseScatterplot(x = abs(Fro_CorbothTrt_FDRSig_Modclus_MM[[Fro_bothStats_FDRsigMEnm[i]]]),
                     y = abs(Fro_FDRSigCor_Mod_GSeachClusNCont[[Fro_bothStats_FDRsigMEnm[i]]]),
                     xlab = paste("Module Membership in", Fro_bothStats_FDRsigMEnm[i], "module"),
                     ylab = paste("Gene significance in each cases grouped separately and controls (Frontal brain)"),
                     main = paste("Module membership vs. gene significance \n"),
                     cex.main = 1, cex.lab = 0.5, cex.axis = 1, pch = 20, cex = 0.5, 
                     col = newFro_bothStats_FDRsigMEnm[i], abline = T, abline.color = newFro_bothStats_FDRsigMEnm[i])
  text(x = abs(Fro_CorbothTrt_FDRSig_Modclus_MM[[Fro_bothStats_FDRsigMEnm[i]]]),
       y = abs(Fro_FDRSigCor_Mod_GSeachClusNCont[[Fro_bothStats_FDRsigMEnm[i]]]),
       labels = Fro_alCasesVsCont_WGCNA_GeneInf3[Fro_CorbothTrt_FDRSig_Mod_ClusIDs[[Fro_bothStats_FDRsigMEnm[i]]], "geneName"], 
       cex = 0.5, pos = 3, font = 2, col = "mediumvioletred")
  abline(h = 0.6, v = 0.8, col = "chartreuse3", lty = 2)
  dev.off()
}


# [F] identify annotation and the modules associated with the key FTD genes ##########################################################################################################; 
KeyFTDGenes <- c("C9orf72", "CHCHD10", "CHMP2B", "FUS", "GRN", "HNRNPA1", "HNRNPA2B1", "MAPT", "OPTN", "PRKAR1B",  "PSEN1", "SIGMAR1",  "SQSTM1", "TARDBP", "TREM2", "UBQLN2", "VCP")
length(KeyFTDGenes)
# [1] 17;

# 1. extract rows containing key ftd genes;
Fro_KeyFTDGenes_asocModInf_ls <- lapply(seq(KeyFTDGenes), function(x) Fro_alCasesVsCont_WGCNA_GeneInf3[which(toupper(Fro_alCasesVsCont_WGCNA_GeneInf3$geneName) == toupper(KeyFTDGenes)[x]), 1:23])
Fro_KeyFTDGenes_asocModInf_df <- do.call(rbind, lapply(Fro_KeyFTDGenes_asocModInf_ls, data.frame, stringsAsFactors=FALSE))
dim(Fro_KeyFTDGenes_asocModInf_df)
# [1] 38 23;

# 2. distribution of key genes in different modules;
table(with(Fro_KeyFTDGenes_asocModInf_df, paste(geneName, MrgdMECols, sep = "_")))
# C9orf72_blue  CHCHD10_lightgreen  CHMP2B_greenyellow          FUS_purple            GRN_pink       HNRNPA1_green      HNRNPA1_purple     HNRNPA2B1_green    HNRNPA2B1_orange 
# 2                   1                   1                   1                   1                   1                   1                   3                   1 
# HNRNPA2B1_purple       HNRNPA2B1_red HNRNPA2B1_turquoise           MAPT_blue         MAPT_orange      MAPT_turquoise          OPTN_green     OPTN_lightgreen            OPTN_red 
# 1                   2                   1                   1                   1                   1                   1                   1                   1 
# OPTN_turquoise      PRKAR1B_orange   PRKAR1B_turquoise         PSEN1_brown       SIGMAR1_black  SQSTM1_darkmagenta    SQSTM1_royalblue    SQSTM1_turquoise          TARDBP_red 
# 2                   2                   1                   1                   1                   1                   1                   3                   1 
# TREM2_pink    UBQLN2_turquoise           VCP_black 
# 1                   1                   1

# 2.1 list modules containing key genes;
Fro_KeyFTDGenes_ModNm_ls <- vector("list", length = length(unique(Fro_KeyFTDGenes_asocModInf_df$MrgdMECols)))
names(Fro_KeyFTDGenes_ModNm_ls) <- unique(Fro_KeyFTDGenes_asocModInf_df$MrgdMECols)
length(Fro_KeyFTDGenes_ModNm_ls)
# [1] 13;

# @distribute key genes corresponding to each module;
for(i in seq(Fro_KeyFTDGenes_ModNm_ls)){
  Fro_KeyFTDGenes_ModNm_ls[[i]] <- paste(lapply(seq(Fro_KeyFTDGenes_ModNm_ls), function(x) 
    unique(Fro_KeyFTDGenes_asocModInf_df$geneName[which(Fro_KeyFTDGenes_asocModInf_df$MrgdMECols == names(Fro_KeyFTDGenes_ModNm_ls)[x])]))[[i]], collapse = "/")
}

unlist(Fro_KeyFTDGenes_ModNm_ls)
# blue;            lightgreen;     greenyellow; purple; 
# "C9orf72/MAPT";  "CHCHD10/OPTN"; "CHMP2B"     "FUS/HNRNPA1/HNRNPA2B1";
# pink;         green;                   orange;                   red; 
# "GRN/TREM2"; "HNRNPA1/HNRNPA2B1/OPTN"; "HNRNPA2B1/MAPT/PRKAR1B"; "HNRNPA2B1/OPTN/TARDBP";
# turquoise;                                    brown;   black;         darkmagenta; 
# "HNRNPA2B1/MAPT/OPTN/PRKAR1B/SQSTM1/UBQLN2"; "PSEN1"; "SIGMAR1/VCP";  "SQSTM1";
# royalblue; 
# "SQSTM1";

# 2.2 crosscheck;
unique(unlist(lapply(seq(Fro_KeyFTDGenes_asocModInf_ls), function(x)
  identical(Fro_KeyFTDGenes_asocModInf_ls[[x]], Fro_KeyFTDGenes_asocModInf_df[rownames(Fro_KeyFTDGenes_asocModInf_ls[[x]]),]))))
# [1] TRUE;
unique(unlist(lapply(seq(length(KeyFTDGenes)), function(x)
  identical(Fro_KeyFTDGenes_asocModInf_ls[[x]], Fro_alCasesVsCont_WGCNA_GeneInf3[which(Fro_alCasesVsCont_WGCNA_GeneInf3$geneName == KeyFTDGenes[x]), 1:23] ) )))
# [1] TRUE;
identical(Fro_KeyFTDGenes_asocModInf_df, 
          Fro_alCasesVsCont_WGCNA_GeneInf3[unlist(lapply(seq(KeyFTDGenes), function(x) which(toupper(Fro_alCasesVsCont_WGCNA_GeneInf3$geneName) == toupper(KeyFTDGenes)[x]))), 1:23])
# [1] TRUE;


# [G] Intramodular analysis: key FTD-genes with high GS and MM #######################################################################;

# 1. names of the modules containing key FTD genes;
Fro_KeyFTDGenes_asocModNms <- unique(Fro_KeyFTDGenes_asocModInf_df$MrgdMECols)
length(Fro_KeyFTDGenes_asocModNms)
# [1] 13;

# 1) Extract cluster IDs relevant to each modules _____________________________________________________________________________________;  
Fro_KeyFTDGenes_asocMod_ClusIDs <- Extract_ClusIDs(ModNames = Fro_KeyFTDGenes_asocModNms, ModAnnoDF = Fro_alCasesVsCont_WGCNA_GeneInf3)
length(Fro_KeyFTDGenes_asocMod_ClusIDs)  
# [1] 13;
unlist(lapply(Fro_KeyFTDGenes_asocMod_ClusIDs, length))
# blue  lightgreen greenyellow      purple        pink       green      orange         red   turquoise       brown       black darkmagenta   royalblue 
# 5475         187         316         326         363         651         798         566        6393         918         480         101         178;

# 2) Extract MMvalues for the cluster-IDs in to each modules _____________________________________________________________________________; 
Fro_KeyFTDGenes_asocMod_MM <- Extract_MM(ModNames = names(Fro_KeyFTDGenes_asocMod_ClusIDs), ModAnnoDF = Fro_alCasesVsCont_WGCNA_GeneInf3)
length(Fro_KeyFTDGenes_asocMod_MM)
# [1] 13;

# @crosscheck;
identical(names(Fro_KeyFTDGenes_asocMod_MM), Fro_KeyFTDGenes_asocModNms)
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_KeyFTDGenes_asocMod_MM)), function(x) 
  identical(Fro_KeyFTDGenes_asocMod_MM[[x]], 
            Fro_alCasesVsCont_WGCNA_GeneInf3[which(Fro_alCasesVsCont_WGCNA_GeneInf3$MrgdMECols == names(Fro_KeyFTDGenes_asocMod_MM)[x]),
                                             match(paste("MM", names(Fro_KeyFTDGenes_asocMod_MM)[x], sep = ""), colnames(Fro_alCasesVsCont_WGCNA_GeneInf3))] ) )))
# [1] TRUE;

# @crosscheck using clusID names;
identical(names(Fro_KeyFTDGenes_asocMod_MM), names(Fro_KeyFTDGenes_asocMod_ClusIDs))
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_KeyFTDGenes_asocMod_MM)), function(x)
  identical(Fro_KeyFTDGenes_asocMod_MM[[x]], Fro_alCasesVsCont_WGCNA_GeneInf3[Fro_KeyFTDGenes_asocMod_ClusIDs[[x]],
                                                                              match(paste("MM", names(Fro_KeyFTDGenes_asocMod_MM)[x], sep = ""), colnames(Fro_alCasesVsCont_WGCNA_GeneInf3))] ) )))
# [1] TRUE;

# @crosscheck;
identical(names(Fro_KeyFTDGenes_asocMod_ClusIDs), Fro_KeyFTDGenes_asocModNms)
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_KeyFTDGenes_asocMod_ClusIDs)), function(x)
  identical(Fro_KeyFTDGenes_asocMod_ClusIDs[[x]], as.character(with(Fro_alCasesVsCont_WGCNA_GeneInf3, Clus_ID[which(MrgdMECols == names(Fro_KeyFTDGenes_asocMod_ClusIDs)[x])]) )) )))
# [1] TRUE;


# 3) Extract GSvalues for significant modules _______________________________________________________________________________________________; 

# 3.1 all case combined and cont - extract their GS values;
Fro_KeyFTDGenes_asocMod_GSolClusNCont <- Extract_GS_both(ModNames = names(Fro_KeyFTDGenes_asocMod_ClusIDs), ModAnnoDF = Fro_alCasesVsCont_WGCNA_GeneInf3, TraitColNm = "GSStat_CaseCont")
length(Fro_KeyFTDGenes_asocMod_GSolClusNCont)
# [1] 13;  

# @crosscheck;
identical(names(Fro_KeyFTDGenes_asocMod_GSolClusNCont), Fro_KeyFTDGenes_asocModNms)
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_KeyFTDGenes_asocMod_GSolClusNCont)), function(x)
  identical(Fro_KeyFTDGenes_asocMod_GSolClusNCont[[x]], Fro_alCasesVsCont_WGCNA_GeneInf3[which(Fro_alCasesVsCont_WGCNA_GeneInf3$MrgdMECols == names(Fro_KeyFTDGenes_asocMod_GSolClusNCont)[x]), "GSStat_CaseCont"]) )))
# [1] TRUE;

# @crosscheck using the identified clus-ids;
identical(names(Fro_KeyFTDGenes_asocMod_GSolClusNCont), names(Fro_KeyFTDGenes_asocMod_ClusIDs))
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_KeyFTDGenes_asocMod_GSolClusNCont)), function(x)
  identical(Fro_KeyFTDGenes_asocMod_GSolClusNCont[[x]], Fro_alCasesVsCont_WGCNA_GeneInf3[Fro_KeyFTDGenes_asocMod_ClusIDs[[x]], "GSStat_CaseCont"]) )))
# [1] TRUE;


# 3.2 each case and cont - extract their GS values;
Fro_KeyFTDGenes_asocMod_GSeachClusNCont <- Extract_GS_both(ModNames = names(Fro_KeyFTDGenes_asocMod_ClusIDs), ModAnnoDF = Fro_alCasesVsCont_WGCNA_GeneInf3, TraitColNm = "GSStat_MutCat")
length(Fro_KeyFTDGenes_asocMod_GSeachClusNCont)
# [1] 13;  

# @crosscheck;
identical(names(Fro_KeyFTDGenes_asocMod_GSeachClusNCont), Fro_KeyFTDGenes_asocModNms)
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_KeyFTDGenes_asocMod_GSeachClusNCont)), function(x)
  identical(Fro_KeyFTDGenes_asocMod_GSeachClusNCont[[x]], Fro_alCasesVsCont_WGCNA_GeneInf3[which(Fro_alCasesVsCont_WGCNA_GeneInf3$MrgdMECols == names(Fro_KeyFTDGenes_asocMod_GSeachClusNCont)[x]), "GSStat_MutCat"]) )))
# [1] TRUE;

# @crosscheck using the identified clus-ids;
identical(names(Fro_KeyFTDGenes_asocMod_GSeachClusNCont), names(Fro_KeyFTDGenes_asocMod_ClusIDs))
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_KeyFTDGenes_asocMod_GSeachClusNCont)), function(x) identical(Fro_KeyFTDGenes_asocMod_GSeachClusNCont[[x]],
                                                                                                 Fro_alCasesVsCont_WGCNA_GeneInf3[Fro_KeyFTDGenes_asocMod_ClusIDs[[x]], "GSStat_MutCat"]) )))
# [1] TRUE;



# [H] Data prep for scatter plot generation ############################################################################;

# 1. names for the pdf files;
Fro_KeyFTDGenes_asocMod_MMvsGSnms <- unlist(lapply(seq(Fro_KeyFTDGenes_asocModNms), 
                                                   function(x) paste("Fro_KeyFTDGenes_asocMod_",  Fro_KeyFTDGenes_asocModNms[x], "_MMvsGS_min50Cpt1.pdf", sep = "")))
length(Fro_KeyFTDGenes_asocMod_MMvsGSnms)
# [1] 13;

# 2. use the function (Fro_KeyFTDGenes_asocMod_AnnoDF) that extract annodf sepearately for each mod in a ls;

# @execute the function;
Fro_KeyFTDGenes_AnnoDF_perModls <- Fro_KeyFTDGenes_asocMod_AnnoDF(ModNames = Fro_KeyFTDGenes_asocModNms, ModAnnoDF = Fro_alCasesVsCont_WGCNA_GeneInf3)
unlist(lapply(Fro_KeyFTDGenes_AnnoDF_perModls, nrow))
# blue  lightgreen greenyellow      purple        pink       green      orange         red   turquoise       brown       black darkmagenta   royalblue; 
# 5475         187         316         326         363         651         798         566        6393         918         480         101         178;

# @crosscheck;
identical(names(Fro_KeyFTDGenes_AnnoDF_perModls), Fro_KeyFTDGenes_asocModNms)
# [1] TRUE;
identical(Fro_KeyFTDGenes_asocModNms,
          unlist(lapply(seq(length(Fro_KeyFTDGenes_AnnoDF_perModls)), function(x) unique(Fro_KeyFTDGenes_AnnoDF_perModls[[x]]$MrgdMECols))))
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_KeyFTDGenes_AnnoDF_perModls)), function(x)
  identical(Fro_KeyFTDGenes_AnnoDF_perModls[[x]], Fro_alCasesVsCont_WGCNA_GeneInf3[which(Fro_alCasesVsCont_WGCNA_GeneInf3$MrgdMECols == names(Fro_KeyFTDGenes_AnnoDF_perModls)[x]), 
                                                                                   colnames(Fro_KeyFTDGenes_AnnoDF_perModls[[x]])]) )))
# [1] TRUE;

# 3. use the function (maxXYthresh_perMod) to get the range of x (MM) and y (GS) values;
# trial <- maxXYthresh_perMod(ModAnnoDF = Fro_KeyFTDGenes_AnnoDF_perModls, GStrt1 = "GSStat_CaseCont", GStrt2 = "GSStat_MutCat")

# @execute the function;
Fro_KeyFTDGenes_maxXYthresh <- maxXYthresh_perMod(ModAnnoDF = Fro_KeyFTDGenes_AnnoDF_perModls, GStrt1 = "GSStat_CaseCont", GStrt2 = "GSStat_MutCat")
length(Fro_KeyFTDGenes_maxXYthresh)
# [1] 13;

# -crosscheck;
identical(names(Fro_KeyFTDGenes_maxXYthresh), names(Fro_KeyFTDGenes_AnnoDF_perModls))
# [1] TRUE;
unique(unlist(lapply(seq(length(Fro_KeyFTDGenes_maxXYthresh)), function(x)
  identical(as.numeric(Fro_KeyFTDGenes_maxXYthresh[[x]]), 
            lapply(seq(length(Fro_KeyFTDGenes_maxXYthresh)), function(x)
              c(as.numeric(substr(max(abs(with(Fro_alCasesVsCont_WGCNA_GeneInf3, GSStat_CaseCont[MrgdMECols == names(Fro_KeyFTDGenes_maxXYthresh)[x]]))), 1, 3)),
                as.numeric(substr(max(abs(with(Fro_alCasesVsCont_WGCNA_GeneInf3, GSStat_MutCat[MrgdMECols == names(Fro_KeyFTDGenes_maxXYthresh)[x]]))), 1, 3)),
                as.numeric(substr(max(abs(Fro_alCasesVsCont_WGCNA_GeneInf3[,match(paste("MM", names(Fro_KeyFTDGenes_maxXYthresh)[x], sep = ""), colnames(Fro_alCasesVsCont_WGCNA_GeneInf3))])), 1, 3)) ))[[x]]) )))
# [1] TRUE;

# 4. use the function (maxXYthresh_perMod) that labels only the data points corresponding to the key genes and the hub genes only (***The combined function);

# 4.1 Trait: All cases combined and controls __________________________________________________;
Fro_allCaseNCont_GenLab_alCaseCont <- NonHubKeyGens_convBlnk(KeyGens_AnnDFls = Fro_KeyFTDGenes_AnnoDF_perModls, KeyGens_XYthres_ls = Fro_KeyFTDGenes_maxXYthresh, GStrt = "GSStat_CaseCont", keyGeNms = KeyFTDGenes)
length(Fro_allCaseNCont_GenLab_alCaseCont)
# [1] 13;

# -crosscheck;
identical(Fro_allCaseNCont_GenLab_alCaseCont, Fro_KeyFTDGenes_AnnoDF_perModls)
# [1] FALSE; 
# NOTE: correct meaning that uninteresting geneNames were labelled blanks unlike the origninal file;

# --crosscheck the two list DFs but excluding the geneName column;
unique(unlist(lapply(seq(length(Fro_KeyFTDGenes_AnnoDF_perModls)), function(x) 
  identical(Fro_KeyFTDGenes_AnnoDF_perModls[[x]][-12], Fro_allCaseNCont_GenLab_alCaseCont[[x]][-12]) )))
# [1] TRUE;

# @test if the function returned the key FTD genes;
length(unlist(lapply(lapply(seq(Fro_allCaseNCont_GenLab_alCaseCont), function(x) toupper(Fro_allCaseNCont_GenLab_alCaseCont[[x]]$geneName)), intersect, toupper(KeyFTDGenes))))
# [1] 30;

# 4.2 Trait: each individual cases grouped separately and the controls ________________________;
Fro_allCaseNCont_GenLab_eachCaseCont <- NonHubKeyGens_convBlnk(KeyGens_AnnDFls = Fro_KeyFTDGenes_AnnoDF_perModls, KeyGens_XYthres_ls = Fro_KeyFTDGenes_maxXYthresh, GStrt = "GSStat_MutCat", keyGeNms = KeyFTDGenes)
length(Fro_allCaseNCont_GenLab_eachCaseCont)
# [1] 13;

# -crosscheck;
identical(Fro_allCaseNCont_GenLab_eachCaseCont, Fro_KeyFTDGenes_AnnoDF_perModls)
# [1] FALSE; 
# NOTE: correct meaning that uninteresting geneNames were labelled blanks unlike the origninal file;

# --crosscheck the two list DFs but excluding the geneName column;
unique(unlist(lapply(seq(length(Fro_KeyFTDGenes_AnnoDF_perModls)), function(x) 
  identical(Fro_KeyFTDGenes_AnnoDF_perModls[[x]][-12], Fro_allCaseNCont_GenLab_eachCaseCont[[x]][-12]) )))
# [1] TRUE;

# @test if the function returned the key FTD genes;
length(unlist(lapply(lapply(seq(Fro_allCaseNCont_GenLab_eachCaseCont), function(x) toupper(Fro_allCaseNCont_GenLab_eachCaseCont[[x]]$geneName)), intersect, toupper(KeyFTDGenes))))
# [1] 30;

# !light yellow dots are invisible, remedy that!;
# Fro_KeyFTDGenes_asocModNms_cols <- ifelse(Fro_KeyFTDGenes_asocModNms == "lightyellow", "goldenrod", Fro_KeyFTDGenes_asocModNms)
Fro_KeyFTDGenes_asocModNms_cols <- Fro_KeyFTDGenes_asocModNms

# 5. save the plots;
for(i in seq(Fro_KeyFTDGenes_asocModNms)){
  pdf(Fro_KeyFTDGenes_asocMod_MMvsGSnms[i], width = 12, height = 7)
  par(mfrow = c(1,2))
  #P1: all cases combined and controls;
  verboseScatterplot(x = abs(Fro_KeyFTDGenes_asocMod_MM[[Fro_KeyFTDGenes_asocModNms[i]]]),
                     y = abs(Fro_KeyFTDGenes_asocMod_GSolClusNCont[[Fro_KeyFTDGenes_asocModNms[i]]]),
                     xlab = paste("Module Membership in", Fro_KeyFTDGenes_asocModNms[i], "module"),
                     ylab = paste("Gene significance in all cases combined and controls"),
                     main = paste("Frontal brain: module membership vs. gene significance \n"),
                     cex.main = 1, cex.lab = 0.5, cex.axis = 1, pch = 20, cex = 0.5,
                     col = Fro_KeyFTDGenes_asocModNms_cols[i], abline = T, abline.color = Fro_KeyFTDGenes_asocModNms_cols[i])
  text(x = abs(Fro_KeyFTDGenes_asocMod_MM[[Fro_KeyFTDGenes_asocModNms[i]]]),
       y = abs(Fro_KeyFTDGenes_asocMod_GSolClusNCont[[Fro_KeyFTDGenes_asocModNms[i]]]),
       labels = Fro_allCaseNCont_GenLab_alCaseCont[[Fro_KeyFTDGenes_asocModNms[i]]][, "geneName"], cex = 0.5, pos = 3, font = 2, col = "mediumvioletred")
  abline(h = Fro_KeyFTDGenes_maxXYthresh[[Fro_KeyFTDGenes_asocModNms[i]]][["GSStat_CaseCont"]]-0.1,
         v = Fro_KeyFTDGenes_maxXYthresh[[Fro_KeyFTDGenes_asocModNms[i]]][[grep("MM", names(Fro_KeyFTDGenes_maxXYthresh[[Fro_KeyFTDGenes_asocModNms[i]]]))]]-0.1, col = "chartreuse3", lty = 2)
  #P2: each cases grouped individually and controls;
  verboseScatterplot(x = abs(Fro_KeyFTDGenes_asocMod_MM[[Fro_KeyFTDGenes_asocModNms[i]]]),
                     y = abs(Fro_KeyFTDGenes_asocMod_GSeachClusNCont[[Fro_KeyFTDGenes_asocModNms[i]]]),
                     xlab = paste("Module Membership in", Fro_KeyFTDGenes_asocModNms[i], "module"),
                     ylab = paste("Gene significance in each individual cases grouped separately and controls"),
                     main = paste("Frontal brain: module membership vs. gene significance \n"),
                     cex.main = 1, cex.lab = 0.5, cex.axis = 1, pch = 20, cex = 0.5,
                     col = Fro_KeyFTDGenes_asocModNms_cols[i], abline = T, abline.color = Fro_KeyFTDGenes_asocModNms_cols[i])
  text(x = abs(Fro_KeyFTDGenes_asocMod_MM[[Fro_KeyFTDGenes_asocModNms[i]]]),
       y = abs(Fro_KeyFTDGenes_asocMod_GSeachClusNCont[[Fro_KeyFTDGenes_asocModNms[i]]]),
       labels = Fro_allCaseNCont_GenLab_eachCaseCont[[Fro_KeyFTDGenes_asocModNms[i]]][, "geneName"], cex = 0.5, pos = 3, font = 2, col = "mediumvioletred")
  abline(h = Fro_KeyFTDGenes_maxXYthresh[[Fro_KeyFTDGenes_asocModNms[i]]][["GSStat_MutCat"]]-0.1,
         v = Fro_KeyFTDGenes_maxXYthresh[[Fro_KeyFTDGenes_asocModNms[i]]][[grep("MM", names(Fro_KeyFTDGenes_maxXYthresh[[Fro_KeyFTDGenes_asocModNms[i]]]))]]-0.1, col = "chartreuse3", lty = 2)
  dev.off()
}



# ##############################################################################################################;
# III. Save this workspace;
# ______________________________________________________________________________________________________________;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc")
save.image("WGCNA_S3_FronSamps_AllCasesNConts_CAGEseq.RData")
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc/WGCNA_S3_FronSamps_AllCasesNConts_CAGEseq.RData")
# END __________________________________________________________________________________________________________________________________________________________________________________________;
