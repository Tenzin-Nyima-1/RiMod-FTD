# Date: 15 May, 2017;
# Project: RIMOD-FTD;
# Aim: differential gene expression (DGE) analysis of frontal samples: MAPT cases vs controls;
# NOTE: the annotation of the count-table was carried using RIKEN's CAT gtf file as the gene model;
# NOTE 2: the same script could be adapted for the count table from: all 7 brain regions; and, temporal brain; 
# Workflow (main steps);
# Step 1: Upload the workspace of pre-DGE;
# Step 2: Data distribution check;
# Step 3: DGE-analysis;


# source("http://bioconductor.org/biocLite.R")
# biocLite()
options(stringsAsFactors=FALSE)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(S4Vectors)
library(edgeR)
library(stringi)
library(stringr)
library(plotrix)
library(sva)
library(RColorBrewer)
library(pheatmap)
library(muStat)
library(plot3D)


# ###################################################################################################;
# I. Upload the workspace that performed pre-DGE analysis: crosschecks, and organized frontal samples;
# ###################################################################################################;
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc_fro/S6_FronSamps_preDGEanaly_crschk_orgnz.RData")

# 1. the ordered annotated table;
dim(Ann_FronOnly_peakGR3k_ordFinal)
# [1] 31000    66;

# 2. library info for the frontal samples only;
dim(updt_GoogleSht_FronOnly_mtchd)
# [1] 58 22;



# #######################################;
# II. Data preparation for edgeR analysis;
# _______________________________________;


# [A] Extract the count-table;
Fron_MAPTvsCont_counts <- Ann_FronOnly_peakGR3k_ordFinal[, grep("control|MAPT", colnames(Ann_FronOnly_peakGR3k_ordFinal), ignore.case = T)]
dim(Fron_MAPTvsCont_counts)
# [1] 31000    28;

# @crosscheck: if we're missing or including wrong samples;
Mut_states <- c("control", "MAPT")
identical(updt_GoogleSht_FronOnly_mtchd$GIVENSAMPLENAME[unlist(lapply(seq(Mut_states), function(x) grep(Mut_states[x], updt_GoogleSht_FronOnly_mtchd$MUTATED.GENE, ignore.case = TRUE)))],
          stri_split_fixed(stri_split_fixed(colnames(Fron_MAPTvsCont_counts), "sample_", simplify = TRUE)[,2], "_no_chrM_", simplify = TRUE)[,1])
# [1] TRUE; NOTE: all the cont & MAPT samples were correctly selected;



# [B] Define groups for each study category; 
Fron_MAPTvsCont_Gps <- str_sub(stri_split_fixed(colnames(Fron_MAPTvsCont_counts), "_no_chrM_", simplify = TRUE)[,2], 1, 4)
Fron_MAPTvsCont_GpsOrd <- factor(Fron_MAPTvsCont_Gps, levels = unique(Fron_MAPTvsCont_Gps))
length(Fron_MAPTvsCont_GpsOrd)
# [1] 28;

# @another layer of crosscheck;
identical(c(table(updt_GoogleSht_FronOnly_mtchd$MUTATED.GENE[grep("cont|MAPT", updt_GoogleSht_FronOnly_mtchd$MUTATED.GENE)])), 
          c(table(stri_split_fixed(colnames(Fron_MAPTvsCont_counts), "_no_chrM_", simplify = TRUE)[,2])))
# [1] TRUE; NOTE: this confirms that the correct #cont/mapt samples are present;


# [C] Gene info;
Fron_MAPTvsCont_genes <- Ann_FronOnly_peakGR3k_ordFinal[,-(grep("sample", colnames(Ann_FronOnly_peakGR3k_ordFinal)))]
dim(Fron_MAPTvsCont_genes)
# [1] 31000    17;

# 1. check if the duplicated gene columns are identical;
identical(rownames(Fron_MAPTvsCont_genes), as.character(Fron_MAPTvsCont_genes$Clus_ID))
# [1] TRUE;
identical(as.character(Fron_MAPTvsCont_genes$seqnames), stri_split_fixed(Fron_MAPTvsCont_genes$Clus_ID, "_", simplify = T)[,1])
# [1] TRUE;
identical(as.character(Fron_MAPTvsCont_genes$geneChr), stri_split_fixed(Fron_MAPTvsCont_genes$Clus_ID, "_", simplify = T)[,1])
# [1] TRUE;
identical(as.character(Fron_MAPTvsCont_genes$strand), stri_split_fixed(Fron_MAPTvsCont_genes$Clus_ID, "_", simplify = T)[,4])
# [1] TRUE;
identical(as.character(Fron_MAPTvsCont_genes$Strand), stri_split_fixed(Fron_MAPTvsCont_genes$Clus_ID, "_", simplify = T)[,4])
# [1] TRUE;
identical(as.character(Fron_MAPTvsCont_genes$geneStrand), stri_split_fixed(Fron_MAPTvsCont_genes$Clus_ID, "_", simplify = T)[,4])
# [1] TRUE;
# NOTE: So, keep only one strand info; and, discard-seqnames,score;  

# 2. the ordered genes file;
Fron_MAPTvsCont_GenesOrd <- Fron_MAPTvsCont_genes[,c("Clus_ID", "seqnames", "start", "end", "strand", "width", "annotation", "distanceToTSS", "geneId", "transcriptId", "geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand")]
dim(Fron_MAPTvsCont_GenesOrd)
# [1] 31000    15;

# 3. read the F6 geneinfo(from chang on 12 Jan 2017) file, in order to get the alternative converted ids for each gene-id;
F6_CATgene_info <- read.delim("/home/student-ag-javier/Documents/Fantom_6/Datas/F6_CAT_GTFfile/F6_CAT_GeneInf_RecvOn12Jan17/F6_CAT.gene.info.tsv",  header = TRUE, sep = "\t", quote = "")
dim(F6_CATgene_info)
#[1] 124047     41;

# 4. Get additional gene info from the above CAT file;

# i) cbind geneInfo using F6 file;
which(duplicated(F6_CATgene_info$geneID))
# integer(0) # NOTE: there isn't any duplicated gene-ids so its save to match with them;
Fron_MAPTvsCont_GenesOrd_F6 <- cbind(Fron_MAPTvsCont_GenesOrd, 
                                     F6_CATgene_info[match(Fron_MAPTvsCont_GenesOrd$geneId, F6_CATgene_info$geneID), 
                                                     c("geneID", "geneName", "geneType", "CAT_geneClass", "CAT_DHS_type", "trnscptID", "HGNC_name", "HGNC_symbol", "entrez_ID", "refseq_ID")])
dim(Fron_MAPTvsCont_GenesOrd_F6)
# [1] 31000    25;

# ii) crosscheck;

# @right gene-ids matched?;
identical(Fron_MAPTvsCont_GenesOrd_F6$geneId, Fron_MAPTvsCont_GenesOrd_F6$geneID)
# [1] TRUE;

# @right TX-ids?;
Fron_MAPTvsCont_GenesOrd_F6_TXchk <- unlist(lapply(lapply(seq(nrow(Fron_MAPTvsCont_GenesOrd_F6)), 
                                                          function(x) grep(Fron_MAPTvsCont_GenesOrd_F6$transcriptId[x], Fron_MAPTvsCont_GenesOrd_F6$trnscptID[x])), length))
length(Fron_MAPTvsCont_GenesOrd_F6_TXchk)
# [1] 31000;
unique(Fron_MAPTvsCont_GenesOrd_F6_TXchk)
# [1] 1;
# NOTE: both gene-ids and the TX ids matched and confirmed that they are in the same order, so remove the extra columns from the data-frame;

# iii) Cleared crosscheck, so remove the extra columns;
Fron_MAPTvsCont_F6Genes_fin <- Fron_MAPTvsCont_GenesOrd_F6[,-(match(c("geneID", "trnscptID"), colnames(Fron_MAPTvsCont_GenesOrd_F6)))]
dim(Fron_MAPTvsCont_F6Genes_fin)
# [1] 31000    23;
range(Fron_MAPTvsCont_F6Genes_fin$width)
# [1]    2 2915;



# ###########################;
# III. Perform edgeR analysis;
# ___________________________;

# 1. create a dgelist object;
Fron_MAPTvsCont_DGls <- DGEList(counts = Fron_MAPTvsCont_counts, group = Fron_MAPTvsCont_GpsOrd, genes = Fron_MAPTvsCont_F6Genes_fin)
dim(Fron_MAPTvsCont_DGls)
# [1] 31000    28;

# @a crosscheck;
identical(rownames(Fron_MAPTvsCont_DGls), as.character(Fron_MAPTvsCont_DGls$genes$Clus_ID))
# [1] TRUE;


# 2. filter out lowly expressed TXs;

# i) setting 1cpm threshold to estimate raw count range for that value;
range((1 * Fron_MAPTvsCont_DGls$samples$lib.size)/1e6)
# [1]  6.76344 32.21440;

# ii) get the #sample in both cont and case categories;
table(Fron_MAPTvsCont_DGls$samples$group)
# cont MAPT 
# 16   12

# iii) since the lowest #sample is 12, set the filter for TX with <1cpm in that many samples;
Fron_MAPTvsCont_1in12_keepIDs <- rowSums(cpm(Fron_MAPTvsCont_DGls) > 1) >= 12
summary(Fron_MAPTvsCont_1in12_keepIDs)
# Mode   FALSE    TRUE    NA's 
# logical    7578   23422       0

# iv) now select only those TXs that passed the threshold;
Fron_MAPTvsCont_DGls_1in12 <- Fron_MAPTvsCont_DGls[Fron_MAPTvsCont_1in12_keepIDs, , keep.lib.sizes=FALSE]
dim(Fron_MAPTvsCont_DGls_1in12)
# [1] 23422    28;

# @crosscheck if the library sizes were kept similar;
all(Fron_MAPTvsCont_DGls_1in12$samples$lib.size == Fron_MAPTvsCont_DGls$samples$lib.size)
# [1] FALSE; NOTE: correct, meaning that lib-size were adjusted after the filtering;


# 3. normalize using RLE with filtered transcripts only;
Fron_MAPTvsCont_DGls_1in12_nrmd <- calcNormFactors(Fron_MAPTvsCont_DGls_1in12, method = "RLE")
dim(Fron_MAPTvsCont_DGls_1in12_nrmd)
# [1] 23422    28;


# 4. prepare the data for mds plot;  
Fron_MAPTvsCont_mds_df <- data.frame(SampNms = stri_split_fixed(stri_split_fixed(colnames(Fron_MAPTvsCont_DGls_1in12_nrmd), "sample_", simplify = T)[,2], "_no_chrM_", simplify = T)[,1],
                                     Mutation = stri_split_fixed(colnames(Fron_MAPTvsCont_DGls_1in12_nrmd), "_no_chrM_", simplify = T)[,2], 
                                     Samp_Mut = NA, 
                                     Category = Fron_MAPTvsCont_DGls_1in12_nrmd$samples$group,
                                     Col = NA,
                                     Pch = 8, stringsAsFactors = FALSE)
# @rownames;
rownames(Fron_MAPTvsCont_mds_df) <- colnames(Fron_MAPTvsCont_DGls_1in12_nrmd)
# @Samp_Mut;
Fron_MAPTvsCont_mds_df$Samp_Mut <- paste(str_sub(Fron_MAPTvsCont_mds_df$SampNms, 1, -5), str_sub(Fron_MAPTvsCont_mds_df$Mutation, 1, 2), sep = "_")
# @Colors;
Fron_MAPTvsCont_mds_df$Col <- ifelse(Fron_MAPTvsCont_mds_df$Category == "cont", "forestgreen", "red")

# @legend;
Fron_MAPTvsCont_mds_df_lgnd <- Fron_MAPTvsCont_mds_df[!duplicated(Fron_MAPTvsCont_mds_df$Category),]

# 6. save the MDS plot:#A] method = lfc-pw; #B] method = lfc-com; #C] method = bcv; #D] PCA: 5k genes;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Frontal_RikAnno/ana_DEG/MAPT_vs_Cont")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;
#A] method = lfc - pairwise ____________;
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;
# top = default 500;
Fron_MAPTvsCont_MDS_pairWlfc <- plotMDS(Fron_MAPTvsCont_DGls_1in12_nrmd, gene.selection = "pairwise", method = "logFC")

# plot 1: to check for ourselves;
pdf("lab_Fron_MAPTvsCont_MDS_pairWlfc.pdf")
plotMDS(Fron_MAPTvsCont_MDS_pairWlfc, col = Fron_MAPTvsCont_mds_df$Col, labels = Fron_MAPTvsCont_mds_df$Samp_Mut, cex = 0.5)
mtext(text = "MDS plot: MAPT vs control samples from frontal brain", adj = 0, side = 3, cex = 1, line = 0.5)
legend("topleft", legend = Fron_MAPTvsCont_mds_df_lgnd$Category, cex = 0.7, col = Fron_MAPTvsCont_mds_df_lgnd$Col, pch = Fron_MAPTvsCont_mds_df_lgnd$Pch, pt.bg = Fron_MAPTvsCont_mds_df_lgnd$Col, box.lwd = 0.5)
dev.off()

# plot 2: for presentation purposes;
pdf("ppt_Fron_MAPTvsCont_MDS_pairWlfc.pdf")
plotMDS(Fron_MAPTvsCont_MDS_pairWlfc, col = Fron_MAPTvsCont_mds_df$Col, pch = Fron_MAPTvsCont_mds_df$Pch)
mtext(text = "MDS plot: MAPT vs control samples from frontal brain", adj = 0, side = 3, cex = 1, line = 0.5)
legend("topleft", legend = Fron_MAPTvsCont_mds_df_lgnd$Category, cex = 0.7, col = Fron_MAPTvsCont_mds_df_lgnd$Col, pch = Fron_MAPTvsCont_mds_df_lgnd$Pch, pt.bg = Fron_MAPTvsCont_mds_df_lgnd$Col, box.lwd = 0.5)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;
#B] method = lfc - common ____________;
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;
# top = default 500;
Fron_MAPTvsCont_MDS_commlfc <- plotMDS(Fron_MAPTvsCont_DGls_1in12_nrmd, gene.selection = "common", method = "logFC")

# plot 1: to check for ourselves;
pdf("lab_Fron_MAPTvsCont_MDS_commlfc.pdf")
plotMDS(Fron_MAPTvsCont_MDS_commlfc, col = Fron_MAPTvsCont_mds_df$Col, labels = Fron_MAPTvsCont_mds_df$Samp_Mut, cex = 0.5)
mtext(text = "MDS plot: MAPT vs control samples from frontal brain", adj = 0, side = 3, cex = 1, line = 0.5)
legend("topright", legend = Fron_MAPTvsCont_mds_df_lgnd$Category, cex = 0.7, col = Fron_MAPTvsCont_mds_df_lgnd$Col, pch = Fron_MAPTvsCont_mds_df_lgnd$Pch, pt.bg = Fron_MAPTvsCont_mds_df_lgnd$Col, box.lwd = 0.5)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;
#C] method = bcv _______________________;
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;
# top = default 500;
Fron_MAPTvsCont_MDS_bcv <- plotMDS(Fron_MAPTvsCont_DGls_1in12_nrmd, method = "bcv")

# plot 1: to check for ourselves;
pdf("lab_Temp_MAPTvsCont_MDS_bcv.pdf")
plotMDS(Fron_MAPTvsCont_MDS_bcv, col = Fron_MAPTvsCont_mds_df$Col, labels = Fron_MAPTvsCont_mds_df$Samp_Mut, cex = 0.5)
mtext(text = "MDS plot: MAPT vs control samples from frontal brain", adj = 0, side = 3, cex = 1, line = 0.5)
legend("topright", legend = Fron_MAPTvsCont_mds_df_lgnd$Category, cex = 0.7, col = Fron_MAPTvsCont_mds_df_lgnd$Col, pch = Fron_MAPTvsCont_mds_df_lgnd$Pch, pt.bg = Fron_MAPTvsCont_mds_df_lgnd$Col, box.lwd = 0.5)
dev.off()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;
#D] Generate a PCA plot for the top 4000 genes;
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;

# 1. convert the count into logcpm to avoid count bias due to highly expressed genes;
Fron_MAPTvsCont_count_logcpm <- cpm(Fron_MAPTvsCont_DGls_1in12_nrmd, prior.count=2, log=TRUE)
dim(Fron_MAPTvsCont_count_logcpm)
# [1] 23422    28;

# 2. get the rowVars of the logcpm count and select the expression of the top 5k genes;
Fron_MAPTvsCont_count_logcpm_RV <- rowVars(as.matrix(Fron_MAPTvsCont_count_logcpm))
Fron_MAPTvsCont_count_logcpm_RV_ord <- order(Fron_MAPTvsCont_count_logcpm_RV, decreasing = TRUE)[1:5000]
length(Fron_MAPTvsCont_count_logcpm_RV_ord)
#[1] 5000;

# 3. get the count of the 5k most variant genes;
Fron_MAPTvsCont_count_logcpm_5k <- Fron_MAPTvsCont_count_logcpm[Fron_MAPTvsCont_count_logcpm_RV_ord, ]
dim(Fron_MAPTvsCont_count_logcpm_5k)
# [1] 5000   28;

# 4. get the PCA of the 5k count;
Fron_MAPTvsCont_count_logcpm_5k_PC <- prcomp(t(Fron_MAPTvsCont_count_logcpm_5k), scale. = FALSE)
names(Fron_MAPTvsCont_count_logcpm_5k_PC)
#[1] "sdev"     "rotation" "center"   "scale"    "x";

# get the %variance;
Fron_MAPTvsCont_count_logcpm_5k_PC_perVar <- round((Fron_MAPTvsCont_count_logcpm_5k_PC$sdev^2)*100/sum(Fron_MAPTvsCont_count_logcpm_5k_PC$sdev^2), 2)
# just a check if the both give the same result;
identical(Fron_MAPTvsCont_count_logcpm_5k_PC_perVar, round(100*Fron_MAPTvsCont_count_logcpm_5k_PC$sdev^2/sum(Fron_MAPTvsCont_count_logcpm_5k_PC$sdev^2),2))
#[1] TRUE;

# 5. prepare a data-frame for generating PCA plot;
Fron_MAPTvsCont_count_logcpm_5k_PC_DF <- cbind(Fron_MAPTvsCont_count_logcpm_5k_PC$x[,1:4], Fron_MAPTvsCont_mds_df)

# 6. generate pca plot;
pdf("Fron_MAPTvsCont_RVlcpm_PCA.pdf", width = 8, height = 9)
with(Fron_MAPTvsCont_count_logcpm_5k_PC_DF, text3D(PC1, PC2, PC3, col = Col, theta = 60, phi = 20,
                                                   xlab = paste0("PC1: ",Fron_MAPTvsCont_count_logcpm_5k_PC_perVar[1], "% variance"),
                                                   ylab = paste0("PC2: ",Fron_MAPTvsCont_count_logcpm_5k_PC_perVar[2], "% variance"),
                                                   zlab = paste0("PC3: ",Fron_MAPTvsCont_count_logcpm_5k_PC_perVar[3], "% variance"),
                                                   main = "Frontal-brain: MAPT vs control samples",
                                                   labels = Samp_Mut, cex = 0.6, bty = "g", ticktype = "detailed", d = 2, adj = 0.5, font = 2))
dev.off()


# 7. following implements batch analysis using sva along with DE analysis using edgeR analysis_________________________________________________;

# 7.1 design matrix;
Fron_MAPTvsCont_DsgMat <- Fron_MAPTvsCont_mds_df[,c("SampNms", "Mutation", "Category")]
Fron_MAPTvsCont_DsgMat$Category <- relevel(Fron_MAPTvsCont_DsgMat$Category, ref = "cont")
dim(Fron_MAPTvsCont_DsgMat)
# [1] 28  3;

# @a tiny crosscheck;
identical(Fron_MAPTvsCont_DsgMat$Category, Fron_MAPTvsCont_DGls_1in12_nrmd$samples$group)
# [1] TRUE

# 7.2 model matrixes;

# 7.2.1 full model;
Fron_MAPTvsCont_mainMod <- model.matrix(~Category, data = Fron_MAPTvsCont_DsgMat)
colnames(Fron_MAPTvsCont_mainMod) <- c("Intercept", "MAPT") #instead of "(Intercept) CategoryMAPT";

# @check if the rownames corresponds to the colnames of the DGls object;
identical(rownames(Fron_MAPTvsCont_mainMod), colnames(Fron_MAPTvsCont_DGls_1in12_nrmd))
# [1] TRUE;

# 7.2.2 null model;
Fron_MAPTvsCont_nulMod <- model.matrix(~1, data = Fron_MAPTvsCont_DsgMat)
rownames(Fron_MAPTvsCont_nulMod) <- colnames(Fron_MAPTvsCont_DGls_1in12_nrmd)
colnames(Fron_MAPTvsCont_nulMod) <- "Intercept" #instead of "(Intercept)";

# 7.3 now, run sva;
Fron_MAPTvsCont_sva <- svaseq(cpm(Fron_MAPTvsCont_DGls_1in12_nrmd), Fron_MAPTvsCont_mainMod, Fron_MAPTvsCont_nulMod)
Fron_MAPTvsCont_sva$n.sv
#[1] 5

# @check which methods were implemented, run the following code;
# x_sva_be <- svaseq(cpm(Fron_MAPTvsCont_DGls_1in12_nrmd), Fron_MAPTvsCont_mainMod, Fron_MAPTvsCont_nulMod, method = "irw", numSVmethod = "be")
# identical(x_sva_be, x_sva) #[1] TRUE;
# xx.nsv <- num.sv(cpm(Fron_MAPTvsCont_DGls_1in12_nrmd), Fron_MAPTvsCont_mainMod, method = "be") #2;
# test_sva <- svaseq(Fron_MAPTvsCont_DGls_1in12_nrmd$counts, Fron_MAPTvsCont_mainMod, Fron_MAPTvsCont_nulMod) #4;
# y <- svaseq(Fron_MAPTvsCont_DGls_1in12_nrmd$counts, Fron_MAPTvsCont_mainMod, Fron_MAPTvsCont_nulMod) #4;

# 7.4 sva design matrix;
Fron_MAPTvsCont_mainMod_sva <- cbind(Fron_MAPTvsCont_mainMod, Fron_MAPTvsCont_sva$sv)
colnames(Fron_MAPTvsCont_mainMod_sva)[-(seq(ncol(Fron_MAPTvsCont_mainMod)))] <- paste("SV", seq(ncol(Fron_MAPTvsCont_sva$sv)), sep = "")

# # 7.5 sva plots and correlation among the covariates;
# ft <- updt_GoogleSht_FronOnly_mtchd[match(Fron_MAPTvsCont_mds_df$SampNms, updt_GoogleSht_FronOnly_mtchd$GIVENSAMPLENAME),]
# ft$RIN <- as.numeric(ft$RIN)
# ft$PH <- as.numeric(ft$PH)
# par(mfrow = c(2, 3))
# for(i in seq(Fron_MAPTvsCont_sva$n.sv)){
#   #stripchart(Fron_MAPTvsCont_sva$sv[,i] ~ ft$PH, vertical = TRUE, main = paste0("SV", i))
#   stripchart(Fron_MAPTvsCont_sva$sv[,i] ~ Fron_MAPTvsCont_DGls_1in12_nrmd$samples$group, vertical = TRUE, main = paste0("SV", i))
#   abline(h = 0)
# }
# cor(cbind(Fron_MAPTvsCont_mainMod_sva[,-1], ft[,c("Number_of_input_reads", "Ratio_of_uniquely_mapped_reads", "RIN", "AGE", "PMD.MIN.")]))
# library(psych)
# corr.test(cbind(Fron_MAPTvsCont_mainMod_sva[,-1], ft[,c("Number_of_input_reads", "Ratio_of_uniquely_mapped_reads", "RIN", "AGE", "PMD.MIN.")]))
# manovTest <- manova(Fron_MAPTvsCont_mainMod_sva[,-c(1:2)] ~ Fron_MAPTvsCont_mainMod_sva[,"Cases"])
# fron_sv_glm <- vector(mode = "list", length = 5)
# names(fron_sv_glm) <- colnames(Fron_MAPTvsCont_mainMod_sva[,-c(1:2)])
# for(i in seq(length(fron_sv_glm))){
#   fron_sv_glm[[i]] <- summary(glm(Fron_MAPTvsCont_sva$sv[,i] ~ ft$LANE + ft$AGE + as.numeric(ft$LINKERS) + ft$PMD.MIN. + ft$GENDER + ft$Number_of_input_reads + Fron_MAPTvsCont_mainMod_sva[,"MAPT"]))
# }


# 8. estimate dispersion;
Fron_MAPTvsCont_SVA_disp <- estimateDisp(Fron_MAPTvsCont_DGls_1in12_nrmd, Fron_MAPTvsCont_mainMod_sva, robust = TRUE)
dim(Fron_MAPTvsCont_SVA_disp)
# [1] 23422    28;
#common disp and bio var values;
Fron_MAPTvsCont_SVA_disp$common.dispersion
# [1] 0.1000018;
sqrt(Fron_MAPTvsCont_SVA_disp$common.dispersion)
# [1] 0.3162306;

# @plot;
pdf("Fron_MAPTvsCont_SVA_disp.pdf")
plotBCV(Fron_MAPTvsCont_SVA_disp)
dev.off()

# 9. fit glm QLF to the dispersion object;
Fron_MAPTvsCont_SVA_glmQF <- glmQLFit(Fron_MAPTvsCont_SVA_disp, design = Fron_MAPTvsCont_mainMod_sva, robust=TRUE)
dim(Fron_MAPTvsCont_SVA_glmQF)
# [1] 23422     7;
dim(Fron_MAPTvsCont_SVA_glmQF$counts)
# [1] 23422    28;

# @plot;
pdf("Fron_MAPTvsCont_SVA_glmQF.pdf")
plotQLDisp(Fron_MAPTvsCont_SVA_glmQF)
dev.off()

# 10. QLFtest;
Fron_MAPTvsCont_SVA_glmQlFtest <- glmQLFTest(Fron_MAPTvsCont_SVA_glmQF, coef = "MAPT")

# 11. Differentially expressed (DE) clusters;
Fron_MAPTvsCont_SVA_glmQlFtest_toptags <- topTags(Fron_MAPTvsCont_SVA_glmQlFtest, n=Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
dim(Fron_MAPTvsCont_SVA_glmQlFtest_toptags)
# [1] 5371   28;

# 12. Now, extract DE- protein coding mRNAs;
Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod <- Fron_MAPTvsCont_SVA_glmQlFtest_toptags$table[which(Fron_MAPTvsCont_SVA_glmQlFtest_toptags$table$'CAT_geneClass' == "coding_mRNA"),]
dim(Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod)
# [1] 4837   28;

# 13. and, get the promoters only;
Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod_promt <- Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod[grep("Promoter", Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod$annotation),]
dim(Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod_promt)
# [1] 4780   28;

# 13.1 Up-regulated;
Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod_promtUp <- Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod_promt[which(Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod_promt$logFC > 0), ]
dim(Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod_promtUp)
# [1] 1652   28;
range(Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod_promtUp$logFC)
# [1] 0.3464392 2.9932309;

# 13.2 Down-regulated;
Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod_promtDn <- Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod_promt[which(Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod_promt$logFC < 0), ]
dim(Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod_promtDn)
# [1] 3128   28;
range(Fron_MAPTvsCont_SVA_glmQlFtest_toptags_protCod_promtDn$logFC)
# [1] -3.791878 -0.351334;



# #######################;
# IV. Save this workspace;
# #######################;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc_fro")
# save.image("S7_FronSamps_MAPTvsCont_RikAnn_SVA_DEG.RData")
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc_fro/S7_FronSamps_MAPTvsCont_RikAnn_SVA_DEG.RData")
# END ______________________________________________________________________________________________________________________________________________________________________________________;
