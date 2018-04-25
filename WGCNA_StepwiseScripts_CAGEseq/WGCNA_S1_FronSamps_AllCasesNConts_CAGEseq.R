# Date: 20-21 March, 2018;
# Project: RIMOD-FTD;
# Aim: Weighted Gene Coexpression Network Analysis(WGCNA) of frontal samples: all cases and control samples;
# NOTE: the annotation of the count-table was carried using RIKEN's CAT gtf file as the gene model;
# NOTE 2: the same script could be adapted for the count table from: all 7 brain regions; and, temporal brain; 
# Workflow (main steps);
# Step 1: Data exploration: outlier detection and data cleaning;



# source("http://bioconductor.org/biocLite.R") 
# biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
# install.packages("WGCNA")
options(stringsAsFactors=FALSE)
library(AnnotationDbi)
library(edgeR)
library(stringi)
library(stringr)
library(WGCNA)
library(edgeR)
library(impute)
library(GO.db)
library(preprocessCore)



# ################################################################################################;
# I. Upload the workspace that performed pre-DGE analysis cleaning & organizing of frontal samples;
# ________________________________________________________________________________________________;
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc_fro/S6_FronSamps_preDGEanaly_crschk_orgnz.RData")

# 1. recall the required objects; 

# @ordered annotated table;
dim(Ann_FronOnly_peakGR3k_ordFinal)
# [1] 31000    66;

# @library info for the frontal samples only;
dim(updt_GoogleSht_FronOnly_mtchd)
# [1] 58 22;


# 2. Couple checks if the sporadics are removed properly;

# crosscheck1: sporadic;
Fro_Spo_Samps <- updt_GoogleSht_FronOnly_mtchd$GIVENSAMPLENAME[which(updt_GoogleSht_FronOnly_mtchd$MUTATED.GENE == "Sporadic")]
length(Fro_Spo_Samps)
# [1] 9;
unlist(lapply(seq(Fro_Spo_Samps), function(x) grep(Fro_Spo_Samps[x], colnames(Ann_FronOnly_peakGR3k_ordFinal))))
# integer(0);

# crosscheck2: non-sporadic;
Fro_NonSpo_Samps <- updt_GoogleSht_FronOnly_mtchd$GIVENSAMPLENAME[-(which(updt_GoogleSht_FronOnly_mtchd$MUTATED.GENE == "Sporadic"))]
length(Fro_NonSpo_Samps)
# [1] 49;
unique(unlist(lapply(lapply(seq(Fro_NonSpo_Samps), function(x) grep(Fro_NonSpo_Samps[x], colnames(Ann_FronOnly_peakGR3k_ordFinal))), length)))
# [1] 1;

# crosscheck3: intersect between the 2 lists - sporadic and non-sporadic;
intersect(Fro_Spo_Samps, Fro_NonSpo_Samps)
# character(0);
unique(unlist(lapply(lapply(seq(Fro_Spo_Samps), function(x) grep(Fro_Spo_Samps[x], colnames(Ann_FronOnly_peakGR3k_ordFinal))), length)))
# [1] 0;
# Note: meaning that none of the sporadic samples are present in the count table;

# @Final;
# NOTE: meaning that all sporadic are excluded, and all non sporadic ones are included;

# 4. clean space;
# @keeping only the required objects but remove the rest;
length(ls())
# [1] 4;
rm(list = ls()[-(match(c("Ann_FronOnly_peakGR3k_ordFinal", "updt_GoogleSht_FronOnly_mtchd"), ls()))])
ls()
# [1] "Ann_FronOnly_peakGR3k_ordFinal"; "updt_GoogleSht_FronOnly_mtchd";
######################################################################################################



# ########################################;
# II. Data preparation using edgeR library;
# ________________________________________;


# [A] Extract the count-table;
Fron_AlCasesVsCont_counts <- Ann_FronOnly_peakGR3k_ordFinal[, grep("sample", colnames(Ann_FronOnly_peakGR3k_ordFinal), ignore.case = T)]
dim(Fron_AlCasesVsCont_counts)
# [1] 31000    49;


# [B] Define groups for each study category; 
Fron_AlCasesVsCont_Gps <- stri_split_fixed(colnames(Fron_AlCasesVsCont_counts), "_no_chrM_", simplify = TRUE)[,2]
Fron_AlCasesVsCont_GpsOrd <- factor(ifelse(Fron_AlCasesVsCont_Gps == "control", "cont", "case"), levels = c("cont", "case"))

# @infos;
length(Fron_AlCasesVsCont_GpsOrd)
# [1] 49;
identical(length(Fron_AlCasesVsCont_Gps), length(Fron_AlCasesVsCont_GpsOrd))
# [1] TRUE;
table(Fron_AlCasesVsCont_GpsOrd)
# Fron_AlCasesVsCont_GpsOrd
# cont case 
# 16   33  


# [C] Gene info file;
Fron_AlCasesVsCont_genes <- Ann_FronOnly_peakGR3k_ordFinal[,-(grep("sample", colnames(Ann_FronOnly_peakGR3k_ordFinal)))]
dim(Fron_AlCasesVsCont_genes)
# [1] 31000    17;

# 1. check if the duplicated gene columns are identical;
identical(rownames(Fron_AlCasesVsCont_genes), as.character(Fron_AlCasesVsCont_genes$Clus_ID))
# [1] TRUE;
identical(as.character(Fron_AlCasesVsCont_genes$seqnames), stri_split_fixed(Fron_AlCasesVsCont_genes$Clus_ID, "_", simplify = T)[,1])
# [1] TRUE;
identical(as.character(Fron_AlCasesVsCont_genes$geneChr), stri_split_fixed(Fron_AlCasesVsCont_genes$Clus_ID, "_", simplify = T)[,1])
# [1] TRUE;
identical(as.character(Fron_AlCasesVsCont_genes$strand), stri_split_fixed(Fron_AlCasesVsCont_genes$Clus_ID, "_", simplify = T)[,4])
# [1] TRUE;
identical(as.character(Fron_AlCasesVsCont_genes$Strand), stri_split_fixed(Fron_AlCasesVsCont_genes$Clus_ID, "_", simplify = T)[,4])
# [1] TRUE;
identical(as.character(Fron_AlCasesVsCont_genes$geneStrand), stri_split_fixed(Fron_AlCasesVsCont_genes$Clus_ID, "_", simplify = T)[,4])
# [1] TRUE;
# NOTE: So, keep only one strand info; and, discard-seqnames,score;  

# 2. the ordered genes file;
Fron_AlCasesVsCont_GenesOrd <- Fron_AlCasesVsCont_genes[,c("Clus_ID", "seqnames", "start", "end", "strand", "width", "annotation", "distanceToTSS", "geneId", "transcriptId", "geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand")]
dim(Fron_AlCasesVsCont_GenesOrd)
# [1] 31000    15;

# 3. read the F6 geneinfo(from chang on 12 Jan 2017) file, in order to get the alternative converted ids for each gene-id;
setwd("/media/student-ag-javier/Data/Tenzin/Fantom_6/Datas/F6_CAT_GTFfile/F6_CAT_GeneInf_RecvOn12Jan17")
F6_CATgene_info <- read.delim("F6_CAT.gene.info.tsv",  header = TRUE, sep = "\t", quote = "")
dim(F6_CATgene_info)
#[1] 124047     41;

# 4. Get additional gene info from the above CAT file;

# i) cbind geneInfo using F6 file;
which(duplicated(F6_CATgene_info$geneID))
# integer(0) # NOTE: there isn't any duplicated gene-ids so its save to match with them;
Fron_AlCasesVsCont_GenesOrd_F6 <- cbind(Fron_AlCasesVsCont_GenesOrd, 
                                        F6_CATgene_info[match(Fron_AlCasesVsCont_GenesOrd$geneId, F6_CATgene_info$geneID), 
                                                        c("geneID", "geneName", "geneType", "CAT_geneClass", "CAT_DHS_type", "trnscptID", "HGNC_name", "HGNC_symbol", "entrez_ID", "refseq_ID")])
dim(Fron_AlCasesVsCont_GenesOrd_F6)
# [1] 31000    25;

# ii) crosscheck;

# @right gene-ids matched?;
identical(Fron_AlCasesVsCont_GenesOrd_F6$geneId, Fron_AlCasesVsCont_GenesOrd_F6$geneID)
# [1] TRUE;

# @right TX-ids?;
Fron_AlCasesVsCont_GenesOrd_F6_TXchk <- unlist(lapply(lapply(seq(nrow(Fron_AlCasesVsCont_GenesOrd_F6)), function(x) grep(Fron_AlCasesVsCont_GenesOrd_F6$transcriptId[x], Fron_AlCasesVsCont_GenesOrd_F6$trnscptID[x])), length))
length(Fron_AlCasesVsCont_GenesOrd_F6_TXchk)
# [1] 31000;
unique(Fron_AlCasesVsCont_GenesOrd_F6_TXchk)
# [1] 1;

# NOTE: both gene-ids and the TX ids matched and confirmed that they are in the same order, so remove the extra columns from the data-frame;

# iii) Cleared crosscheck, so remove the extra columns;
Fron_AlCasesVsCont_F6Genes_fin <- Fron_AlCasesVsCont_GenesOrd_F6[,-(match(c("geneID", "trnscptID"), colnames(Fron_AlCasesVsCont_GenesOrd_F6)))]
dim(Fron_AlCasesVsCont_F6Genes_fin)
# [1] 31000    23;
range(Fron_AlCasesVsCont_F6Genes_fin$width)
# [1]    2 2915;



# ######################################################################################################;
# III. Using edgeR remove the lowly expressed genes and prepare the data suitable for the WGCNA analysis;
# ______________________________________________________________________________________________________;

# 1. create a dgelist object;
Fron_AlCasesVsCont_DGls <- DGEList(counts = Fron_AlCasesVsCont_counts, group = Fron_AlCasesVsCont_GpsOrd, genes = Fron_AlCasesVsCont_F6Genes_fin)
dim(Fron_AlCasesVsCont_DGls)
# [1] 31000    49;

# @a crosscheck;
identical(rownames(Fron_AlCasesVsCont_DGls), as.character(Fron_AlCasesVsCont_DGls$genes$Clus_ID))
# [1] TRUE;

# 2. Select geneTYPE:  "protein_coding" ones only ("from GENCODE, plus the novel lncRNA knockdown targets with tagged as novel gene");
table(Fron_AlCasesVsCont_DGls$genes$geneType)
# antisense; IG_V_gene; lincRNA; miRNA; misc_RNA; __na; polymorphic_pseudogene; processed_transcript;
# 1012       8          899      58     17        1741  8                       225 ;
# protein_coding; pseudogene; rRNA; sense_intronic; sense_overlapping; snoRNA; snRNA; TR_C_gene ;
# 26158         ; 752         9     17              43                 30      19     1 ;
# TR_J_gene; TR_V_gene; 
# 2        ;  1; 

# @crosscheck;
identical(which(Fron_AlCasesVsCont_DGls$genes$geneType == "protein_coding"), grep("protein_", Fron_AlCasesVsCont_DGls$genes$geneType, ignore.case = T))
# [1] TRUE;

# @selection;
Fron_AlCasesVsCont_DGls <- Fron_AlCasesVsCont_DGls[which(Fron_AlCasesVsCont_DGls$genes$geneType == "protein_coding"), ]
dim(Fron_AlCasesVsCont_DGls)
# [1] 26158    49;

# @another crosscheck;
table(Fron_AlCasesVsCont_DGls$genes$geneType)
# protein_coding;
# 26158;

# 3. Select annotation = "promoters" only;
table(substring(Fron_AlCasesVsCont_DGls$genes$annotation, 1, 5))
# Dista Downs Exon  Intro Promo;
# 44     2   405   110 25597;

# @crosscheck;
identical(grep("Promoter", Fron_AlCasesVsCont_DGls$genes$annotation, ignore.case = T), grep("Promoter", Fron_AlCasesVsCont_DGls$genes$annotation))
# [1] TRUE;

# @selection;
Fron_AlCasesVsCont_DGls <- Fron_AlCasesVsCont_DGls[grep("Promoter", Fron_AlCasesVsCont_DGls$genes$annotation, ignore.case = T),]
dim(Fron_AlCasesVsCont_DGls)
# [1] 25597    49;

# @another crosscheck;
identical(rownames(Fron_AlCasesVsCont_DGls), as.character(Fron_AlCasesVsCont_DGls$genes$Clus_ID))
# [1] TRUE;
table(Fron_AlCasesVsCont_DGls$genes$annotation)
# Promoter (1-2kb) Promoter (<=1kb) Promoter (2-3kb);
# 515            24813              269;


# 4. filter out lowly expressed TXs;

# i) setting 1cpm threshold to estimate raw count range for that value;
range((1 * Fron_AlCasesVsCont_DGls$samples$lib.size)/1e6)
# [1]  3.081283 32.214399;

# ii) get the #sample in both cont and case categories;
table(Fron_AlCasesVsCont_DGls$samples$group)
# cont case;
# 16   33 ;

# iii) since the lowest #sample is 16, set the filter for TX with <1cpm in that many samples;
Fron_AlCasesVsCont_1in16_keepIDs <- rowSums(cpm(Fron_AlCasesVsCont_DGls) > 1) >= 16
summary(Fron_AlCasesVsCont_1in16_keepIDs)
# Mode   FALSE    TRUE    NA's ;
# logical    4546   21051       0 ;

# iv) now select only those TXs that passed the threshold;
Fron_AlCasesVsCont_DGls_1in16 <- Fron_AlCasesVsCont_DGls[Fron_AlCasesVsCont_1in16_keepIDs, , keep.lib.sizes=FALSE]
dim(Fron_AlCasesVsCont_DGls_1in16)
# [1] 21051    49;

# @crosscheck if the library sizes were kept similar;
all(Fron_AlCasesVsCont_DGls_1in16$samples$lib.size == Fron_AlCasesVsCont_DGls$samples$lib.size)
# [1] FALSE; NOTE: correct, meaning that lib-sizes were adjusted after the filtering;



# ################################################################;
# IV. WGCNA: Data exploration: outlier detection and data cleaning;
# ________________________________________________________________;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Fron_WGCNA/Unsigned/All_CasesVSCont/Step1")

# 1. normalize using RLE with filtered transcripts only;
Fron_AlCasesVsCont_DGls_1in16_nrmd <- calcNormFactors(Fron_AlCasesVsCont_DGls_1in16, method = "RLE")
dim(Fron_AlCasesVsCont_DGls_1in16_nrmd)
# [1] 21051    49;

# @crosschecks;
identical(rownames(Fron_AlCasesVsCont_DGls_1in16_nrmd), as.character(Fron_AlCasesVsCont_DGls_1in16_nrmd$genes$Clus_ID))
# [1] TRUE;
identical(Fron_AlCasesVsCont_DGls_1in16_nrmd$samples$lib.size, Fron_AlCasesVsCont_DGls_1in16$samples$lib.size)
# [1] TRUE;
identical(as.numeric(colSums(Fron_AlCasesVsCont_DGls_1in16_nrmd$counts)), Fron_AlCasesVsCont_DGls_1in16_nrmd$samples$lib.size)
# [1] TRUE;

# 2. since wgcna was designed for a continuous data, transform the data that is suitable for such analysis;
Fro_alCasesVsCont_normd_logcpm <- cpm(Fron_AlCasesVsCont_DGls_1in16_nrmd, prior.count=2, log=TRUE)
dim(Fro_alCasesVsCont_normd_logcpm)
# [1] 21051    49;

# 3. prepare the data for wgcna;
Fro_alCasesVsCont_normd_logcpm_Trs <- as.data.frame(t(Fro_alCasesVsCont_normd_logcpm))
dim(Fro_alCasesVsCont_normd_logcpm_Trs)
# [1]    49 21051;

# @crosscheck;
identical(colnames(Fro_alCasesVsCont_normd_logcpm_Trs), as.character(Fron_AlCasesVsCont_DGls_1in16_nrmd$genes$Clus_ID))
# [1] TRUE

# 4. check the data for missing values;
Fro_alCasesVsCont_normd_logcpm_Trs_gsg <- goodSamplesGenes(Fro_alCasesVsCont_normd_logcpm_Trs, verbose = 3)
Fro_alCasesVsCont_normd_logcpm_Trs_gsg$allOK
# [1] TRUE;


# 5. all samples/genes are ok, so check for outliers;
Fro_alCasesVsCont_SmpTree <- hclust(dist(Fro_alCasesVsCont_normd_logcpm_Trs), method = "average")
# NOTE: the distance was measured using euclidean method and the clustering was performed using the method "average";

# 5.1 explore the plot and remove if any sample is an outlier;
pdf("Fro_alCasesVsCont_SmpTree_normdlogcpm", width = 12, height = 9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(Fro_alCasesVsCont_SmpTree, main = "Outlier diagnosis", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# @specify the cut(line)  threshold;
abline(h = 250, col = "red");
dev.off()

# 5.2 Now, select the clusters below the line;
Fro_alCasesVsCont_clust <- cutreeStatic(Fro_alCasesVsCont_SmpTree, cutHeight = 250, minSize = 3)
table(Fro_alCasesVsCont_clust)
# Fro_alCasesVsCont_clust
# 0  1 
# 1 48 
Fro_alCasesVsCont_SmpTree$labels[which(Fro_alCasesVsCont_clust == 0)]
# [1] "sample_05108_fro_no_chrM_MAPT";

# 5.3 cluster 1 contains all the required samples;
reqFro_alCasesVsCont_normd_logcpm_Trs <- Fro_alCasesVsCont_normd_logcpm_Trs[Fro_alCasesVsCont_clust==1, ]
dim(reqFro_alCasesVsCont_normd_logcpm_Trs)
# [1]    48 21051;
(Fro_alCasesVsCont_nGenes <- ncol(reqFro_alCasesVsCont_normd_logcpm_Trs))
# [1] 21051;
(Fro_alCasesVsCont_nSamples <- nrow(reqFro_alCasesVsCont_normd_logcpm_Trs))
# [1] 48;


# 6. Map clinical traits to the data;

# 6.1 extract the sample ids;
reqFro_AlCaseVSCont_SampIDs <- stri_split_fixed(stri_split_fixed(rownames(reqFro_alCasesVsCont_normd_logcpm_Trs), "sample_", simplify = T)[,2], "_no", simplify = T)[,1]
length(reqFro_AlCaseVSCont_SampIDs)
# [1] 48;

# 6.2 recall the clinical data and arrange it as per the order of the samples;
dim(updt_GoogleSht_FronOnly_mtchd)
# [1] 58 22;

# @trait file with the required samples only;
reqFroAlCaseVsCont_TraitInfo <- updt_GoogleSht_FronOnly_mtchd[match(reqFro_AlCaseVSCont_SampIDs, updt_GoogleSht_FronOnly_mtchd$GIVENSAMPLENAME), 
                                                              c("CASE.CONTROL", "MUTATED.GENE", "GIVENSAMPLENAME", "FLOWCELL", "LANE", "Number_of_input_reads", "Number_of_Uniquely_mapped_reads", "Ratio_of_uniquely_mapped_reads", "RIN", "LINKERS", "AGE", "GENDER", "PMD.MIN.", "PH")]
rownames(reqFroAlCaseVsCont_TraitInfo) <- rownames(reqFro_alCasesVsCont_normd_logcpm_Trs)
dim(reqFroAlCaseVsCont_TraitInfo)
# [1] 48 14; 

# --@crosscheck;
# the same order?;
identical(stri_split_fixed(rownames(reqFroAlCaseVsCont_TraitInfo), "chrM_", simplify = T)[,2], stri_split_fixed(reqFroAlCaseVsCont_TraitInfo$MUTATED.GENE, "_", simplify = T)[,1])
# [1] TRUE;
identical(stri_split_fixed(stri_split_fixed(rownames(reqFroAlCaseVsCont_TraitInfo), "sample_", simplify = T)[,2], "_no_", simplify = T)[,1], reqFroAlCaseVsCont_TraitInfo$GIVENSAMPLENAME)
# [1] TRUE;
identical(stri_split_fixed(stri_split_fixed(rownames(reqFroAlCaseVsCont_TraitInfo), "sample_", simplify = T)[,2], "_no_", simplify = T)[,1], reqFro_AlCaseVSCont_SampIDs)
# [1] TRUE;

# ---@after checking the consistencies, remove the col -given sample name;
finFroAlCaseVsCont_TraitInfo <- reqFroAlCaseVsCont_TraitInfo[, -(grep("GIVENSAMPLENAME", colnames(reqFroAlCaseVsCont_TraitInfo)))]
dim(finFroAlCaseVsCont_TraitInfo)
# [1] 48 13;

# ----@convert the trait class into numeric;
finFroAlCaseVsCont_TraitInfo$CASE.CONTROL <- with(finFroAlCaseVsCont_TraitInfo, ifelse(CASE.CONTROL == "control", 0, 1))

# -----@convert the mutated category into numeric and in line with case.control status;
finFroAlCaseVsCont_TraitInfo_MutCat <- unique(finFroAlCaseVsCont_TraitInfo$MUTATED.GENE)
finFroAlCaseVsCont_TraitInfo_MutCat_num <- 0:(length(finFroAlCaseVsCont_TraitInfo_MutCat)-1)
# @replace;
for(i in seq(finFroAlCaseVsCont_TraitInfo_MutCat)) {
  finFroAlCaseVsCont_TraitInfo$MUTATED.GENE[grep(finFroAlCaseVsCont_TraitInfo_MutCat[i], finFroAlCaseVsCont_TraitInfo$MUTATED.GENE, ignore.case = T)] <- finFroAlCaseVsCont_TraitInfo_MutCat_num[i]
}
finFroAlCaseVsCont_TraitInfo$MUTATED.GENE <- as.numeric(finFroAlCaseVsCont_TraitInfo$MUTATED.GENE)

# ------@now convert the classes for the rest of the cols;
for(i in seq(ncol(finFroAlCaseVsCont_TraitInfo))){
  if(is.numeric(finFroAlCaseVsCont_TraitInfo[,i])) {
    finFroAlCaseVsCont_TraitInfo[,i] <- as.numeric(finFroAlCaseVsCont_TraitInfo[,i])
  } else {
    finFroAlCaseVsCont_TraitInfo[,i] <- as.numeric(factor(finFroAlCaseVsCont_TraitInfo[,i]))
  }
}

# -------@Crosschecks if the transformation was done correctly;
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

collectGarbage();


# 8. Recluster the required samples and map each to the corresponding traits;
reqFro_alCasesVsCont_SmpTree <- hclust(dist(reqFro_alCasesVsCont_normd_logcpm_Trs), method = "average")
# NOTE: cluster method = average; and, distance=euclidean; 

# @convert the trait features to colors: white means low, red means high, grey means missing entry;
identical(rownames(reqFro_alCasesVsCont_normd_logcpm_Trs), rownames(finFroAlCaseVsCont_TraitInfo))
# [1] TRUE;
finFroAlCaseVsCont_TraitCols <- numbers2colors(finFroAlCaseVsCont_TraitInfo, signed = FALSE)

# -@plot;
pdf("reqFro_alCasesVsCont_SmpTraitTree_2_Avg", width = 15, height = 12)
par(cex = 0.5);
plotDendroAndColors(reqFro_alCasesVsCont_SmpTree, finFroAlCaseVsCont_TraitCols, groupLabels = names(finFroAlCaseVsCont_TraitInfo), main = "Frontal sample dendrogram and trait heatmap")
dev.off()



# #####################;
# V. Save the workspace;
# _____________________;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc")
save.image("WGCNA_S1_FronSamps_AllCasesNConts_CAGEseq.RData")
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc/WGCNA_S1_FronSamps_AllCasesNConts_CAGEseq.RData")
# END ____________________________________________________________________________________________________________________________________________________________________________________________;