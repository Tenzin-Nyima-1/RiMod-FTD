# Date: 23 Jan, 2017;
# Project: RIMOD-FTD;
# Aim: file prep/crosscheck prior to differential gene expression analysis;
# NOTE: the annotation of the count-table was carried using RIKEN's CAT gtf file as the gene model;
# NOTE 2: the same script could be adapted for the count table from: all 7 brain regions; and, temporal brain; 
# Workflow (main steps);
# Step 1: Upload the workspace that performed annotaion of the count table;
# Step 2: crosscheck if the reordered annotated count table matches the original count table;
# Step 3: when cleared, to save space, keep only the required objects;
# Step 4: save this workspace: used for DGE analysis;


options(stringsAsFactors=FALSE)
library(stringi)
library(stringr)
library(plotrix)
library(readxl)
library(Hmisc)
library(reshape2)


# ################################################################################################;
# I. Upload the workspace that annotated fontal count table using F6CAT gtf file as the gene-model;
# ################################################################################################;
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc/FronSamps_Annotbedfile_F6CATfile.RData")



# ##################################;
# II. all samples including sporadic;
# ##################################;


# [A] Google sheet with library info;

# 1. For all the samples;
dim(updt_GoogleSht_6Dec2016)
# [1] 247  32;

# 2. For the frontal samples only;
dim(updt_GoogleSht_FronOnly_mtchd)
# [1] 58 22;

# 3. Crosschecks: while subsetting lib info for frontal samples, check if the right files were chosen and if any of the infos were left out;
all(updt_GoogleSht_6Dec2016[grep("frontal", updt_GoogleSht_6Dec2016$REGION),colnames(updt_GoogleSht_FronOnly_mtchd)][rownames(updt_GoogleSht_FronOnly_mtchd),]
    == updt_GoogleSht_FronOnly_mtchd)
# [1] TRUE;
# NOTE: a subset of the frontal library was correctly generated;



# [B] Count tables;

# 1. Main frontal sample count table;  
dim(FronOnly_counTab)
#[1] 31000    58;

# 2. Frontal table with column names appended with study category names;
dim(FronOnly_counTab_new)
#[1] 31000    58;

# 3.1 Crosscheck: if the counts in the above 2 tables are similar;
all(FronOnly_counTab == FronOnly_counTab_new)
# [1] TRUE;
identical(stri_split_fixed(colnames(FronOnly_counTab), "_no_chrM", simplify = T)[,1], stri_split_fixed(colnames(FronOnly_counTab_new), "_no_chrM_", simplify = T)[,1])
# [1] TRUE;
# NOTE: this confirmed that the above two count tables and the order of their colnames are correct;

# 3.2 crosscheck if the study category name appended at the colnames are correct;
identical(stri_split_fixed(colnames(FronOnly_counTab_new), "_no_chrM_", simplify = T)[,2], 
          updt_GoogleSht_FronOnly_mtchd$MUTATED.GENE[match(stri_split_fixed(stri_split_fixed(colnames(FronOnly_counTab_new), "_no_chrM_", simplify = T)[,1], "sample_", simplify = T)[,2], updt_GoogleSht_FronOnly_mtchd$GIVENSAMPLENAME)])
# [1] TRUE;
# NOTE: this confirmed that the study category names appended to the colnames were correctly inferenced from the library info file;



# [C] NOTE: since the crosschecks were cleared, from here on, use the following two files for inference;

# 1. Google sheet with info for all the libraries from frontal including the sporadic samples;
dim(updt_GoogleSht_FronOnly_mtchd)
# [1] 58 22;

# 2. main frontal count table with study categories specified in the colnames, and including the sporadic samples;
dim(FronOnly_counTab_new)
#[1] 31000    58;



# ###############################;
# III. excluding sporadic samples;
# ###############################;


# [A] checking if the right sporadic sample were excluded;

# 1. Get the given sample names of the sporadic samples;
identical(updt_GoogleSht_FronOnly_mtchd[grep("spor", updt_GoogleSht_FronOnly_mtchd$MUTATED.GENE, ignore.case = T),], 
          updt_GoogleSht_FronOnly_mtchd[grep("spor", updt_GoogleSht_FronOnly_mtchd$DISEASE.CODE, ignore.case = T),])
# [1] TRUE;
fron_spoGivSampNms <- updt_GoogleSht_FronOnly_mtchd$GIVENSAMPLENAME[grep("spor", updt_GoogleSht_FronOnly_mtchd$MUTATED.GENE, ignore.case = T)]
length(fron_spoGivSampNms)
# [1] 9;

# 2. The frontal count table but excluding the sporadic samples
# @checking if the colnames of the main frontal count table after removing sporadic ones are simlar to the colnames of frontal counts without the sporadic samples;
identical(colnames(FronOnly_counTab_new[, -(grep(paste(fron_spoGivSampNms, collapse = "|"), colnames(FronOnly_counTab_new), ignore.case = TRUE))]), 
          colnames(FronOnly_counTab_minSpo))
# [1] TRUE;
# @since the colnames are in the right order, checking if the counts are same in the main and the subsetted ones;
identical(FronOnly_counTab_new[, -(grep(paste(fron_spoGivSampNms, collapse = "|"), colnames(FronOnly_counTab_new), ignore.case = T))], 
          FronOnly_counTab_minSpo)
# [1] TRUE;
# NOTE: yes, the right sporadic samples were removed and subsetted;

# 3. the above frontal count table wthout the sporadic samples;
dim(FronOnly_counTab_minSpo)
#[1] 31000    49;

# @a tiny crosscheck if the sporadic samples are present;
grep(paste(fron_spoGivSampNms, collapse = "|"), colnames(FronOnly_counTab_minSpo), ignore.case = TRUE)
# integer(0);

# @check the disease code of the samples;
unique(updt_GoogleSht_FronOnly_mtchd$DISEASE.CODE[match(stri_split_fixed(stri_split_fixed(colnames(FronOnly_counTab_minSpo), "sample_", simplify = T)[,2], "_no_chrM_", simplify = T)[,1], updt_GoogleSht_FronOnly_mtchd$GIVENSAMPLENAME)])
# [1] "FTD-C9"   "FTD-MAPT" "control"  "FTD-GRN";
unique(updt_GoogleSht_FronOnly_mtchd$MUTATED.GENE[match(stri_split_fixed(stri_split_fixed(colnames(FronOnly_counTab_minSpo), "sample_", simplify = T)[,2], "_no_chrM_", simplify = T)[,1], updt_GoogleSht_FronOnly_mtchd$GIVENSAMPLENAME)])
# [1] "C9orf72" "MAPT"    "control" "GRN";
# NOTE: sporadic samples were correctly excluded from the current frontal count table;



# [B] The annotated count table;

# 1. annotated count table;
dim(FronOnly_peak_GR_3k_DF_ord)
# [1] 31000    66;

# 2. crosscheck, if the counts in the annotated table match the sporadic sample removed main count table;
identical(colnames(FronOnly_counTab_minSpo), colnames(FronOnly_peak_GR_3k_DF_ord)[grep("sample", colnames(FronOnly_peak_GR_3k_DF_ord))])
# [1] TRUE;
identical(FronOnly_counTab_minSpo, FronOnly_peak_GR_3k_DF_ord[,grep("sample", colnames(FronOnly_peak_GR_3k_DF_ord))])
# [1] TRUE;
identical(FronOnly_counTab_minSpo, FronOnly_peak_GR_3k_DF_ord[,colnames(FronOnly_counTab_minSpo)])
# [1] TRUE;
# NOTE: yes, the counts are correct;

# 2.1 checking if the rownames are similar to clusid cols;
identical(rownames(FronOnly_peak_GR_3k_DF_ord), as.character(FronOnly_peak_GR_3k_DF_ord$Clus_ID))
# [1] TRUE;
identical(as.integer(stri_split_fixed(as.character(FronOnly_peak_GR_3k_DF_ord$Clus_ID), "_", simplify = T)[,3]), FronOnly_peak_GR_3k_DF_ord$end)
# [1] TRUE;
identical(as.character(stri_split_fixed(as.character(FronOnly_peak_GR_3k_DF_ord$Clus_ID), "_", simplify = T)[,4]), as.character(FronOnly_peak_GR_3k_DF_ord$strand))
# [1] TRUE;
identical(as.character(stri_split_fixed(as.character(FronOnly_peak_GR_3k_DF_ord$Clus_ID), "_", simplify = T)[,4]), as.character(FronOnly_peak_GR_3k_DF_ord$Strand))
# [1] TRUE;
identical(as.character(stri_split_fixed(as.character(FronOnly_peak_GR_3k_DF_ord$Clus_ID), "_", simplify = T)[,1]), as.character(FronOnly_peak_GR_3k_DF_ord$seqnames))
# [1] TRUE;
# NOTE: clusid row names and clus ids details in the columns are correct;

# 2.2 crosscheck the disease code of the sporadic excluded count table;
unique(updt_GoogleSht_FronOnly_mtchd$DISEASE.CODE[match(stri_split_fixed(stri_split_fixed(colnames(FronOnly_peak_GR_3k_DF_ord)[grep("sample", colnames(FronOnly_peak_GR_3k_DF_ord))], "sample_", simplify = T)[,2], "_no_chrM_", simplify = T)[,1], updt_GoogleSht_FronOnly_mtchd$GIVENSAMPLENAME)])
# [1] "FTD-C9"   "FTD-MAPT" "control"  "FTD-GRN";
unique(updt_GoogleSht_FronOnly_mtchd$MUTATED.GENE[match(stri_split_fixed(stri_split_fixed(colnames(FronOnly_peak_GR_3k_DF_ord)[grep("sample", colnames(FronOnly_peak_GR_3k_DF_ord))], "sample_", simplify = T)[,2], "_no_chrM_", simplify = T)[,1], updt_GoogleSht_FronOnly_mtchd$GIVENSAMPLENAME)])
# [1] "C9orf72" "MAPT"    "control" "GRN";



# [C] NOTE: since the crosschecks were cleared, from here on, use the annotated tabler for DEG;
dim(FronOnly_peak_GR_3k_DF_ord)
# [1] 31000    66;



# #################################################;
# V. order the columns of the annotated count table;
# #################################################;
dim(FronOnly_peak_GR_3k_DF_ord)
# [1] 31000    66;

# 1. the unique study category names;
unique(stri_split_fixed(colnames(FronOnly_peak_GR_3k_DF_ord)[grep("sample", colnames(FronOnly_peak_GR_3k_DF_ord), ignore.case = T)], "_no_chrM_", simplify = T)[,2])
# [1] "C9orf72" "MAPT"    "control" "GRN";

# 2. now order the table;
Ann_FronOnly_peakGR3k_ordFinal <- FronOnly_peak_GR_3k_DF_ord[, as.character(c(colnames(FronOnly_peak_GR_3k_DF_ord)[-(grep("sample", colnames(FronOnly_peak_GR_3k_DF_ord), ignore.case = T))], 
                                                                              colnames(FronOnly_peak_GR_3k_DF_ord)[grep("cont", colnames(FronOnly_peak_GR_3k_DF_ord), ignore.case = T)], 
                                                                              colnames(FronOnly_peak_GR_3k_DF_ord)[grep("MAPT", colnames(FronOnly_peak_GR_3k_DF_ord), ignore.case = T)], 
                                                                              colnames(FronOnly_peak_GR_3k_DF_ord)[grep("C9orf72", colnames(FronOnly_peak_GR_3k_DF_ord), ignore.case = T)], 
                                                                              colnames(FronOnly_peak_GR_3k_DF_ord)[grep("GRN", colnames(FronOnly_peak_GR_3k_DF_ord), ignore.case = T)]))]
dim(Ann_FronOnly_peakGR3k_ordFinal)
# [1] 31000    66;

# 3. Crosscheck if the orderd is similar to the unordered annotated table;
identical(FronOnly_peak_GR_3k_DF_ord, Ann_FronOnly_peakGR3k_ordFinal[,colnames(FronOnly_peak_GR_3k_DF_ord)])
# [1] TRUE;
# NOTE: yes, they are similar;



# ######################################;
# VI. final crosscheck cleared files are;
# ######################################;

# 1. Google sheet with info for all the libraries from frontal including the sporadic samples;
dim(updt_GoogleSht_FronOnly_mtchd)
# [1] 58 22;

# sample count;
sort(table(updt_GoogleSht_FronOnly_mtchd$MUTATED.GENE[-(grep("Spo", updt_GoogleSht_FronOnly_mtchd$MUTATED.GENE))]))
# GRN    MAPT C9orf72 control 
# 8      12      13      16

# # 2. frontal count table wthout the sporadic samples;
# dim(FronOnly_counTab_minSpo)
# #[1] 31000    49;

# # 3. the annotated table;
# dim(FronOnly_peak_GR_3k_DF_ord)
# # [1] 31000    66;

# 4. the ordered annotated table;
dim(Ann_FronOnly_peakGR3k_ordFinal)
# [1] 31000    66;



# ###################################;
# VII. keep only the required objects;
# ###################################;
length(ls())
# [1] 33;
rm(list = ls()[-(match(c("Ann_FronOnly_peakGR3k_ordFinal", "updt_GoogleSht_FronOnly_mtchd"), ls()))])
ls()
# [1] "Ann_FronOnly_peakGR3k_ordFinal" "updt_GoogleSht_FronOnly_mtchd";



# #########################;
# VIII. save this workspace;
# #########################;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc_fro")
save.image("S6_FronSamps_preDGEanaly_crschk_orgnz.RData")
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc_fro/S6_FronSamps_preDGEanaly_crschk_orgnz.RData")
# END_______________________________________________________________________________________________________________________________________________________________________________________;

