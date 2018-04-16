# Date: 23 Jan, 2017;
# Project: RiMod-FTD;
# Aim: to generate BED files from the count table from frontal brain;
# NOTE: for now, sporadic samples are to be excluded from the analysis; 
# -also the same script could be adapted for the count table from: all 7 brain regions; and, temporal brain; 
# --the integrity of the count-tables were checked using md5sum, check script: "md5sum_countTabs_clusNwrksp.sh";
# Workflow;
# Step 1: upload and explore the metadata; use its info to identify the cases/cont status of the samples;
# Step 2: generate a bed file of the count table: assigning the required col-classes to the key cols;


library(stringi)
library(stringr)
library(xlsx)
library(muStat)


# #################################################;
# I. Read the Google sheet containing the data info;
# #################################################;

# 1. upload the Google sheet: identify and exclude the sporadic samples;
updt_GoogleSht_6Dec2016 <- read.xlsx("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Nov16_Analysis_FTD/Workspace_Files/FTD_Brain_6Dec16_AftrDelRedonSam.xlsx", sheetIndex = 1, header = TRUE, stringsAsFactors=FALSE)
dim(updt_GoogleSht_6Dec2016)
#[1] 247  32;

# 2. check the #samples in each subgp of cases/conts;
table(updt_GoogleSht_6Dec2016$MUTATED.GENE)
#C9orf72 C9orf72_homozygote            control                GRN               MAPT                 NA 
#50                  3                 66                 41                 61                 26
sum(table(updt_GoogleSht_6Dec2016$MUTATED.GENE))
#[1] 247 # Note: so, that includes even the outliers;


# 2.1 identify which mutation category does the NA samples belong;
# @if NAs (mutation category) have disease code as sporadic, then refill it;
length(which(updt_GoogleSht_6Dec2016$MUTATED.GENE == "NA"))
#[1] 26;
length(grep("spo", updt_GoogleSht_6Dec2016$DISEASE.CODE[which(updt_GoogleSht_6Dec2016$MUTATED.GENE == "NA")], ignore.case = T))
#[1] 26;


# 2.2 since all NAs (as mutation category) has disease code as sporadic, substitute NA with sporadic;
updt_GoogleSht_6Dec2016$MUTATED.GENE[grep("spo", updt_GoogleSht_6Dec2016$DISEASE.CODE, ignore.case = T)] <- "Sporadic"

# @crosschecks; 

# -#mutation categories after substituting the NAs;
table(updt_GoogleSht_6Dec2016$MUTATED.GENE)
#C9orf72 C9orf72_homozygote            control                GRN               MAPT           Sporadic 
#50                  3                 66                 41                 61                 26 

# --for the substituted mutation catgories, check if both disease-code and the mutated gene has the same disease label as sporadic;
spo_x <- updt_GoogleSht_6Dec2016[grep("spo", updt_GoogleSht_6Dec2016$MUTATED.GENE, ignore.case = T), c("GIVENSAMPLENAME", "DISEASE.CODE", "MUTATED.GENE")]
dim(spo_x) 
# [1] 26  3;
unique(unlist(lapply(lapply(seq(nrow(spo_x)), function(y) grep("spo", spo_x[y,], ignore.case = T)), length)))
#[1] 2; NOTE: counting reports 2, meaning that both cols had "spo" as the label, meaning replacement was done correctly;



# #######################################################;
# II. Upload the count table generated using the pipeline;
# #######################################################;

# 1. the count table;
FronOnly_counTab <- read.table("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Nov16_Analysis_FTD/Server_Files/Files_after_05Dec2016/Count_Tables/CounTab_Fron_CTSS_files/countTable_Fron_CTSSfiles.txt")
dim(FronOnly_counTab)
#[1] 31000    58;


# 2. identify and exclude the sporadic samples;
FronOnly_counTab_colNams <- colnames(FronOnly_counTab)
length(FronOnly_counTab_colNams)
# [1] 58;

# 2.1) as a crosscheck, generate a dataframe of the original and the split colnames;
FronOnly_counTab_colNams_splt <- stri_split_fixed(stri_split_fixed(FronOnly_counTab_colNams, "sample_", simplify = T)[,2], "_no_chrM", simplify = T)[,1]
FronOnly_colNm_chkspt <- data.frame(ori_colNm = colnames(FronOnly_counTab), splt_colNm = FronOnly_counTab_colNams_splt, stringsAsFactors = F)
unique(unlist(lapply(lapply(seq(FronOnly_counTab_colNams_splt), function(x) grep(FronOnly_counTab_colNams_splt[x], FronOnly_colNm_chkspt[x,])), length)))
#[1] 2; NOTE: counting the split names at every rows result a length of 2 matches at every row, therefore, original and the split colnames matched;

# 2.2) reorder metadata in the order of the countable, plus exclude the unwanted columns;
updt_GoogleSht_FronOnly_mtchd <- updt_GoogleSht_6Dec2016[match(FronOnly_counTab_colNams_splt, updt_GoogleSht_6Dec2016$GIVENSAMPLENAME),
                                                         -(seq(grep("Pipeline_Passed", colnames(updt_GoogleSht_6Dec2016)), ncol(updt_GoogleSht_6Dec2016)))]
dim(updt_GoogleSht_FronOnly_mtchd)
#[1] 58 22;
unique(updt_GoogleSht_FronOnly_mtchd$REGION)
#[1] "frontal";

# 2.3) check if the above ordered metadata is correctly ordered as per the count table;
identical(updt_GoogleSht_FronOnly_mtchd$GIVENSAMPLENAME, FronOnly_counTab_colNams_splt)
# [1] TRUE; NOTE: correctly matched!;
identical(updt_GoogleSht_FronOnly_mtchd$GIVENSAMPLENAME, stri_split_fixed(stri_split_fixed(colnames(FronOnly_counTab), "sample_", simplify = T)[,2], "_no_chrM", simplify = T)[,1])
# [1] TRUE;

# 2.4) since both metadata and the count table are in the same order, suffix case/cont status for each column name in the count-table;
FronOnly_counTab_new <- FronOnly_counTab
dim(FronOnly_counTab_new)
#[1] 31000    58;
colnames(FronOnly_counTab_new) <- paste(colnames(FronOnly_counTab_new), stri_split_fixed(updt_GoogleSht_FronOnly_mtchd$MUTATED.GENE, "_", simplify = T)[,1], sep = "_")
all(FronOnly_counTab_new == FronOnly_counTab)
# [1] TRUE;

# 2.5) crosscheck the order of the renamed colname;
FronOnly_colNm_chkspt_fin <- data.frame(Googlsht_givSam=updt_GoogleSht_FronOnly_mtchd$GIVENSAMPLENAME, ori_colNm=colnames(FronOnly_counTab), newC_colNm=colnames(FronOnly_counTab_new), stringsAsFactors = F)
unique(unlist(lapply(lapply(seq(nrow(FronOnly_colNm_chkspt_fin)), function(x) grep(FronOnly_colNm_chkspt_fin$Googlsht_givSam[x], FronOnly_colNm_chkspt_fin[x,])), length)))
#[1] 3; NOTE: the unique length of the matched string is 3, therefore, the order is correct;

# 2.7) finally, exclude the sporadic samples from the countable;
FronOnly_counTab_minSpo <- FronOnly_counTab_new[,-(grep("Spor", colnames(FronOnly_counTab_new), ignore.case = T))]
dim(FronOnly_counTab_minSpo)
#[1] 31000    49;

all(FronOnly_counTab_minSpo[,intersect(colnames(FronOnly_counTab_minSpo), colnames(FronOnly_counTab_new))] 
    == FronOnly_counTab_new[,intersect(colnames(FronOnly_counTab_minSpo), colnames(FronOnly_counTab_new))])
# [1] TRUE;



# 3. obtain the start and end positions of the cluster ids;
FronOnly_counTab_minSpo_splitClusNames <- as.data.frame(matrix(unlist(strsplit(rownames(FronOnly_counTab_minSpo), "_")), ncol=4, byrow=T), stringsAsFactors = FALSE)
colnames(FronOnly_counTab_minSpo_splitClusNames) <- c("Chr", "Chr_Start", "Chr_End", "Strand")


# 4. transform the column classes;
FronOnly_counTab_minSpo_splitClusNames_df <- transform(FronOnly_counTab_minSpo_splitClusNames, 
                                                       "Chr" = as.character(Chr), 
                                                       Chr_Start = as.numeric(Chr_Start), 
                                                       Chr_End = as.numeric(Chr_End), 
                                                       Strand = as.factor(Strand))
str(FronOnly_counTab_minSpo_splitClusNames_df)
#'data.frame':	31000 obs. of  4 variables:
# $ Chr      : chr  "chr1" "chr1" "chr1" "chr1" ...
# $ Chr_Start: num  629906 630570 631352 633992 778791 ...
# $ Chr_End  : num  629940 630575 631382 634030 778809 ...
# $ Strand   : Factor w/ 2 levels "-","+": 2 2 2 2 2 2 2 2 2 2 ...

# 5. Attach a "score" column to the dataframe;
FronOnly_counTab_minSpo_splitClusNames_df$Clus_ID <- rownames(FronOnly_counTab_minSpo)
identical(paste(FronOnly_counTab_minSpo_splitClusNames_df$Chr, FronOnly_counTab_minSpo_splitClusNames_df$Chr_Start, FronOnly_counTab_minSpo_splitClusNames_df$Chr_End, FronOnly_counTab_minSpo_splitClusNames_df$Strand, sep = "_"), 
          rownames(FronOnly_counTab_minSpo)) 
# TRUE;

FronOnly_counTab_minSpo_splitClusNames_df$Score <- 0
colnames(FronOnly_counTab_minSpo_splitClusNames_df)[1] <- "#Chr"

FronOnly_counTab_minSpo_splitClusNames_df_fin <- FronOnly_counTab_minSpo_splitClusNames_df[, c("#Chr", "Chr_Start", "Chr_End", "Clus_ID", "Score", "Strand")]
dim(FronOnly_counTab_minSpo_splitClusNames_df_fin)
#[1] 31000     6;

# 5. Now, cbind with the count data;
FronOnly_counTab_minSpo_splitClusNames_final <- cbind(FronOnly_counTab_minSpo_splitClusNames_df_fin, FronOnly_counTab_minSpo)
dim(FronOnly_counTab_minSpo_splitClusNames_final)
#[1] 31000    55;

# @check if the splitted clusid is similar to the non-splitted ones;
identical(paste(FronOnly_counTab_minSpo_splitClusNames_final[,"#Chr"], FronOnly_counTab_minSpo_splitClusNames_final$Chr_Start, FronOnly_counTab_minSpo_splitClusNames_final$Chr_End, FronOnly_counTab_minSpo_splitClusNames_final$Strand, sep="_"), 
          FronOnly_counTab_minSpo_splitClusNames_final$Clus_ID)
# [1] TRUE;
identical(rownames(FronOnly_counTab_minSpo_splitClusNames_final), FronOnly_counTab_minSpo_splitClusNames_final$Clus_ID)
# [1] TRUE;

# 6. Now, save the count table as a BED file;
write.table(FronOnly_counTab_minSpo_splitClusNames_final, "/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Frontal_RikAnno/genBEDfile/FronOnly_counTab_minSpo_final.bed", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)



# ########################;
# III. Save this workspace;
# ########################;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc")
save.image("S4_Gen_BEDfiles_FronSamps.RData")
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc/S4_Gen_BEDfiles_FronSamps.RData")
# End ______________________________________________________________________________________________________________________________________________________________________;
