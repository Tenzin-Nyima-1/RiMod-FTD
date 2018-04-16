# Date: 24 January, 2017;
# Project: RiMod-FTD;
# Aim: to annotate the frontal count table (BED file format);
# NOTE: the same script could be adapted for the count table from: all 7 brain regions; and, temporal brain; 
# Workflow;
# Step 1: make a TxDb object from the FANTOM-CAT-GTF file;
# Step 2: annotation of the count tablefrom the frontal-49-samples using the TxDb;


library(GenomicFeatures)
library(org.Hs.eg.db)


################################################################################################;
# [I] make TxDb object from the GTF file________________________________________________________;
################################################################################################;

# 1. locate the gtf file;
fantom6_inpDir <- "/home/student-ag-javier/Documents/Fantom_6/Datas/F6_CAT_GTFfile/"
fantom6_inp_ls_files <- list.files(fantom6_inpDir, full.names = TRUE)
# [1] "/home/student-ag-javier/Documents/Fantom_6/Datas/F6_CAT_GTFfile//F6_CAT_GeneInf_RecvOn12Jan17" 
# [2] "/home/student-ag-javier/Documents/Fantom_6/Datas/F6_CAT_GTFfile//F6_CAT.transcript.gtf" #*;
# [3] "/home/student-ag-javier/Documents/Fantom_6/Datas/F6_CAT_GTFfile//F6_CAT.transcript.gtf.gz"     
# [4] "/home/student-ag-javier/Documents/Fantom_6/Datas/F6_CAT_GTFfile//F6_old_new_chck.txt"         
# [5] "/home/student-ag-javier/Documents/Fantom_6/Datas/F6_CAT_GTFfile//New_check_gtfFile"       

#2. convert the GTF into TXDb object;
fantom6_TxDb <- makeTxDbFromGFF(fantom6_inp_ls_files[2], format = "gtf", dataSource = "FANTOM6",  organism = "Homo sapiens")

#3. Got the following warnings;
# Import genomic features from the file as a GRanges object ... OK
#Prepare the 'metadata' data frame ... OK
#Make the TxDb object ... OK
#Warning messages:
#1: RSQLite::dbGetPreparedQuery() is deprecated, please switch to DBI::dbGetQuery(params = bind.data). 
#2: Named parameters not used in query: internal_chrom_id, chrom, length, is_circular 
# 3: Named parameters not used in query: internal_id, name, type, chrom, strand, start, end 
# 4: Named parameters not used in query: internal_id, name, chrom, strand, start, end 
# 5: Named parameters not used in query: internal_tx_id, exon_rank, internal_exon_id, internal_cds_id 
# 6: Named parameters not used in query: gene_id, internal_tx_id 
#
# However such warnings were also reported on this bioconductor help file as well;
# https://bioconductor.org/packages/devel/bioc/vignettes/AnnotationHub/inst/doc/AnnotationHub-HOWTO.html;



################################################################################################;
# [II] annotation of the count table using the above TxDb_______________________________________;
################################################################################################;
library(ChIPseeker)
library(stringi)

#1. Read bed file;
BED_FronOnly_inpDir <- "/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Frontal_RikAnno/genBEDfile"
BED_FronOnly_lsFiles <- list.files(BED_FronOnly_inpDir, full.names = TRUE)
#[1] "/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Frontal_RikAnno/genBEDfile/FronOnly_counTab_minSpo_final.bed";

#2. Read the peakfile file;
FronOnly_peak <- readPeakFile(BED_FronOnly_lsFiles[1])

#3. read genomic ranges;
FronOnly_peak_GR <- GRanges(FronOnly_peak, strand = Rle(strand(read.table(BED_FronOnly_lsFiles[1])[,6])))

#4. define promoters around 3K(since most reads are located around that region ~90%);
promoter_3k <- getPromoters(TxDb=fantom6_TxDb, upstream=3000, downstream=3000, by="transcript")
length(promoter_3k)
#[1] 379952;

#4.1 get a tag-matrix around that promoter region;
FronOnly_peak_GR_tagMat_3k <- getTagMatrix(FronOnly_peak_GR, windows=promoter_3k)
dim(FronOnly_peak_GR_tagMat_3k)
#[1] 90994  6001;

#4.2 visualize them using a heatmap;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Frontal_RikAnno/annotation")

pdf("FronOnly_tagMat_heatmap_3k.pdf")
tagHeatmap(FronOnly_peak_GR_tagMat_3k,  xlim=c(-3000, 3000), xlab="Genomic range", ylab="Transcripts", color="red")
dev.off()

#4.3 visualize them using average profile;
pdf("FronOnly_tagMat_AvgProf_3k.pdf")
plotAvgProf(FronOnly_peak_GR_tagMat_3k,  xlim=c(-3000, 3000), xlab="Genomic range", ylab="Read Count Frequency")
dev.off()
# NOTE: almost 95.5% of the reads are within 2kb around TSS, and only a small 1.07% is between 2-3kb, also read count is not close to 0 at 3kb hence this range was chosen;


#5. Annotation within the range of +/- 3kbases;
FronOnly_peak_GR_3k <- annotatePeak(FronOnly_peak_GR, tssRegion=c(-3000, 3000), TxDb=fantom6_TxDb, level = "transcript", assignGenomicAnnotation = TRUE, genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"), annoDb="org.Hs.eg.db", addFlankGeneInfo = FALSE, flankDistance = 5000, verbose = TRUE)

#5.1 Generate summary plots for the annotation performed;
pdf("AnnoPie_FronOnly_3k.pdf")
plotAnnoPie(FronOnly_peak_GR_3k)
dev.off()
pdf("AnnoBar_FronOnly_3k.pdf")
plotAnnoBar(FronOnly_peak_GR_3k)
dev.off()
pdf("UpsetPlot_FronOnly_3k.pdf")
upsetplot(FronOnly_peak_GR_3k)
dev.off()
pdf("Vennpie_FronOnly_3k.pdf")
vennpie(FronOnly_peak_GR_3k)
dev.off()
pdf("Dist2TSS_FronOnly_3k.pdf")
plotDistToTSS(FronOnly_peak_GR_3k, title="Distribution of CAGE clusters relative to TSS")
dev.off()


#6. convert the annotated file in s4 class into a dataframe, to save as table;
FronOnly_peak_GR_3k_DF <- as.data.frame(FronOnly_peak_GR_3k)
rownames(FronOnly_peak_GR_3k_DF) <- as.character(FronOnly_peak_GR_3k_DF[,"V4"])
identical(rownames(FronOnly_peak_GR_3k_DF), as.character(FronOnly_peak_GR_3k_DF[,"V4"]))
#[1] TRUE;

#7. to identify colnames of the counts, upload the saved workspace used for generating bed file;
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc/S4_Gen_BEDfiles_FronSamps.RData")
# keep only the relevant R objects;
rm(list = ls()[grep("All7RegionSamps|TempOnly", ls())])

# check if the columns are same;
all(FronOnly_peak_GR_3k_DF[,grep("V", colnames(FronOnly_peak_GR_3k_DF))] == 
      FronOnly_counTab_minSpo_splitClusNames_final[,seq(grep("Clus_ID", colnames(FronOnly_counTab_minSpo_splitClusNames_final)), ncol(FronOnly_counTab_minSpo_splitClusNames_final))])
#[1] TRUE;

#8. since the count-tables matched, assign the colnames;
colnames(FronOnly_peak_GR_3k_DF)[grep("V", colnames(FronOnly_peak_GR_3k_DF))] <- 
  colnames(FronOnly_counTab_minSpo_splitClusNames_final[,seq(grep("Clus_ID", colnames(FronOnly_counTab_minSpo_splitClusNames_final)), ncol(FronOnly_counTab_minSpo_splitClusNames_final))])
dim(FronOnly_peak_GR_3k_DF)
# [1] 31000    66;

#9. organize by rearranging the columns;
FronOnly_peak_GR_3k_DF_ord <- FronOnly_peak_GR_3k_DF[, c(colnames(FronOnly_peak_GR_3k_DF)[-(grep("sample", colnames(FronOnly_peak_GR_3k_DF)))], colnames(FronOnly_peak_GR_3k_DF)[grep("sample", colnames(FronOnly_peak_GR_3k_DF))])]

# check if the columns are dupplicated, if yes, keep just one of them;
all(as.character(FronOnly_peak_GR_3k_DF_ord$seqnames) == as.character(FronOnly_peak_GR_3k_DF_ord$geneChr))
#[1] TRUE;
all(as.character(FronOnly_peak_GR_3k_DF_ord$strand) == as.character(FronOnly_peak_GR_3k_DF_ord$geneStrand))
#[1] TRUE;
all(as.character(FronOnly_peak_GR_3k_DF_ord$Strand) == as.character(FronOnly_peak_GR_3k_DF_ord$geneStrand))
#[1] TRUE;

#10. final reorder: exclude redundant annotation columns;
FronOnly_reqOrdCols <- c("Clus_ID", "width", "annotation", "distanceToTSS", "geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand", "geneId", "transcriptId", "control", "MAPT", "C9orf72", "GRN")  
FronOnly_peak_GR_3k_DF_ord_fin <- FronOnly_peak_GR_3k_DF_ord[,unlist(lapply(seq(FronOnly_reqOrdCols), function(x) grep(FronOnly_reqOrdCols[x], colnames(FronOnly_peak_GR_3k_DF_ord))))]
dim(FronOnly_peak_GR_3k_DF_ord_fin)
# [1] 31000    60;


#11. Crosschecks;

#@ within the annotated count-table;
identical(rownames(FronOnly_peak_GR_3k_DF_ord_fin), as.character(FronOnly_peak_GR_3k_DF_ord_fin$Clus_ID))
#[1] TRUE;

#@ count tables;
identical(FronOnly_peak_GR_3k_DF_ord[,grep("sample", colnames(FronOnly_peak_GR_3k_DF_ord))], FronOnly_counTab_minSpo)
#[1] TRUE;
identical(rownames(FronOnly_peak_GR_3k_DF_ord), rownames(FronOnly_counTab_minSpo))
#[1] TRUE;

#@ properly ordered count tables;
FronOnly_reqord_mainCount <- c("control", "MAPT", "C9orf72", "GRN")
identical(FronOnly_counTab_minSpo[,unlist(lapply(seq(FronOnly_reqord_mainCount), function(x) grep(FronOnly_reqord_mainCount[x], colnames(FronOnly_counTab_minSpo))))], 
          FronOnly_peak_GR_3k_DF_ord_fin[,grep("sample", colnames(FronOnly_peak_GR_3k_DF_ord_fin))])
#[1] TRUE
identical(rownames(FronOnly_counTab_minSpo), rownames(FronOnly_peak_GR_3k_DF_ord_fin))
#[1] TRUE;


#12. write this table;
getwd()
#[1] "/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Frontal_RikAnno/annotation";
write.table(FronOnly_peak_GR_3k_DF_ord_fin, "Annot_FronOnly_peak_GR_3k_DF_ord_fin.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)


###################################################################;
# [III] Save this workspace _______________________________________;
###################################################################;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc/")
save.image("S5_FronSamps_Annobedfile_F6CAT.RData")
#load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc/S5_FronSamps_Annobedfile_F6CAT.RData")
# END ___________________________________________________________________________________________________________________________________________________________________________;
