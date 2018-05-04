# Date: 2 May, 2018;
# Project: RIMOD-FTD;
# Aim: to perform hypergeometric test to identify over-represented GO-BP-MF-CC in all the interesting modules; 
# -(WGCNA) module: significantly correlated to the bio-trait, also those containing key genes;
# NOTE: the annotation of the count-table was carried using RIKEN's CAT gtf file as the gene model;
# NOTE 2: the same script could be adapted for the count table from: all 7 brain regions; and, temporal brain; 
# Workflow (main steps);
# Step 6: Unsigned - functional annotation using WGCNA library;


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
library(GSEABase)



# [I] create a list of gene-set collections ###############################################################################################; 

# *GO;
GO_BP <- GOGeneSets(species="Hs", ontologies=c("BP"))
GO_MF <- GOGeneSets(species="Hs", ontologies=c("MF"))
GO_CC <- GOGeneSets(species="Hs", ontologies=c("CC"))

# *KEGG;
PW_KEGG <- KeggGeneSets(species="Hs")

# *MSigDB;
setwd("/home/student-ag-javier/Documents/Project-FTD/R_CAGE/HTSanalyzeR/MSigDb")

# -a) genesets from the curated (C2) MSigDB all;
MSigDB_c2_all <- getGmt(con="c2.all.v6.1.symbols.gmt", geneIdType = SymbolIdentifier(), collectionType= BroadCollection(category="c2"))
MSigDB_c2_all_entrez <- mapIdentifiers(MSigDB_c2_all, EntrezIdentifier("org.Hs.eg.db"))
MSigDB_c2_all_entrez_gs <- geneIds(MSigDB_c2_all_entrez)
names(MSigDB_c2_all_entrez_gs) <- names(MSigDB_c2_all_entrez)

# -b) genesets from the MSigDB canonical all;
MSigDB_c2_CanonAll <- getGmt(con="c2.cp.v6.1.symbols.gmt", geneIdType = SymbolIdentifier(), collectionType= BroadCollection(category="c2"))
MSigDB_c2_CanonAll_entrez <- mapIdentifiers(MSigDB_c2_CanonAll, EntrezIdentifier("org.Hs.eg.db"))
MSigDB_c2_CanonAll_entrez_gs <- geneIds(MSigDB_c2_CanonAll_entrez)
names(MSigDB_c2_CanonAll_entrez_gs) <- names(MSigDB_c2_CanonAll_entrez)

# @organize the GSls into a common list;
ListGSC <- list(GO_BP=GO_BP, GO_MF=GO_MF, GO_CC=GO_CC, PW_KEGG=PW_KEGG, MSigAll = MSigDB_c2_all_entrez_gs, MSigCanonAll = MSigDB_c2_CanonAll_entrez_gs)



# [II] load the workspace from All cases combined vs. control samples ###################################################################################################################;
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/R_svdWrkSpc_fro/Fron_AllCasesVsCont_RikAnn_SVA_DEG.RData")



# [III] Prepare the universe gene list ###############################################################################################################################;

# 1. Extract all the expressed genes;
Fron_AlCasesVsCont_SVA_toptags_All <- topTags(Fron_AlCasesVsCont_SVA_glmQlFtest, n=Inf, adjust.method = "BH", sort.by = "PValue", p.value = 1)
dim(Fron_AlCasesVsCont_SVA_toptags_All)
# [1] 24459    28;

# 2. select only the protein coding and the promoters;
Fron_AlCasesVsCont_SVA_toptags_All_mRNA <- Fron_AlCasesVsCont_SVA_toptags_All$table[grep("coding_mRNA", Fron_AlCasesVsCont_SVA_toptags_All$table$CAT_geneClass, ignore.case = T),]
dim(Fron_AlCasesVsCont_SVA_toptags_All_mRNA)
# [1] 21418    28;
identical(grep("coding_mRNA", Fron_AlCasesVsCont_SVA_toptags_All$table$CAT_geneClass, ignore.case = T), which(Fron_AlCasesVsCont_SVA_toptags_All$table$CAT_geneClass == "coding_mRNA"))
# [1] TRUE;

# 3. now, select only the promoters;
Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt <- Fron_AlCasesVsCont_SVA_toptags_All_mRNA[grep("Promoter", Fron_AlCasesVsCont_SVA_toptags_All_mRNA$annotation, ignore.case = T),]
dim(Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt)
# [1] 21132    28;

# @Just a check;
unique(Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt$CAT_geneClass)
# [1] "coding_mRNA"
unique(Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt$annotation)
# [1] "Promoter (<=1kb)" "Promoter (2-3kb)" "Promoter (1-2kb)"

# @housekeeping;
which.na(Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt$entrez_ID)
# integer(0);
unique(Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt$entrez_ID[grep("na", Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt$entrez_ID)])
# [1] "__na";
length(Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt$entrez_ID[grep("na", Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt$entrez_ID)])
# [1] 351;

# 4. Cleaning: since there are na entrez, have to remove them;
Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt_cln <- Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt[-(which(Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt$entrez_ID == "__na")),]
dim(Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt_cln)
# [1] 20781    28;
identical(which(Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt$entrez_ID == "__na"), grep("na", Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt$entrez_ID))
# [1] TRUE;

# 5. Finally, select the universe genelist with gene-ids as entrez and values as logfc;
Fron_AlCasesVsCont_Univ <- Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt_cln$logFC
names(Fron_AlCasesVsCont_Univ) <- Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt_cln$entrez_ID
length(Fron_AlCasesVsCont_Univ)
# [1] 20781;

# @crosschecks;
identical(as.numeric(Fron_AlCasesVsCont_Univ), as.numeric(Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt_cln$logFC))
# [1] TRUE;
identical(as.character(names(Fron_AlCasesVsCont_Univ)), as.character(Fron_AlCasesVsCont_SVA_toptags_All_mRNA_promt_cln$entrez_ID))
# [1] TRUE;
grep("na", names(Fron_AlCasesVsCont_Univ), ignore.case = T)
# integer(0);



# [IV] Prepare the hits gene list ###############################################################################################################################;
#
# _@load the workspace that analyzed calculated GS and MM;
# #load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc/WGCNA_S4_FronSamps_AllCasesNConts_CAGEseq.RData")
#
# __@call all the functions I designed to get values specific for each cluster IDs repestively for each module;
# source("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_Scripts/TenzinNyima_functions4WGCNA_analysis.R")

# 1. combine all the modules that are interesting;
Fro_intMod_Nms <- c(Fro_bothStats_FDRsigMEnm, Fro_KeyFTDGenes_asocModNms)  
length(Fro_intMod_Nms)  
# [1] 15;

# 2. Put the dataframes of the modules in a list;
Fro_olCasesNcont_sigMod_DFls <- Fro_KeyFTDGenes_asocMod_AnnoDF(ModNames = Fro_intMod_Nms, ModAnnoDF = Fro_alCasesVsCont_WGCNA_GeneInf3)
length(Fro_olCasesNcont_sigMod_DFls)
# [1] 15;

# @check the module sizes;
unlist(lapply(Fro_olCasesNcont_sigMod_DFls, nrow))
# skyblue3; lightyellow; blue; lightgreen; greenyellow; purple; pink; green; orange; red; turquoise; brown; black; darkmagenta; royalblue;
# 91        281          5475  187         316          326     363   651    798     566  6393       918    480    101          178;

# 2.1 use the funtion(exEntrezIds_DFinLs_updt) to extract the entrez-ids only but without tha "__na's";

# 2.2 extract the entrez ids;
Fro_olCasesNcont_sigMod_entrezls <- exEntrezIds_DFinLs_updt(ModAnnoDFinLs = Fro_olCasesNcont_sigMod_DFls)
length(Fro_olCasesNcont_sigMod_entrezls)
# [1] 15;

# @crosschecks;
identical(names(Fro_olCasesNcont_sigMod_DFls), names(Fro_olCasesNcont_sigMod_entrezls))
# [1] TRUE;
unique(unlist(lapply(seq(Fro_olCasesNcont_sigMod_entrezls), function(x)
  identical(Fro_olCasesNcont_sigMod_entrezls[[x]],
            lapply(seq(Fro_olCasesNcont_sigMod_DFls), function(y)
              if(length(which(Fro_olCasesNcont_sigMod_DFls[[y]]$entrez_ID == "__na")) == 0){
                Fro_olCasesNcont_sigMod_DFls[[y]]$entrez_ID
              } else {Fro_olCasesNcont_sigMod_DFls[[y]]$entrez_ID[-(which(Fro_olCasesNcont_sigMod_DFls[[y]]$entrez_ID == "__na"))]} ) [[x]] ) )))
# [1]  TRUE;



# [V] Perform hypergeometric test ###############################################################################################################################;

# 1. Initialize and pre-process;
Fro_olCasesNCont_gscals_pp <- vector("list", length = length(Fro_olCasesNcont_sigMod_entrezls))
names(Fro_olCasesNCont_gscals_pp) <- names(Fro_olCasesNcont_sigMod_entrezls)

# @get a preprocessed ls;
for(i in seq(length(Fro_olCasesNCont_gscals_pp))){
  Fro_olCasesNCont_gscals_pp[[i]] <- new("GSCA", listOfGeneSetCollections=ListGSC, geneList=Fron_AlCasesVsCont_Univ, hits=Fro_olCasesNcont_sigMod_entrezls[[i]])  
  Fro_olCasesNCont_gscals_pp[[i]] <- preprocess(Fro_olCasesNCont_gscals_pp[[i]], species="Hs", initialIDs="Entrez.gene", keepMultipleMappings=TRUE, duplicateRemoverMethod="max", orderAbsValue=FALSE)
}

# 2. now, perform the analysis;
Fro_olCasesNCont_gscals_pp_ana_annot <- vector("list", length = length(Fro_olCasesNCont_gscals_pp))
names(Fro_olCasesNCont_gscals_pp_ana_annot) <- names(Fro_olCasesNCont_gscals_pp)

# @get the analyzed and annotated ls;
for(i in seq(Fro_olCasesNCont_gscals_pp_ana_annot)){
  Fro_olCasesNCont_gscals_pp_ana_annot[[i]] <- analyze(Fro_olCasesNCont_gscals_pp[[i]], para=list(pValueCutoff=0.05, pAdjustMethod="BH", minGeneSetSize=10), doGSOA=TRUE, doGSEA=FALSE)
  Fro_olCasesNCont_gscals_pp_ana_annot[[i]] <- appendGSTerms(Fro_olCasesNCont_gscals_pp_ana_annot[[i]], keggGSCs="PW_KEGG", goGSCs=c("GO_BP", "GO_MF", "GO_CC"))
}

# 3. report the results in a html;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Fron_WGCNA/Unsigned/All_CasesVSCont/Step6_HTSanalyzeR_Results")

# @initialize exp names;
Fro_AllCasesNCont_ExpNms <- paste("newPara_Fro_AllCasesNCont_", names(Fro_olCasesNCont_gscals_pp_ana_annot), sep = "")

# @report the analysis;
for(i in seq(Fro_AllCasesNCont_ExpNms)){
  report(Fro_olCasesNCont_gscals_pp_ana_annot[[i]], 
         experimentName=Fro_AllCasesNCont_ExpNms[i], 
         species="Hs", allSig=TRUE, keggGSCs="PW_KEGG", goGSCs=c("GO_BP", "GO_MF", "GO_CC"), 
         reportDir=paste(Fro_AllCasesNCont_ExpNms[i], "_HypGeom", sep = ""))  
}



# [VI] Report the #Sig processes ###############################################################################################################################;
# to check the summary;
# summarize(Fro_olCasesNCont_gscals_pp_ana_annot[[1]])
# attributes(Fro_olCasesNCont_gscals_pp_ana_annot$lightyellow)$result$HyperGeo.results$MSigCanonAll;

# 1. use the function(GSOA_sigRes_lsTOSumm) to report the output - #significant processes identified in each module;

# 2. use the above function to report the #sig process;
Fro_olCasesNcont_sigMod_gscaSumm <- GSOA_sigRes_lsTOSumm(Fro_olCasesNCont_gscals_pp_ana_annot)
dim(Fro_olCasesNcont_sigMod_gscaSumm)
# [1] 15 6;

# @include a column specifying the size of each module;
identical(names(Fro_olCasesNcont_sigMod_DFls), rownames(Fro_olCasesNcont_sigMod_gscaSumm))
# [1] TRUE;
Fro_olCasesNcont_sigMod_gscaSummUpdt <- cbind(Size = unlist(lapply(Fro_olCasesNcont_sigMod_DFls, nrow)), Fro_olCasesNcont_sigMod_gscaSumm)
dim(Fro_olCasesNcont_sigMod_gscaSummUpdt)
# [1] 15  7;



# [V] Save this workspace ###############################################################################################################################;
setwd("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc")
save.image("WGCNA_S6_FronSamps_AllCasesNConts_CAGEseq.RData")
# load("/home/student-ag-javier/Documents/Project-FTD/Correct_RefGenome_Res_Count_Data/Analy_4m_20Jan17_RikenGeneModel/Scripts/Fron_WGCNA_RWrkSpc/WGCNA_S6_FronSamps_AllCasesNConts_CAGEseq.RData")
# END ____________________________________________________________________________________________________________________________________________________________________________________________;
