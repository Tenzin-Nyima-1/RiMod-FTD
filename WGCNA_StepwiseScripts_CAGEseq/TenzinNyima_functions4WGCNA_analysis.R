# Date: 13 April 2018;
# By: Tenzin Nyima;
# NOTE: contains all the functions I created to run the tasks relevant to WGCNA analysis;
# NOTE2: they were written to suit my requirements, you can change them to suit your requirements;
# For: WGCNA analysis;



# #########################################################################;
# Function 1: extracts the cluster-ids of each individual module separately;
# _________________________________________________________________________;
Extract_ClusIDs <- function(ModNames, ModAnnoDF){
  #S1: generate an empty ls with same length as #mods;
  exClusIDs <- vector("list", length = length(ModNames))
  names(exClusIDs) <- ModNames
  #S2: assign clus-ids to respective modules;
  for(i in seq(length(exClusIDs))){
    exClusIDs[[i]] <- 
      lapply(seq(length(exClusIDs)), function(x) as.character(with(ModAnnoDF, Clus_ID[which(MrgdMECols == names(exClusIDs)[x])])))[[i]]
  }
  return(exClusIDs)
}


# ###############################################################################;
# Function 2: extracts the module-values of cluster-ids in each individual module;
# _______________________________________________________________________________;
Extract_MM <- function(ModNames, ModAnnoDF){
  #S1: generate an empty ls with same length as #mods;
  extMM <- vector("list", length = length(ModNames))
  names(extMM) <- ModNames
  #S2: get MMvalue for each clus in a module;
  for(i in seq(length(extMM))){
    extMM[[i]] <- ModAnnoDF[which(ModAnnoDF$MrgdMECols == names(extMM)[i]), 
                            match(paste("MM", names(extMM)[i], sep = ""), colnames(ModAnnoDF))]
  }
  return(extMM)
}


# ########################################################################################################################;
# Function 3: a unique function that extracts GS-values of both traits: all cases combined and cont; and, each case n cont;
# ________________________________________________________________________________________________________________________;
Extract_GS_both <- function(ModNames, ModAnnoDF, TraitColNm){
  #S1: generate an empty ls with same length as #mods;
  extGS <- vector("list", length = length(ModNames))
  names(extGS) <- ModNames
  #S2: get GSvalue for the sig-module;
  for(i in seq(length(extGS))){
    extGS[[i]] <- ModAnnoDF[which(ModAnnoDF$MrgdMECols == names(extGS)[i]), TraitColNm]
  }
  return(extGS)
}


# #################################################################################;
# Function 4: extracts annotation data frame separately for each module in the list;
# _________________________________________________________________________________;
Fro_KeyFTDGenes_asocMod_AnnoDF <- function(ModNames, ModAnnoDF){
  #S1: generate an empty ls with same length as #mods;
  KeyFTDGenes_AnnoDF <- vector("list", length = length(ModNames))
  names(KeyFTDGenes_AnnoDF) <- ModNames
  #S2: get the required col numbers;
  KeyFTDGenes_colNo <- c(1:19,
                         match("GSStat_CaseCont", colnames(ModAnnoDF)),
                         match("GSStat_MutCat", colnames(ModAnnoDF)))
  KeyFTDGenes_MMcolNo <- match(paste("MM", ModNames, sep = ""), 
                               colnames(ModAnnoDF))
  #S3: get the associated data frame;
  for(i in seq(length(KeyFTDGenes_AnnoDF))){
    KeyFTDGenes_AnnoDF[[i]] <- lapply(seq(length(KeyFTDGenes_AnnoDF)), function(x) 
      ModAnnoDF[which(ModAnnoDF$MrgdMECols == names(KeyFTDGenes_AnnoDF)[i]), c(KeyFTDGenes_colNo, KeyFTDGenes_MMcolNo[i])] )[[i]]
  }
  return(KeyFTDGenes_AnnoDF)
}


# ##########################################################;
# Function 5: extracts the range of x (MM) and y (GS) values;
# __________________________________________________________;
maxXYthresh_perMod <- function(ModAnnoDF, GStrt1, GStrt2){
  maxXYthres <- vector("list", length = length(ModAnnoDF))
  names(maxXYthres) <- names(ModAnnoDF)
  for(i in seq(length(maxXYthres))){
    maxXYthres[[i]] <- lapply(seq(length(ModAnnoDF)), function(x) c(
      as.numeric(substr(max(abs(ModAnnoDF[[names(maxXYthres)[i]]][, GStrt1])), 1, 3)), 
      as.numeric(substr(max(abs(ModAnnoDF[[names(maxXYthres)[i]]][, GStrt2])), 1, 3)),
      as.numeric(substr(max(abs(ModAnnoDF[[names(maxXYthres)[i]]][, grep("MM", colnames(ModAnnoDF[[names(maxXYthres)[i]]]))])), 1, 3)) ))[[i]]
    names(maxXYthres[[i]]) <- c(GStrt1, GStrt2, paste("MM", names(maxXYthres)[i], sep = ""))
  }
  return(maxXYthres)
}


# ########################################################################################################################;
# Function 6: labels only the data points corresponding to the key genes and the hub genes only (***The combined function);
# ________________________________________________________________________________________________________________________;
NonHubKeyGens_convBlnk <- function(KeyGens_AnnDFls,  KeyGens_XYthres_ls, GStrt, keyGeNms){
  #S1: get the index of nonHubGenes;
  NonHub_Ind <- vector("list", length = length(KeyGens_AnnDFls))
  names(NonHub_Ind) <- names(KeyGens_AnnDFls)
  for(i in seq(length(NonHub_Ind))){
    NonHub_Ind[[i]] <- lapply(seq(length(NonHub_Ind)), function(x)
      which(abs(KeyGens_AnnDFls[[x]][,GStrt]) > c(KeyGens_XYthres_ls[[x]][GStrt] - 0.1) &
              abs(KeyGens_AnnDFls[[x]][,grep("MM", names(KeyGens_AnnDFls[[x]]))]) > c(KeyGens_XYthres_ls[[x]][grep("MM", names(KeyGens_XYthres_ls[[x]]))] - 0.1) ))[[i]]
  }
  #S2: get the index of key FTD Genes;
  NonKey_Ind <- vector("list", length = length(KeyGens_AnnDFls))
  names(NonKey_Ind) <- names(KeyGens_AnnDFls)
  for(i in seq(length(NonKey_Ind))){
    NonKey_Ind[[i]] <- lapply(seq(length(NonKey_Ind)), function(x)
      c(na.omit(match(toupper(keyGeNms), toupper(KeyGens_AnnDFls[[x]][,"geneName"])))) )[[i]] }
  #S3: combine the above indexes into one;
  both_ind <- vector("list", length = length(KeyGens_AnnDFls))
  names(both_ind) <- names(KeyGens_AnnDFls)
  if(identical(names(NonHub_Ind), names(NonKey_Ind)))  for(i in seq(length(both_ind))) {
    both_ind[[i]] <- lapply(seq(length(both_ind)), function(x) unique(c(NonHub_Ind[[x]], NonKey_Ind[[x]])))[[i]] } else 
      warning("mismatches in the names/length of hub/non-key geneIndx")
  #S4: substitute unintertings ones as blanks in the main DFls;
  subsKeyGens_AnnDFls <- KeyGens_AnnDFls 
  for(i in seq(length(subsKeyGens_AnnDFls))){
    subsKeyGens_AnnDFls[[names(both_ind)[i]]][setdiff(seq(nrow(subsKeyGens_AnnDFls[[i]])), both_ind[[i]]),"geneName"] <- "" }
  return(subsKeyGens_AnnDFls)
}


# ##################################################################################################;
# Function 7: re-write a function to extract annotation DF for each modules separately but in a list;
# __________________________________________________________________________________________________;
Fro_KeyFTDGenesasocMod_AnnoIntkDF <- function(ModNames, ModAnnoDF){
  #S1: generate an empty ls with same length as #mods;
  KeyFTDGenes_AnnoDF <- vector("list", length = length(ModNames))
  names(KeyFTDGenes_AnnoDF) <- ModNames
  #S2: get the associated data frame;
  for(i in seq(length(KeyFTDGenes_AnnoDF))){
    KeyFTDGenes_AnnoDF[[i]] <- lapply(seq(length(KeyFTDGenes_AnnoDF)), function(x)
      ModAnnoDF[which(ModAnnoDF$MrgdMECols == names(KeyFTDGenes_AnnoDF)[i]), ] )[[i]]
  }
  return(KeyFTDGenes_AnnoDF)
}


# #################################################################################;
# Function 8: labels only those genes with top 10% intramodular connectivity values;
# _________________________________________________________________________________;
NonHubGens_convBlnk <- function(KeyGens_AnnDFls, keyGeNms){
  #S1: add an empty col to col sps gp of labels;
  for(i in seq(length(KeyGens_AnnDFls))){
    KeyGens_AnnDFls[[i]]$LablCol <- NA
  }
  #S2: get the index of nonHubGenes;
  NonHub_Ind <- vector("list", length = length(KeyGens_AnnDFls))
  names(NonHub_Ind) <- names(KeyGens_AnnDFls)
  for(i in seq(length(NonHub_Ind))){
    NonHub_Ind[[i]] <- lapply(seq(length(NonHub_Ind)), function(x) 
      which(KeyGens_AnnDFls[[names(NonHub_Ind)[x]]][,"kWithin"] > quantile(KeyGens_AnnDFls[[names(NonHub_Ind)[x]]][,"kWithin"], seq(0, 1, 0.10))[["90%"]]) )[[i]]
  }
  #S3: substitute labcol column;
  for(i in seq(length(KeyGens_AnnDFls))){
    KeyGens_AnnDFls[[names(NonHub_Ind)[i]]][NonHub_Ind[[i]],"LablCol"] <- "mediumvioletred" }
  #S4: get the index of key FTD Genes;
  NonKey_Ind <- vector("list", length = length(KeyGens_AnnDFls))
  names(NonKey_Ind) <- names(KeyGens_AnnDFls)
  for(i in seq(length(NonKey_Ind))){
    NonKey_Ind[[i]] <- lapply(seq(length(NonKey_Ind)), function(x)
      c(na.omit(match(toupper(keyGeNms), toupper(KeyGens_AnnDFls[[x]][,"geneName"])))) )[[i]] }
  #S5: substitute labcol column;
  for(i in seq(length(KeyGens_AnnDFls))){
    KeyGens_AnnDFls[[names(NonKey_Ind)[i]]][NonKey_Ind[[i]],"LablCol"] <- "blue" }
  #S6: combine the above indexes into one;
  both_ind <- vector("list", length = length(KeyGens_AnnDFls))
  names(both_ind) <- names(KeyGens_AnnDFls)
  if(identical(names(NonHub_Ind), names(NonKey_Ind)))  for(i in seq(length(both_ind))) {
    both_ind[[i]] <- lapply(seq(length(both_ind)), function(x) unique(c(NonHub_Ind[[x]], NonKey_Ind[[x]])))[[i]] } else 
      warning("mismatches in the names/length of hub/non-key geneIndx")
  #S7: substitute unintertings ones as blanks in the main DFls;
  for(i in seq(length(KeyGens_AnnDFls))){
    KeyGens_AnnDFls[[names(both_ind)[i]]][setdiff(seq(nrow(KeyGens_AnnDFls[[i]])), both_ind[[i]]),"geneName"] <- "" }
  return(KeyGens_AnnDFls)
}


# ##############################################################################;
# Function 9: incase the key-gene color is not clear use this function to rename;
# ______________________________________________________________________________;
LabTxt_colChng <- function(KeyGens_genLab_AnnDFls, setCol, chgCol){
  colVec <- KeyGens_genLab_AnnDFls[[match(setCol, names(KeyGens_genLab_AnnDFls))]]$LablCol
  KeyGens_genLab_AnnDFls[[match(setCol, names(KeyGens_genLab_AnnDFls))]]$LablCol <- ifelse(colVec == setCol, chgCol, colVec)
  return(KeyGens_genLab_AnnDFls)
}


# ##################################################################;
# Function 10: extracts the entrez-ids only but without tha "__na's";
# __________________________________________________________________;
exEntrezIds_DFinLs_updt <- function(ModAnnoDFinLs){
  exEntrezIds <- vector("list", length = length(ModAnnoDFinLs))  
  names(exEntrezIds) <- unlist(lapply(seq(length(ModAnnoDFinLs)), function(x) unique(ModAnnoDFinLs[[x]]$MrgdMECols)))
  for(i in seq(length(exEntrezIds))){
    exEntrezIds[[i]] <- lapply(seq(ModAnnoDFinLs), function(x) ModAnnoDFinLs[[x]]$entrez_ID)[[i]]
    if(length(which(exEntrezIds[[i]] == "__na")) == 0){
      exEntrezIds[[i]] <- exEntrezIds[[i]]
    } else{
      exEntrezIds[[i]] <- exEntrezIds[[i]][-(which(exEntrezIds[[i]] == "__na"))]}
  }
  return(exEntrezIds)
}


# #################################################################################;
# Function 11: report the output - #significant processes identified in each module;
# _________________________________________________________________________________;
GSOA_sigRes_lsTOSumm <- function(aGSCAls){
  #S1: save the names of the overRep procs;
  procNms <- names(attributes(aGSCAls[[1]])$result$HyperGeo.results)
  #S2: generate an empty ls to save the #sig processes;
  sigOvRepProc <- vector("list", length = length(aGSCAls))
  names(sigOvRepProc) <- names(aGSCAls)
  for(i in seq(sigOvRepProc)){
    sigOvRepProc[[i]] <- unlist(lapply(seq(procNms), function(x) 
      length(which(attributes(aGSCAls[[i]])$result$HyperGeo.results[[procNms[x]]][,"Adjusted.Pvalue"] < 0.05)) )) }
  #S3: convert the outputs in ls into a dataframe;
  sigOvRepProc_df <- as.data.frame(Reduce(rbind, sigOvRepProc))
  rownames(sigOvRepProc_df) <- names(aGSCAls)
  colnames(sigOvRepProc_df) <- procNms
  return(sigOvRepProc_df)
}

# _______________________________________________________________________________________________________________________;
# END ___________________________________________________________________________________________________________________;
# _______________________________________________________________________________________________________________________;