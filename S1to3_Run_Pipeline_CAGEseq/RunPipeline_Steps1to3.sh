# Date: 21 March 2018;
# Project: RiMOD-FTD;
# Data: CAGEseq transcriptomic data;
# Aim: generate count table using in-house pipeline;
# NOTE: sample script for the first few steps of analyses;
# NOTE 2: should specify the right directory where the pipeline is kept;
# Workflow;
# Step 1-1: check the integrity of files;
# Step 1-2: generate barcodes file per sample in a pool using this script: "S1_2_GenBarcodes_files.Rscript";
# Step 1-2: get the CAGE-artefact file into the folder containing input files: fastq and barcodes txt file;
# NOTE: now the input dir contains: a fastq file, barcode txt file and the CAGE artefact file;
# Step 2-1: run the pipeline to generate CTSS files;
# Step 2-2: if required remove chrM counts from the CTSS files using this script: rmChrM_CTSS.sh;
# Step 3-1: finally, generate count table;




# [I] Step 1: prepration of the files required for running the pipeline #############################;
# the three required files(input) are;
# a) a fastq file;
# b) a txt file specifying barcodes coresponding to each sample in a pool;
# c) a cage-artefact file;

# 1. check the integrity of fastq files;
# @assuming fastq files are in the current working directory;
md5sum -c md5sums >>chk_md5sums.txt

# 2. generate a file specifying barcodes coresponding to each sample in a pool;
# - using this script "S1_2_GenBarcodes_files.Rscript";
R CMD BATCH S1_2_GenBarcodes_files.Rscript

# 3. get the CAGE-artefact file into the dir containing the input files, if it doesn't exist in the pwd;
rsync -r -avz --progress --ignore-existing ./input/CAGE_artefacts_201507.txt ./Pools/
# _____________________________________________________________________________________________________;




# [II] Step 2: run the pipeline to generate CTSS files #################################################;
# Note: the following script assumes that all the three required files (in step 1) are in the dir;

# 1. run the pipeline to generate CTSS files;
nohup ./cage_pipeline_scripts/python_scripts/dataset_cage_final.py -d ./Pools/ -r hg38 &

# 2. if required, run this script to remove chrM counts from the CTSS files: avoids heteroplasmic bias;
./S2_2_rmChrM_CTSS.sh
# _____________________________________________________________________________________________________;




# [III] Step 3: finally, generate count table #########################################################;
# Note: this run assumes that all the CTSS (with or without chrM) files are in the specified dir;

# 1. generate count table using this script: ;
sbatch ./S3_1_Run_r05_script.sh
# _____________________________________________________________________________________________________;

# The end _____________________________________________________________________________________________;
