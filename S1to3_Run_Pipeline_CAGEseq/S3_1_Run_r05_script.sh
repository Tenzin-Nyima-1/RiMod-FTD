#!/bin/bash
#SBATCH -A nyimat
#SBATCH -J FTD
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tenzin.nyima@dzne.de
#SBATCH --export=ALL
#SBATCH --cpus-per-task 64
#SBATCH -N 1
#SBATCH --nodelist=compute-01

module load R

srun /home/nyimat/cage_pipeline_scripts/R_scripts/r05_make_count_table_CAGEr.R BSgenome.Hsapiens.UCSC.hg38 /home/nyimat/All_ctss_files_no_chrM/AllCTSS_4m_7Regions/

wait
