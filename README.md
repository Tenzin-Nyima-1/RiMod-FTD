# RiMod-FTD

**Goal:** to investigate the differences and similarities between sub-groups of patients with frontotemporal dementia (FTD), and their differences with the control group. Such differences were investigated across brain regions that are most affected by FTD. 

**Method:** the data was processed in BASH using the in-house pipeline. The analysis was performed using the [R](https://cran.r-project.org/) and [Bioconductor](http://bioconductor.org/) software tools.
1. Model genome-wide count data (CAGEseq) as a function of case or control status
2. Identified modules containing highly correlated genes using the weighted gene coexpression network analysis
3. Functional annotation of the hits (differentially expressed genes or modules containing highly correlated genes) was performed using the hypergeometric test.

_**NOTE**_
1. The first three steps of the analysis generate the count table using the in-house pipeline. The relevant scripts are in the folder, "S1to3_Run_Pipeline_CAGEseq". 

1. Then the annotation of the count table is performed using the scripts in the folder, "S4to5_Annotate_CounTable_CAGEseq".


**Used the following R packages:**
* [base](https://cran.r-project.org/bin/windows/base/old/3.2.3/)
* [CAGEr](http://bioconductor.org/packages/release/bioc/html/CAGEr.html)
* [stringr](https://cran.r-project.org/web/packages/stringr/index.html)
* [stringi](https://cran.r-project.org/web/packages/stringi/index.html)
* [openxlsx](https://cran.r-project.org/web/packages/openxlsx/index.html)
* [xlsx](https://cran.r-project.org/web/packages/xlsx/index.html)
* [muStat](https://cran.r-project.org/web/packages/muStat/index.html)
* [GenomicFeatures](http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
* [org.Hs.eg.db](http://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)
* [ChIPseeker](http://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)

_(Note: this repository is not complete yet. Additional scripts relevant to the downstream analysis will be included shortly into this repository)_
