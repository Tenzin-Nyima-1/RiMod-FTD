# RiMod-FTD

**Goal:** to investigate the differences and similarities between sub-groups of patients with frontotemporal dementia (FTD), and their differences with the control group. Such differences were investigated across brain regions that are most affected by FTD. 

**Method:** the data was processed in BASH using the in-house pipeline. The analysis was performed using the [R](https://cran.r-project.org/) and [Bioconductor](http://bioconductor.org/) software tools.
1. Model genome-wide count data (CAGEseq) as a function of case or control status
2. Identified modules containing highly correlated genes using the weighted gene coexpression network analysis
3. Functional annotation of the hits (differentially expressed genes or modules containing highly correlated genes) was performed using the hypergeometric test.

_**NOTE**_
1. The first three steps of the analysis generate the count table using the in-house pipeline. The relevant scripts are in the folder, ["S1to3_Run_Pipeline_CAGEseq"](https://github.com/Tenzin-Nyima-1/RiMod-FTD/tree/master/S1to3_Run_Pipeline_CAGEseq). 

1. Then the annotation of the count table is performed using the scripts in the folder, ["S4to5_Annotate_CounTable_CAGEseq"](https://github.com/Tenzin-Nyima-1/RiMod-FTD/tree/master/S4to5_Annotate_CounTable_CAGEseq).

1. The scripts relevant to the differential gene expression (DGE) analysis and the pre-DGE data preparation/crosscheck is here: [" S6to7_preDE_DGEanaly_CAGEseq"](https://github.com/Tenzin-Nyima-1/RiMod-FTD/tree/master/S6to7_preDE_DGEanaly_CAGEseq).

**Used the following R packages:**
**Used the following R packages:**
* [base](https://cran.r-project.org/bin/windows/base/old/3.2.3/)
* [CAGEr](http://bioconductor.org/packages/release/bioc/html/CAGEr.html)
* [stringr](https://cran.r-project.org/web/packages/stringr/index.html); [stringi](https://cran.r-project.org/web/packages/stringi/index.html); [gdata](https://cran.r-project.org/web/packages/gdata/index.html)
* [openxlsx](https://cran.r-project.org/web/packages/openxlsx/index.html); [xlsx](https://cran.r-project.org/web/packages/xlsx/index.html); [readxl](https://cran.r-project.org/web/packages/readxl/index.html)
* [GenomicFeatures](http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html); [org.Hs.eg.db](http://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html); [BSgenome.Hsapiens.UCSC.hg38](http://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html); [TxDb.Hsapiens.UCSC.hg38.knownGene](http://bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg38.knownGene.html); [AnnotationDbi](http://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html); [GO.db](http://bioconductor.org/packages/release/data/annotation/html/GO.db.html); [Biobase](http://bioconductor.org/packages/release/bioc/html/Biobase.html); [KEGG.db](http://bioconductor.org/packages/release/data/annotation/html/KEGG.db.html)   
* [ChIPseeker](http://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
* [S4Vectors](https://bioconductor.org/packages/release/bioc/html/S4Vectors.html); [Hmisc](https://cran.r-project.org/web/packages/Hmisc/index.html); [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html); [muStat](https://cran.r-project.org/web/packages/muStat/index.html)
* [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html); [sva](http://bioconductor.org/packages/release/bioc/html/sva.html); [statmod](CRAN.r-project.org/web/packages/rpart/index.html); [HTSanalyzeR](http://bioconductor.org/packages/release/bioc/html/HTSanalyzeR.html); [snow](https://cran.r-project.org/web/packages/snow/index.html)
* [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html); [pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html); [plot3D](https://cran.r-project.org/web/packages/plot3D/index.html); [plotrix](https://cran.r-project.org/web/packages/plotrix/index.html)

_(Note: this repository is not complete yet. Additional scripts relevant to the downstream analysis will be included shortly into this repository)_
