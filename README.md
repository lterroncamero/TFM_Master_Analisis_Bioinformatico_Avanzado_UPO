<h1 align="center"> TFM Master Analisis Bioinformatico Avanzado UPO</h1>

<h2 align="center"> Unravelling the molecular mechanisms in the early steps of the implantation process<h6 align="justify">

### Alumna: Laura C. Terrón Camero

### Co-tutores: Eduardo Andrés-León, Signe Altmäe

#### Abstract ####

> ##### STUDY QUESTION: Which genes in the luminal epithelium (the first maternal cellular lining in the uterus to contact with embryo) interact in the molecular dialogue with the implanting embryo? How is this molecular embryo-endometrium crosstalk regulated in the case of embryonal chromosomal abnormalities?

> ##### SUMMARY ANSWER: A set of molecules are expressed on the endometrial lining and implantation-competent blastocyst that could interact on the embryo-endometrium crosstalk, and this molecular dialogue seems to be dysregulated in the case of embryonal chromosomal abnormalities such as trisomies in chromosomes 7, 11, 15, 21 and 22.

> ##### WHAT IS KNOWN ALREADY: The molecular dialogue that takes place between the endometrial layer interacting with the more superficial embryonic cell layer is still poorly understood due to ethical and technical reasons, leaving the early stage of the human embryo implantation process mostly unexplored.

> ##### STUDY DESIGN, SIZE, DURATION: Paired samples of transcriptomic pre-receptive and receptive endometrial luminal epithelial cells from five healthy women were used. Five samples of the transcriptome of the trophoblast cells from 5-day-old implantation-competent blastocyst with normal chromosomal endowment and twenty-five samples of trophoblast cells from embryos with chromosomal abnormalities (trisomy on chromosomes 7, 11, 15, 21 or 22) were analysed.

> ##### MATERIALS, SETTING, METHODS: The analysis focused on the early stages of the implantation process. To this end, Agilent microarray analysis of luminal epithelial cells from human endometrial samples from five tissue-matched early secretory phase (LH+2, pre-receptive phase) and mid-secretory phase (LH+7, receptive phase) fertile females was performed. In addition, RNA-seq was carried out on trophoblast cells from 5-day-old euploid embryos and those with chromosomal abnormalities, namely trisomies on chromosome 7, 11, 15, 21 and 22 from a total of 30 samples, five for each condition. Both datasets were analysed separately and functional enrichment analyses of the differentially expressed genes (DEGs) were performed. In addition, Weighted correlation network analysis (WGCNA) algorithm was applied for enrichment analyses. Next, these results were analysed by generation of protein-protein (embryo-endometrium) interaction networks, using STRING/Cytoscape, and subsequent clustering of the generated networks accompanied by functional enrichment analysis of the sub-networks.

> ##### MAIN RESULTS AND THE ROLE OF CHANCE: We identified 52 DEGs with a key role in the regulation of defense response, ion homeostasis and regulation of the extracellular matrix processes involved in the endometrial maturation process into receptive phase within the luminal epithelial cells in healthy fertile women (comparisons of LH+2 vs. LH+7 samples). In addition, 113, 221, 215, 19, and 21 DEGs were identified after comparing the transcriptomics of trophoblast cells of euploid embryos and trophoblast cells in embryos with trisomy 7, 11, 15, 21, 22 respectively, showing notable differences between them. The constructed endometrium-embryo molecular dialogue was characterised, and the molecular differences in the crosstalk in case of chromosomal aberrations was identified.

> ##### LARGE SCALE DATA: Endometrial data are from the work of Evans et al. published in 2014 (doi: 10.1016/j.fertnstert.2014.04.005). Data from embryos have not yet been published.

> ##### LIMITATIONS, REASONS FOR CAUTION: The use of the microarray technique instead of newer techniques such as RNAseq or Single cell RNAseq in the endometrial samples, and the effect of the technical batch present in the embryo samples which has been adjusted as far as possible by bioinformatics techniques might not have enabled to detect all the possible molecular networks.

> #####  WIDER IMPLICATIONS OF THE FINDINGS: Our study provides novel insights into the so-far-understudied molecular interaction between the embryo and the endometrium in the first steps of the implantation process in humans. Knowledge of the DEGs and their possible functions in receptive luminal epithelial cells and chromosomally altered embryos compared to trophoblastic cells of euploid embryos, is fundamental to understanding human reproduction and the possible causes of implantation failure and infertility

<h2 align="center"> Dependencies </h2>
    
* [miARma-seq](https://github.com/eandresleon/miARma-seq)
* [R](https://www.r-project.org/) version 4.1.2 (2021-11-01)
* [RStudio-Posit](https://posit.co/) version 2021.09.2 Build 382
* R packages:
  * [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
  * [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
  * [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)
  * [genefilter](https://bioconductor.org/packages/release/bioc/html/genefilter.html)
  * [ggrepel](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html)
  * [ggfortify](https://cran.r-project.org/web/packages/ggfortify/index.html)
  * [cluster](https://cran.r-project.org/web/packages/cluster/index.html)
  * [factoextra](https://cran.r-project.org/web/packages/factoextra/index.html)
  * [ggplot2](https://ggplot2.tidyverse.org/)
  * [M3C](https://www.bioconductor.org/packages/release/bioc/html/M3C.html)
  * [AnnotationDbi](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)
  * [xlsx](https://cran.r-project.org/web/packages/xlsx/index.html)
  * [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
  * [png](https://cran.r-project.org/web/packages/png/index.html)
  * [curl](https://cran.r-project.org/web/packages/curl/vignettes/intro.html)
  * [corrplot](https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html)
  * [tidyr](https://tidyr.tidyverse.org/)
  * [dplyr](https://www.r-project.org/nosvn/pandoc/dplyr.html)
  * [WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/)
  * [ComBat](https://rdrr.io/bioc/sva/man/ComBat.html)
  
  
  
  
  
  


