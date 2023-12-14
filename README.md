# O5C2-treated BA.5 COVID-19 mouse lung project - DGE analysis pipeline

This repository includes code for the analysis pipeline of O5C2-treated BA.5 COVID-19 mouse lung project, from raw bulk RNA-seq data to differential gene expression output. Codes used in this pipeline includes both shell and R scripts, and the program/tools used are as follow:

Shell:
- [ ] [fastp](https://github.com/OpenGene/fastp)
- [ ] [STAR](https://github.com/alexdobin/STAR)
- [ ] [samtools](https://github.com/samtools/samtools)
- [ ] [qualimap](https://github.com/EagleGenomics-cookbooks/QualiMap)
- [ ] [featureCounts](https://subread.sourceforge.net/)

R:
- [ ] [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [ ] [ggplot2](http://cran.r-project.org/web/packages/ggplot2/index.html)
- [ ] [EnhancedVolcano](https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html)
- [ ] [goseq](https://bioconductor.org/packages/release/bioc/html/goseq.html)
- [ ] [patchwork](http://cran.r-project.org/web/packages/patchwork/index.html)

There are main 2 scripts in this repo:
1. **fastp_STAR_qualimap.sh** - shell script to (1) quality check using *fastp*, (2) mapping to reference genome using *STAR*, (3) check for mapping metrics using *qualimap*, and lastly (4) count the read count of each gene using *featureCounts*.
2. **DESeq2_plot.R** - R scruipt to conduct differential gene expression (DGE) analysis using *DESeq2*, and plot the main figures related to DGE analysis. 
