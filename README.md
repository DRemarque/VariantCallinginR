# VariantCallinginR

### Abstract
Linkage and association studies have gained popularity with the rise of fast and cheap whole genome sequencing. These studies rely on the use of genetic markers - short variable DNA sequences with a known location in the genome - to identify haplotypes and phenotypes related to specific variant alleles. The discovery and analysis of genetic markers within a genome, also known as variant calling, currently requires the use of multiple specialised software controlled via terminal commands. Therefore, the process of variant calling is laborious for research groups without access to a bio-informatician. The goal of this research project was to lower the threshold by creating an easy-to-use variant calling pipeline in R, with which most scientists already have ample experience. The final product is the R package VariantCallinginR, where each step of variant calling is compiled into an R function. Now, the entire variant discovery pipeline from raw data to variant call file can be executed in ten lines of code or less, including regular quality controls. In order to test the functionality of the package, compound variant discovery has been performed for nineteen previously sequenced, high interest Cannabis sativa strains aligned to the CBDRx reference strain. In total, over 11.5 million single nucleotide polymorphisms and 0.8 million microsatellites were discovered within the Cannabis genome.

### Package Installation
In R, use the following command to install the lastest version of VariantCallinginR from the repository:
```{}
install.packages("devtools")
library(devtools)
install_github('VariantCallinginR','DRemarque')
```
