# VariantCallinginR

### Abstract
Linkage and association studies have gained popularity with the rise of fast and cheap whole genome sequencing. These studies rely on the use of genetic markers - short variable DNA sequences with a known location in the genome - to identify haplotypes and phenotypes related to specific variant alleles. The discovery and analysis of genetic markers within a genome, also known as variant calling, currently requires the use of multiple specialised software controlled via terminal commands. Therefore, the process of variant calling is laborious for research groups without access to a bio-informatician. The goal of this research project was to lower the threshold by creating an easy-to-use variant calling pipeline in R, with which most scientists already have ample experience. The final product is the R package VariantCallinginR, where each step of variant calling is compiled into an R function. Now, the entire variant discovery pipeline from raw data to variant call file can be executed in ten lines of code or less, including regular quality controls. In order to test the functionality of the package, compound variant discovery has been performed for nineteen previously sequenced, high interest Cannabis sativa strains aligned to the CBDRx reference strain. In total, over 11.5 million single nucleotide polymorphisms and 0.8 million microsatellites were discovered within the Cannabis genome.

### Overview of variant discovery and corresponding R functions:
![Image of Method](https://octodex.github.com/images/yaktocat.png)

### Package Installation
In R, use the following command to install the lastest version of VariantCallinginR from the repository:
```{r}
# Install access to github
install.packages("devtools")

# Download VariantCallinginR
devtools::install_github('DRemarque/VariantCallinginR',subdir='VariantCallinginR')
```
### Installation additional software  
To access the GATK and BWA-MEM functions of the package, additional software is required.

For Windows systems, this Unix based software requires a virtual Linux environment provided by Docker Toolbox:
1. Download the [Docker Toolbox installation software][ref-1].  
   For Windows, select the latest .exe file and for Mac, select the latest .pkg file.
2. Double click the downloaded file and follow the installation prompts
3. To verify the installation, run the Docker Quickstart terminal once. This will setup the virtual environment.
4. For Windows, check that the Docker Quickstart terminal installed Git Bash (bash.exe) at the default location: C:/Program Files/Git/bin/bash.exe  
   For Mac, check that the terminal.exe is at its default location: /Applications/Utilities/Terminal.app.
Note that Docker Toolbox is not available for Linux systems as these do not require a virtual Linux environment. 

For Linux systems, the base software itself has to be installed. The latest installation software versions and instructions can be found [here][ref-2] for GATK, BWA and Samtools and Picard (to support GATK).

### Citation Information
_To follow soon_

[ref-1]: https://github.com/docker/toolbox/releases
[ref-2]: https://gatk.broadinstitute.org/hc/en-us/articles/360041320571--How-to-Install-all-software-packages-required-to-follow-the-GATK-Best-Practices
