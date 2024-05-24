# SIS_prj

### This script was used in [this study](https://onlinelibrary.wiley.com/doi/10.1002/ajpa.24911) 

The script will estimate allele frequencies on population Y based on the alleles available in population X, it was origianlly used for [this work](https://onlinelibrary.wiley.com/doi/10.1002/ajpa.24911).

### Usage: python script.py popX_file popY_file  

Where popX_file is the prefix of .frq plink file obtained from the reference population.  
While popY_file is the prefix of .frq and .tped plink files obtained from the population Y on which we wanted to estimate the allele frequencies, limited on the alleles shared with population X.
