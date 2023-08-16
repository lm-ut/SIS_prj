#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import math
import sys


REF_name = sys.argv[1]
Target_name = sys.argv[2]


#Loading files: frq file for the reference group & tped and frq file for Target group

FRQ_file = pd.read_csv(REF_name+".frq", delim_whitespace=True)
tped_file = pd.read_csv(Target_name+".tped",delim_whitespace=True, header=None)
FRQ_Target_file = pd.read_csv(Target_name+".frq", delim_whitespace=True)


#Check that the Target and Reference files have the same number of SNPs

if len(FRQ_file) == len(tped_file):
    print("Same length, go")
else:
    print("Different length, check")


#Calculate AVG

def AVG(list_het):
    mysum = 0
    count = 0
    for value in list_het:
        mysum += float(value)
        count += 1
    
    print(mysum/count)
    AVG_est = mysum/count
    return(AVG_est)


# Calculate 2pq
def HET_2pq(insert_float):
    h2pq = 2*insert_float*(1-insert_float) 
    return h2pq

# This is the core estimation: 
# where we extract the number of alleles available per each SNPs in the Reference group and then select that number of alleles in the Target group to then compute the frequencies.

list_h2pq = []

for locus in range(0,len(FRQ_file)):
    NCHR_found = FRQ_file.at[locus,"NCHROBS"]
    alleles_tped = tped_file.iloc[locus,4:len(tped_file.columns)]
    
    # Selecting one NON-ZERO allele each pair 
    p_allele = alleles_tped.iloc[::2]
    p_allele_NoZeros = p_allele[p_allele > 0]
    
    # Selecting only CHRs available in REF file
    p_Target_BasedOn_Reference = p_allele_NoZeros[0:NCHR_found]

    #Summary on dataset 
    nZeros = alleles_tped.isin([0]).sum()
    NoZeroValues = len(p_allele_NoZeros)
    NumOf_p = len(p_allele)
    NumOf_Chr = len(alleles_tped)
    NumOf_p_BasedOnRef = len(p_Target_BasedOn_Reference)
    #print(NumOf_p_BasedOnRef)
    
    # Selecting Series that are not empty to then calculate the het
    if (len(p_Target_BasedOn_Reference)) != 0:    
        p_het = p_Target_BasedOn_Reference.value_counts(normalize=True).tolist()[0]
        h2pq = HET_2pq(p_het)
        list_h2pq.append(h2pq)


# Calculating AVG of Target based on Reference
AVG_est = AVG(list_h2pq)


## To check if any mistakes are made,
# I am estimating the overall het in the tped of the target 
# Considering all availables sites, not only those present in REF as well 
# I will then estimate het on the Target frq file as well and compare with the results obtained in this box

TAR_h2pq = []

for locus in range(0,len(FRQ_Target_file)):
    alleles_tped = tped_file.iloc[locus,4:len(tped_file.columns)]
    
    # Selecting one NON-ZERO allele each pair 
    p_allele = alleles_tped.iloc[::2]
    p_allele_NoZeros = p_allele[p_allele > 0]
    
    if (len(p_allele_NoZeros)) != 0:    
        p_het = p_allele_NoZeros.value_counts(normalize=True).tolist()[0]
        h2pq = HET_2pq(p_het)
        TAR_h2pq.append(h2pq)

# Calculating AVG of Target based on Reference
TAR_AVG_est = AVG(TAR_h2pq)


## HET from frq file
MAF_Cleaned = FRQ_Target_file["MAF"].dropna()

frq_h2pq = []

for maf in MAF_Cleaned:
    h2pq = HET_2pq(maf)
    frq_h2pq.append(h2pq)


# Calculating AVG of Target based on Reference
frq_AVG_est = AVG(frq_h2pq)


with open(Target_name + "_BasedOn" + REF_name + "_HET_AVG.txt", "w") as out_file:
    out_file.write("Heterozigosity of Target in .frq " + f"{frq_AVG_est}\n")
    out_file.write("Heterozigosity of Target by my estimation " + f"{TAR_AVG_est}\n")
    out_file.write("Heterozigosity based on REF " + f"{AVG_est}\n")




