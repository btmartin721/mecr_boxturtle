#!/bin/bash

# Per-individual missing data
i=50

# Per-population missing data
j=25

# Minor allele frequency filter
m=0.05

dataDIR=./data
dataFILE=missInd"$i"_Pop"$j"_maf"$m".onehot.txt
outputDIR=./vae_runs
newDIR=parsedoutput_rf_maf05/vae

awk '{print $1}' $dataDIR/$dataFILE > inds_tmp.txt

mkdir -p $newDIR

for r in {0..99}; do
	runFILE=r"$r".txt
	awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' OFS="\t" inds_tmp.txt "$outputDIR"/"$runFILE" > "$newDIR"/r"$r".txt;
done;

