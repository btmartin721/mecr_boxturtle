#!/bin/bash

plist=(0.0 0.01 0.03 0.05)
for i in {25..100..25}; do
	for j in {25..100..25}; do
		for m in ${plist[@]}; do
			dataDIR=./data
			dataFILE=missInd"$i"_Pop"$j"_maf"$m".onehot.txt
			outputDIR=missInd"$i"_Pop"$j"_maf"$m".pkl
			newDIR=missInd"$i"_Pop"$j"_maf"$m"_output
			awk '{print $1}' $dataDIR/$dataFILE > inds_tmp.txt
			mkdir -p $newDIR
			for r in {0..99}; do
				runFILE=r"$r".txt
				awk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' OFS="\t" inds_tmp.txt "$outputDIR"/"$runFILE" > "$newDIR"/r"$r".txt;
			done;
		done;
	done;
done

