#!/bin/bash

# Individual missing data
i=50

# Population missing data
j=25

# Minor allele frequency filter
m=0.05

# t-SNE perplexity
p=15

PREFIX=gila_maf05

myDIR=parsedoutput_rf_maf05
awk '{print $1}' $myDIR/cmds/gapstat/uml_"$PREFIX"_cmds_pam_gapstatisticK_run1.txt > inds_tmp.txt
# Fix format of tsne hierarchical. It has quotes and a column header, and it shouldn't.
for k in {1..100}; do
	HIER=tsne/hierarchical/p"$p"
	GS=tsne/gapstat/p"$p"
	PAM=tsne/pam/p"$p"
	tsneHIER=$myDIR/$HIER/uml_"$PREFIX"_tsnep"$p"_grps_hierarchcialK_run"$k".txt
	tsneHIER2=$myDIR/$HIER/uml_"$PREFIX"_tsnep"$p"_grps_hierarchcialK_run"$k"_mod.txt
	tsneGS=$myDIR/$GS/uml_"$PREFIX"_tsnep"$p"_grps_gapstatisticK_pam_run"$k".txt
	tsneGS2=$myDIR/$GS/uml_"$PREFIX"_tsnep"$p"_grps_gapstatisticK_pam_run"$k"_mod.txt
	tsnePAM=$myDIR/$PAM/uml_"$PREFIX"_tsnep"$p"_grps_pam_K*_run"$k".txt
	tsnePAM2=$myDIR/$PAM/uml_"$PREFIX"_tsnep"$p"_grps_pam_K*_run"$k"_mod.txt

	# Fix incorrect format of tsne hierarchical.
	sed 's/"//g' $tsneHIER | grep -v "^x$" | tr ' ' '\t' > $tsneHIER2
	mv $tsneHIER2 $tsneHIER
	awk 'FNR==NR{a[NR]=$2;next}{$2=a[FNR]}1' OFS="\t" $tsneHIER inds_tmp.txt > $tsneHIER2
	# Replace first column with actual individual IDs.
	mv $tsneHIER2 $tsneHIER                      
			
	# Replace first column with actual individual IDs.
        awk 'FNR==NR{a[NR]=$2;next}{$2=a[FNR]}1' OFS="\t" $tsneGS inds_tmp.txt > $tsneGS2
        mv $tsneGS2 $tsneGS;

	# Replace first column with actual individual IDs.
	awk 'FNR==NR{a[NR]=$2;next}{$2=a[FNR]}1' OFS="\t" $tsnePAM inds_tmp.txt > $tsnePAM2
	mv $tsnePAM2 $tsnePAM;
done;

			
for r in {1..100}; do
	cmdsHIER=$myDIR/cmds/hierarchical/uml_"$PREFIX"_cmds_RFK_hierarchicalClustering_grps_run"$r".txt
	cmdsHIER2=$myDIR/cmds/hierarchical/uml_"$PREFIX"_cmds_RFK_hierarchicalClustering_grps_run"$r"_mod.txt
	
	isomdsHIER=$myDIR/isomds/hierarchical/uml_"$PREFIX"_isomds_grps_hierarchicalK_run"$r".txt
	isomdsHIER2=$myDIR/isomds/hierarchical/uml_"$PREFIX"_isomds_grps_hierarchicalK_run"$r"_mod.txt
	
	awk 'FNR==NR{a[NR]=$2;next}{$2=a[FNR]}1' OFS="\t" $cmdsHIER inds_tmp.txt > $cmdsHIER2
	mv $cmdsHIER2 $cmdsHIER

	awk 'FNR==NR{a[NR]=$2;next}{$2=a[FNR]}1' OFS="\t" $isomdsHIER inds_tmp.txt > $isomdsHIER2
	mv $isomdsHIER2 $isomdsHIER;
done;

