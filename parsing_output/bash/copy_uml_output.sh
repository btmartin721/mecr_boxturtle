#!/bin/bash


# Make separate output directory for each analysis.
for i in {25..100..25}; do
	        for j in {25..100..25}; do
			for m in `seq 0.01 0.02 0.05`; do
				mkdir -p missInd"$i"_Pop"$j"_maf"$m"/{cmds/{gapstat,hierarchical,pam/{cmds_groups,prox}},isomds/{gapstat,hierarchical,pam},tsne/{gapstat/{p5,p10,p15,p20,p25,p30,p35,p40,p45,p50},hierarchical/{p5,p10,p15,p20,p25,p30,p35,p40,p45,p50},pam/{p5,p10,p15,p20,p25,p30,p35,p40,p45,p50}}};
		done;
	done;
done

# Define destination directories
cmdsGS=cmds/gapstat
cmdsHIER=cmds/hierarchical
cmdsPAMcmds=cmds/pam/cmds_groups
cmdsPAMprox=cmds/pam/prox

isomdsGS=isomds/gapstat
isomdsHIER=isomds/hierarchical
isomdsPAM=isomds/pam

tsneGS=tsne/gapstat
tsneHIER=tsne/hierarchical
tsnePAM=tsne/pam

# Copy results to their own independent directories for each individual analysis.
for i in {25..100..25}; do 
	for j in {25..100..25}; do
		for m in `seq 0.01 0.02 0.05`; do
	        	sourceDIR=../missInd"$i"_Pop"$j"_maf"$m"
       			destDIR=./missInd"$i"_Pop"$j"_maf"$m"	       
			cp $sourceDIR/uml_*_cmds_pam_gapstatisticK_run*.txt $destDIR/$cmdsGS/;
			cp $sourceDIR/uml_*_cmds_RFK_hierarchicalClustering_grps_run*.txt $destDIR/$cmdsHIER/
			cp $sourceDIR/uml_*_cmds_grps_pam_K*_run*.txt $destDIR/$cmdsPAMcmds/
			cp $sourceDIR/uml_*_cmds_pam_proximityScores_K*_run*.txt $destDIR/$cmdsPAMprox/

			cp $sourceDIR/uml_*_isomds_grps_pam_gapstatisticK_run*.txt $destDIR/$isomdsGS/
			cp $sourceDIR/uml_*_isomds_grps_hierarchicalK_run*.txt $destDIR/$isomdsHIER/
			cp $sourceDIR/uml_*_isomds_grps_pam_K*_run*.txt $destDIR/$isomdsPAM/
		
			# For each perplexity
			for p in {5..50..5}; do
				# All missing values in separate directories.
				cp $sourceDIR/uml_*_tsnep"$p"_grps_gapstatisticK_pam_run*.txt $destDIR/$tsneGS/p"$p"
				cp $sourceDIR/uml_*_tsnep"$p"_grps_hierarchcialK_run*.txt $destDIR/$tsneHIER/p"$p"
				cp $sourceDIR/uml_*_tsnep"$p"_grps_pam_K*_run*.txt $destDIR/$tsnePAM/p"$p";

			done;
		done; 
	done;
done


