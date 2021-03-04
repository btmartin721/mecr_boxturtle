#!/bin/bash


# Make separate output directory for each analysis.
mkdir -p parsedoutput_rf_maf05/{cmds/{gapstat,hierarchical,pam/{cmds_groups,prox}},isomds/{gapstat,hierarchical,pam},tsne/{gapstat/p15,hierarchical/p15,pam/p15}};

# Individual missing data
i=50

# Population missing data
j=25

# Minor allele frequency filter
m=0.05

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
sourceDIR=./output_maf05
destDIR=./parsedoutput_rf_maf05	       
cp $sourceDIR/uml_*_cmds_pam_gapstatisticK_run*.txt $destDIR/$cmdsGS/;
cp $sourceDIR/uml_*_cmds_RFK_hierarchicalClustering_grps_run*.txt $destDIR/$cmdsHIER/
cp $sourceDIR/uml_*_cmds_grps_pam_K*_run*.txt $destDIR/$cmdsPAMcmds/
cp $sourceDIR/uml_*_cmds_pam_proximityScores_K*_run*.txt $destDIR/$cmdsPAMprox/
cp $sourceDIR/uml_*_isomds_grps_pam_gapstatisticK_run*.txt $destDIR/$isomdsGS/
cp $sourceDIR/uml_*_isomds_grps_hierarchicalK_run*.txt $destDIR/$isomdsHIER/
cp $sourceDIR/uml_*_isomds_grps_pam_K*_run*.txt $destDIR/$isomdsPAM/

# t-SNE perplexity
p=15

# All missing values in separate directories.
cp $sourceDIR/uml_*_tsnep"$p"_grps_gapstatisticK_pam_run*.txt $destDIR/$tsneGS/p"$p"
cp $sourceDIR/uml_*_tsnep"$p"_grps_hierarchcialK_run*.txt $destDIR/$tsneHIER/p"$p"
cp $sourceDIR/uml_*_tsnep"$p"_grps_pam_K*_run*.txt $destDIR/$tsnePAM/p"$p";


