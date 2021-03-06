#!/bin/bash

POPMAP=popmap/phylogen.popmap.final.txt
ALN=alignment_original/phylogen.ingroup.phy
INTER=aln_inter
FINAL=aln_final

mkdir -p $INTER
mkdir -p $FINAL

plist1=(25 50 75)
plist2=(25 50 75)
plist3=(0.0 0.01 0.03 0.05)
for i in ${plist1[@]}; do
  for j in ${plist2[@]}; do
    for m in ${plist3[@]}; do
	    # Filter by individuals ($i), populations ($j), and minor allele frequency ($m)
  nremover.pl -f $ALN -t phylip -i 0.$i -c 1.0 -o $INTER/missInd"$i"_Pop"$j"_maf"$m".phy -m -s -b -g 0.0 -a "$m";
  phylipFilterPops.pl -p $POPMAP -i $INTER/missInd"$i"_Pop"$j"_maf"$m".phy -n 0.$j -N 1.0 -o $FINAL/missInd"$i"_Pop"$j"_maf"$m".phy;
		done;
	done;
done

# No population-level filtering.
for m in ${plist3[@]}; do
	nremover.pl -f $ALN -t phylip -i 0.25 -c 1.0 -o $FINAL/missInd25_Pop100_maf"$m".phy -m -s -b -g 0.0 -a "$m";
	nremover.pl -f $ALN -t phylip -i 0.50 -c 1.0 -o $FINAL/missInd50_Pop100_maf"$m".phy -m -s -b -g 0.0 -a "$m";
	nremover.pl -f $ALN -t phylip -i 0.75 -c 1.0 -o $FINAL/missInd75_Pop100_maf"$m".phy -m -s -b -g 0.0 -a "$m";
	nremover.pl -f $ALN -t phylip -i 1.0 -c 1.0 -o $FINAL/missInd100_Pop100_maf"$m".phy -m -s -b -g 0.0 -a "$m";
done

# Convert all to STRUCTURE format.
for file in $FINAL/*.phy; do
  NEW="$(basename "$file" .phy).str"
  phylip2structure.pl -i "$file" -p $POPMAP -o $FINAL/$NEW;
  sort -n -k2 $FINAL/$NEW > $FINAL/"$(basename "$NEW" .str).sorted.str";
done

