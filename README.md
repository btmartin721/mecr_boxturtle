# mecr_boxturtle
Collection of scripts for doing and plotting machine learning species delimitation results

## Filtering

To filter the dataset various ways, use the filter_maf.sh script. It requires two external Perl scripts:  

[nremover.pl](https://github.com/tkchafin/scripts)  
[phylipFilterPops.pl](https://github.com/tkchafin/scripts)  

You will also need a PHYLIP-formatted alignment and a population map file.  
The population map (popmap) is a two-column, tab-separated file with individualIDs as the first column and populationIDs as the second. It is formatted like this:

```
ind1  pop1
ind2  pop1
ind3  pop2
ind4  pop2
```  

If you want to change the filtering parameters, you'll need to change the values in the filter_maf.sh script. Just change what's in the ```plist``` variables.  

## Running UML Analyses

To run the UML analyses, we used modified versions of the R and python scripts from [shahanderkarabetian/uml_species_delim](https://github.com/shahanderkarabetian/uml_species_delim) and [sokrypton/sp_deli](https://github.com/sokrypton/sp_deli). If you use them, please cite the original authors: 

