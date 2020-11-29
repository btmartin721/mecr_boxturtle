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

Derkarabetian S, Castillo S, Koo, PK, Ovchinnikov S, Hedin M. 2019. A demonstration of unsupervised machine learning in species delimitation. Molecular Phylogenetics and Evolution. 139: 106562. https://doi.org/10.1016/j.ympev.2019.106562  

Our primary modifications to the scripts from Derkarabetian et al. (2019) are to 1) allow parallel execution of multiple replicates, and 2) perform a t-SNE perplexity grid search.  

Other minor modifications in the R script:
1. Saves tab-separated files for each replicate.  
  + Two-column files resemble popmaps, except they have the clustering assignments for each individual instead of population IDs.  
  + This allows the replicates to be parsed later.  

Other minor modifications in the VAE script:
1. Early stopping callback to reduce overfitting (restores the best model prior to the tolerance period).
2. Saves latent variables from all replicates as pickle objects  
  + This allows the unsupervised clustering algorithm DBSCAN to be applied with our vae_dbscan.py script.  
  
To use the modified R script, supply a commands.txt file. It looks like this:

```

```
  

