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

See *example_files/popmap.txt*

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
3. Also splits the dataset into training/test sets to validate the model loss.  

### cMDS, isoMDS, t-SNE

To use the modified R script, supply a commands.txt file. See *example_files/commands_R.txt*  

The commands.txt file should have one line per replicate.  
You can supply dummy values if you are using multiple filtering parameters, and replace the dummy variables using sed when you run the analysis.  

Run ```Derkarabetian_etal_2019/PCA-DAPC-RF-tSNE_gridSearch_maf.r -h``` to see a list of options you can change in the commands.txt file.

Then, to run the replicates in parallel, use [GNU parallel](https://www.gnu.org/software/parallel/man.html) like so:  

```cat commands.txt | parallel```

I found that sometimes it would use more memory than I had available, so if that's the case you can set some limitations:  

```cat commands.txt | parallel -j 8```

This limits the number of parallel jobs to 8 at a time.  
You can also try setting a memory usage limit, which won't start new jobs if available memory is less than what you specify.  E.g.,  

```cat commands.txt | parallel --memfree 3.0G --delay 30```

Here, if there is less than 3 GB available, parallel won't start a new job until more memory frees up.  
The delay waits for 30 seconds before trying to start a new job.  

Finally, if some random replicates keep failing, you can try the retries option:  

```cat commands.txt | parallel --retries 10```

This will retry the replicate up to 10 times if it fails.  

### VAE

To run VAE with multiple replicates, just run the modified python script.  Use the -h option to see the command-line options.  

```Derkerabetian_etal_2019/sp_deli_clust_commandline_noClust.py -h```

To run the replicates in parallel, use a commands.txt file with one replicate per line, as discussed above with the R script. 
Then pass it to GNU parallel.   

## Preparing UML output 

### cMDS, isoMDS, t-SNE  

1. Run *parsing_output/bash/copy_uml_output.sh*
2. Run *parsing_output/bash/clean_uml_output.sh*

The output for all the cMDS, isoMDS, and t-SNE replicates across all algorithms will be in a single directory.  
You can run the *parsing_output/bash/copy_uml_output.sh* script to neatly organize the output files.
This will move with the replicates for each algorithm into their own directories.  

Once all the output is moved, use the *parsing_output/bash/clean_uml_output.sh* script to fix some inconsistencies in file formats and add switch the 

### VAE

3. Run *parsing_output/python/vae_dbscan.py* to do an automated clustering of the VAE latent variables.  

It uses DBSCAN to do the clustering.
The standard deviation (SD) latent variable is averaged across all replicates, and DBSCAN clusters with an epsilon (eps) of 2 X SD.  

4. Run *parsing_output/bash/add_sampleids.sh*

Use the add_sampleids.sh script to switch the first column in each replicate with sampleIDs.  
This makes for easier plotting later.  

## Parsing UML Output

5. Run the *parsing_output/aligningK/UML_alignK_missDataRuns_maf.R* script to align K across replicates. 

First, you need to align K across replicates. We implemented the CLUMPAK algorithm (Kopelman et al. 2015) using the PopHelper R package (Francis 2017).

To use this script, replicates for a single algorithm should be in its own directory, **with nothing else in it**.  E.g. if you did 100 cMDS gapstat replicates, there should be 100 files in the *cmds/gapstat* directory.  

The script loads the replicates into a single data.frame, aligns K for each clustering algorithm, 
and saves RDS objects for each algorithm so they can be loaded in downstream R scripts. 

The .rds objects are lists of data.frames, with the first data.frame being the K-aligned replicates as columns and the second beign assignment proportions.  

Francis, R. M. (2017). pophelper: an R package and web app to analyse and visualize population structure. Molecular ecology resources, 17(1), 27-32.

Kopelman, N. M., Mayzel, J., Jakobsson, M., Rosenberg, N. A., & Mayrose, I. (2015). Clumpak: a program for identifying clustering modes and packaging population structure inferences across K. Molecular Ecology Resources, 15(5), 1179-1191.


## Plotting UML Output

6. Run *plotting/plot_missData_stats.R* to plot summary statistics.  

Once you have the .rds objects, you can make summary statistic plots. The script makes heatmaps, Tukey-style boxplots, and regression plots.  
You will have to tweak the code in this one. It isn't full automated. 

7. Run *plotting/plotUML_missData_maf.R* to make structure-style barplots

Here, you can make structure-style barplots with assignment proportions among the UML replicates.
It saves one barplot per clustering algorithm and puts them all in one PDF file. 
If you ran multiple filtering parameters you can run them all in the for loop.  

8. Run *plotting/plotDatedTree.R* to add a phylogeny to the structure-style barplots. 

This uses ggtree to add a phylogeny to a page of structure-style barplots. The individuals from the barplots and phylogeny will be aligned.
It makes for easier interpretation. 
You'll have to tweak the code, as we don't have it automated. Currently, it is tailored to our specific dataset and aesthetic choices.  

It's also set up to read a nexus file output from LSD2 divergence dating, as implemented in IQ-TREE.   

Well, that's all our scripts :)
