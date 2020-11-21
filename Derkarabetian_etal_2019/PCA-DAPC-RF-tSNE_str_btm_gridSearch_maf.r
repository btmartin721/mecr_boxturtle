#########################################################################################
#########################################################################################
## Script from:
## Derkarabetian S., Castillo S., Peter K.K., Ovchinnikov S., Hedin M. 
## 2019. 
## "An Empirical Demonstration of Unsupervised Machine Learning in Species Delimitation".
## Molecular Phylogenetics and Evolution, 139: 106562. 
## https://doi.org/10.1016/j.ympev.2019.106562.
#########################################################################################
#########################################################################################

##########################################################
## Script has been modified from the original 
## to run multiple datasets with varying missing 
## data and minor allele frequency filters.
## Also does a grid search for t-SNE perplexity settings.
## Finally, it is modified to perform multiple replicate runs per
## clustering algorithm.
##########################################################

if (!require("optparse")) stop("Error: The required package optparse is not installed")
library("optparse")

# Set command-line arguments
option_list <- list(make_option(c("-s", "--str"), 
                                type="character", 
                                default=NULL, 
                                help="Required; Input structure filename; default = NULL", 
                                metavar="character"),
                    make_option(c("-r", "--run"),
                                type="integer",
                                default=NULL,
                                help="Required; Run number",
                                metavar="integer"),
                    make_option(c("-p", "--prefix"),
                                type="character",
                                default=NULL,
                                help="Required; Specify prefix for output files",
                                metavar="character"),
                    make_option(c("--nsites"),
                                type="integer",
                                default=NULL,
                                help="Required; Number of sites in structure input file",
                                metavar="integer"),
                    make_option(c("--ninds"),
                                type="integer",
                                default=NULL,
                                metavar="integer",
                                help="Required; Number of individuals in structure input file")
                )
                   

opt_parser <- OptionParser(option_list=option_list,
                           description="Rscript to run uml_species_delim")

opt <- parse_args(opt_parser)

#required packages
library("adegenet")
library("randomForest")
library("PCDimension", lib.loc="/home/btm002/anaconda3/lib/R/library")
library("mclust")
library("cluster")
library("MASS")
library("factoextra")
library("tsne")

################################################################################

## Returns list of pam clustering object as element 1 and optimal K integer as element 2
pamK <- function(rftest, numgroups) {
  means <- vector()
  for (i in seq(2, numgroups)) {
    clust_tmp <- pam(rftest, i)
    means[i - 1] <- mean(silhouette(clust_tmp)[, "sil_width"])
  }
  
  pam_clust_bestK <- which(means == max(means)) + 1
  pam_clust_bestmean <- max(means)
  pam <- pam(rftest, pam_clust_bestK)
  
  return(list(pam, pam_clust_bestK))
  
}

required.args <- function(arg, argString) {
  
  print(arg)
  if(is.null(arg)) {
    
    stop(paste0("\n\nError: ", argString, " is a required argument\n"))
  }
  
}

read.infile <- function(infile, header) {
  
  if(!is.null(infile)) {
    
    if (!file.exists(infile)) {
      stop("Error: Could not find required input file; aborting program")
    }
  }
}

MyKmeansFUN <- function(x,k) list(cluster=kmeans(x, k, iter.max=20, nstart = 25))

# Make sure required arguments are supplied at command-line
required.args(opt$prefix, "--prefix")
required.args(opt$str, "--str")
required.args(opt$run, "--run")
required.args(opt$nsites, "--nsites")
required.args(opt$ninds, "--ninds")


## This assumes your input file is structure/adegenet format (.str). You can also import data in .csv format, although
## depending on data type, you will have to convert:
## If importing .csv with raw nucleotides/haplotypes, convert to factor, go to random forest step.
## If importing .csv with nucleotide data in one-hot format, covert to numeric, go to PCA step.

###############################################

# Used 7395 sites and 37 inds
data <-
  read.structure(
    opt$str,
    n.ind = opt$ninds,
    n.loc = opt$nsites,
    onerowperind = FALSE,
    col.lab = 1,
    col.pop = 2,
    col.others = NULL,
    row.marknames = NULL,
    NA.char = "-9",
    pop = NULL,
    ask = FALSE,
    quiet = FALSE
  )

data_scaled <- scaleGen(data, center=FALSE, scale=FALSE, NA.method=c("zero"), nf)

pdf(file = paste0(opt$prefix, "pca-dapc_plots_run", opt$run, ".pdf"), onefile = T)

# PCA, can adjust nf to include more components
pca1 <- dudi.pca(data_scaled, center=TRUE, scale=TRUE, scannf=FALSE, nf=2)
plot(pca1$li, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="PCA", pch=16)
s.label(pca1$li, clabel=0.5, grid=0)

# DAPC (interactive, requires input)
# max.n.clust equal to number of pops
# Chose 37 PCs, 4 clusters.
clusters <- find.clusters(data, max.n.clust=9, n.iter=1e6, n.start=1000, pca.select = "percVar", choose.n.clust = FALSE, perc.pca=100, stat="BIC", criterion="min")

results <- dapc(data, pop = clusters$grp, pca.select = "percVar", perc.pca=100, n.da=nlevels(clusters))

# 3D PCA and DAPC
pca2 <- dudi.pca(data_scaled, center=TRUE, scale=TRUE, scannf=FALSE, nf=3)

df12 <- data.frame(pca2$li$Axis1, pca2$li$Axis2)
df23 <- data.frame(pca2$li$Axis2, pca2$li$Axis3)

# Axes 1 and 2
s.class(
  pca2$li,
  fac = results$grp,
  xax = 1,
  yax = 2,
  col = c("#CCCC00", "darkorange", "#4B0082", "#228b22"),
  label = c(
    "bauri",
    "carolina/coahuila/major",
    "ornata/luteola",
    "triunguis/mexicana"
  )
)

# Axes 2 and 3
s.class(
  pca2$li,
  fac = results$grp,
  xax = 2,
  yax = 3,
  col = c("#CCCC00", "darkorange", "#4B0082", "#228b22"),
  label = c(
    "bauri",
    "carolina/coahuila/major",
    "ornata/luteola",
    "triunguis/mexicana"
  )
)

write.table(paste0(opt$prefix, "_dapc_groups_run", opt$run, ".txt"), x = results$grp, quote = F, sep = "\t", col.names = F)

# Do cross-validation for DAPC
x <- data
mat <-
  as.matrix(tab(x, NA.method="mean"))

grp <- pop(x)

writeLines("\nDoing DAPC cross-validation...\n")

xval <-
  xvalDapc(
    mat,
    clusters$grp,
    n.pca.max = results$n.pca,
    result = "groupMean",
    center = T,
    scale = F,
    n.pca = NULL,
    n.rep = 1000,
    xval.plot = T
  )

saveRDS(xval, file = paste0("xval_run", opt$run, ".rds"))

writeLines("\nDone!\n")

results <-
  dapc(
    data,
    pop = clusters$grp,
    n.pca = xval$DAPC$n.pca,
    n.da = xval$DAPC$n.da
  )

writeLines(paste0("\nUsing ", xval$DAPC$n.pca, " PCA axes and ", xval$DAPC$n.da, " discriminant functions\n"))


assignplot(results)
compoplot(results)
grp_k <- nlevels(clusters$grp)

write.table(x = results$grp, file = paste0(opt$prefix, "_dapc_pops_run", opt$run, ".txt"), quote = F, sep = "\t", col.names = F)
print(paste0("DAPC groups: ", results$grp))

# The next two lines are to color my groups like I want.

# PCA with DAPC groups
plot(
  pca1$li,
  xlab = "Scaling Dimension 1",
  ylab = "Scaling Dimension 2",
  main = "PCA with DAPC clusters",
  col = results$grp,
  pch = 16
)

s.label(pca1$li, clabel = 0.5, grid = 0)

dev.off()

###############################################
###############################################
# into the Random Forest, unsupervised
###############################################
###############################################

rm(xval)
rm(x)
rm(pca2)
rm(clusters)
gc()

# convert genind scaled data to factors for randomForest
data_conv <- as.data.frame(data_scaled)
data_conv[is.na(data_conv)] <- ""
data_conv[sapply(data_conv, is.integer)] <- lapply(data_conv[sapply(data_conv, is.integer)], as.factor)
data_conv[sapply(data_conv, is.character)] <- lapply(data_conv[sapply(data_conv, is.character)], as.factor)
nsamp <- nrow(data_conv)

writeLines("\nDoing Random Forest...\n")

# unsupervised random forest
rftest <- randomForest(data_conv, ntree = 1e4)

saveRDS(rftest, file = paste0("rf_run", opt$run, ".rds"))

writeLines("\nDone!\n")

###############
# classic MDS
###############

writeLines("\nDoing cMDS...")

# cMDS with optimal number of components to retain using broken-stick
# may need to adjust number of dimensions if given error
pdf(file = paste0(opt$prefix, "cmdsPlot_run", opt$run, ".pdf"), onefile = T, width = 40, height = 40)
  cmdsplot1 <- NULL
  while (is.null(cmdsplot1) && nsamp > 0) {
    print("nsamp")
    nsamp <- nsamp - 1
    tryCatch({
      cmdsplot1 <- MDSplot(rftest, results$grp, nsamp)
    },
    error = function(e){
      print("Too many dimensions on MDSplot. Retrying with nsamp - 1"); NULL;},
  warning = function(w){
    print("Too many dimensions on MDSplot. Retrying with nsamp - 1"); NULL;}
  )}
dev.off()

pdf(file = paste0(opt$prefix, "RF_plots_run", opt$run, ".pdf"), onefile = T)
cmdsplot_bstick <- bsDimension(cmdsplot1$eig)
cmdsplot2 <- MDSplot(rftest, results$grp, cmdsplot_bstick)

# cMDS with optimal DAPC k and clusters
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS DAPC optimal K and clusters", col=results$grp, pch=16)

# pam clustering on proximity scores with optimal k from DAPC
DAPC_pam_clust_prox <- pam(rftest$proximity, grp_k)

writeLines("\nDoing PAM with Proximity Scores\n")

# Returns list with element 1 being pam clustering object with best K
# Element 2 is the optimal K as an integer  
cmdsproxK <- pamK(rftest$proximity, 9)

# Write mean silhouette of optimal K to text file  
write.table(x = mean(silhouette(cmdsproxK[[1]])[, "sil_width"]), file = paste0(opt$prefix, "_dapcK_prox_pam_mean_meansilhouette", cmdsproxK[[2]], "_run", opt$run, ".txt"))

# cMDS with optimal k of DAPC and clusters via PAM
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main=paste0("cMDS K=", cmdsproxK[[2]], " and clusters (PAM clustering)"), col=cmdsproxK[[1]]$clustering, pch=16)

s.label(cmdsplot2$points, clabel=0.5, grid=0)

writeLines("Done!\n")

writeLines("Doing PAM with cMDS...\n")

# pam clustering on cMDS output with optimal k from PAM clustering
cmdsK <- pamK(cmdsplot1$points, 9)

# mean silhouette width of k
# can test multiple k values and select optimal k with highest value
write.table(
  x = mean(silhouette(cmdsK[[1]])[, "sil_width"]),
  file = paste0(opt$prefix, "_dapcK_cmds_pam_mean_meansilhouette", cmdsK[[2]], "_run", opt$run, ".txt"),
  col.names = F,
  quote = F,
  row.names = F
)

# Plot PAM optimal K
plot(cmdsplot1$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main=paste0("cMDS K=", cmdsK[[2]], " and clusters (PAM clustering)"), col=cmdsK[[1]]$clustering, pch=16)

writeLines("Done!\n")

writeLines("\nDoing cMDS with Gap Statistic...\n")
# determine optimal k from cMDS using gap statistic with PAM clusters from proximity scores
# can adjust ust.max
cmds_nbclust <- fviz_nbclust(cmdsplot1$points, kmeans, nstart=25, iter.max=50, method = "gap_stat", k.max = 9, nboot = 1000) + labs(subtitle = "Gap Statistic Method")
cmds_nbclust
cmds_nbclust_k <- cmds_nbclust[["layers"]][[4]][["data"]][["xintercept"]]
# pam clustering with optimal k from gap statistic
cmds_nbclust_clust <- pam(cmdsplot1$points, cmds_nbclust_k)

# cMDS with optimal k of RF via gap statistic and clusters via PAM (euc)
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="RF gap statistic optimal K and clusters (PAM clustering)", col=cmds_nbclust_clust$clustering, pch=16)

writeLines("\nDone with gap statistic!\n")

writeLines("Doing cMDS with Hierarchical Clustering...\n")
# determine optimal k from cMDS via hierarchical clustering with BIC
# adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
cmdsplot_clust <- Mclust(cmdsplot2$points, G = 1:9)
mclust_grps_cmdsplot <- as.numeric(cmdsplot_clust$classification)
optimalK.max.cmds <- max(mclust_grps_cmdsplot)
# cMDS with optimal k and clusters of RF via hierarchical clustering
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS RF optimal K and clusters (hierarchical clustering)", col=mclust_grps_cmdsplot, pch=16)
mclust.grps.cmdsplot <- mclust_grps_cmdsplot

writeLines("Done with hierarchical clustering!\n\n")

############################
# Write groups to file for RF proximity scores, cmds, gap statistic and hierarchical clustering optimal K

writeLines("\nWriting cMDS groups to file...\n")

# Proximity Scores
write.table(file = paste0(opt$prefix, "_cmds_pam_proximityScores_K", cmdsproxK[[2]], "_run", opt$run, ".txt"), x = cmdsproxK[[1]]$clustering, quote = F, sep = "\t", col.names = F)

# cmds
write.table(file = paste0(opt$prefix, "_cmds_grps_pam_K", cmdsK[[2]], "_run", opt$run, ".txt"), x = cmdsK[[1]]$clustering, quote = F, sep = "\t", col.names = F)

# Gap statistic
write.table(file = paste0(opt$prefix, "_cmds_pam_gapstatisticK_run", opt$run, ".txt"), x = cmds_nbclust_clust$clustering, quote = F, sep = "\t", col.names = 
F)

# Hierarchical
write.table(file = paste0(opt$prefix, "_cmds_RFK_hierarchicalClustering_grps_run", opt$run, ".txt"), x = mclust_grps_cmdsplot, quote = F, sep = "\t", col.names = F)

s.label(cmdsplot2$points, clabel=0.5, grid=0)

writeLines("Done with cMDS!\n")

rm(cmdsplot_bstick)
rm(cmdsplot2)
rm(cmdsplot1)
rm(cmdsproxK)
rm(cmdsplot_clust)
rm(mclust_grps_cmdsplot)
rm(cmds_nbclust)
rm(cmds_nbclust_clust)
rm(cmds_nbclust_k)
gc()

###############
# isotonic MDS
###############

writeLines("\nDoing isoMDS...\n")

# isoMDS
isomdsplot <- isoMDS(1-rftest$proximity, maxit = 10000, tol = 0.00001)
# "The output of cmdscale on 1 - rf$proximity is returned invisibly" (MDSplot documentation)
plot(isomdsplot$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="isoMDS DAPC optimal K and clusters", col=results$grp, pch=16)

writeLines("\nDoing isoMDS with PAM clustering...\n")

# pam clustering with optimal k from DAPC
isomdsK <- pamK(isomdsplot$points, 9)


# mean silhouette width of k
# can test multiple k values and select k with highest value
write.table(x = mean(silhouette(isomdsK[[1]])[, "sil_width"]), file = paste0(opt$prefix, "_isomds_pamClusters_K", isomdsK[[2]], "_meansilhouette_run", opt$run, ".txt"))

# idoMDS with optimal k of DAPC and clusters via PAM
plot(isomdsplot$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main=paste0("isoMDS K=", isomdsK[[2]], " and clusters (PAM clustering)"), col=isomdsK[[1]]$clustering, pch=16)

writeLines("\nDone!\n")

writeLines("Doing isoMDS with Gap Statistic...\n")

# determine optimal k using gap statistic
# can adjust k.max
isomds_nbclust <- fviz_nbclust(isomdsplot$points, kmeans, nstart=25, iter.max=50, method = "gap_stat", k.max = 9, nboot = 1000) + labs(subtitle = "Gap statistic method")
isomds_nbclust
isomds_nbclust_k <- isomds_nbclust[["layers"]][[4]][["data"]][["xintercept"]]
# pam clustering with optimal k from gap statistic
isomds_nbclust_clust2 <- pam(rftest$proximity, isomds_nbclust_k)
# isoMDS with optimal k of RF via gap statistic and clusters via PAM
plot(isomdsplot$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="isoMDS RF gap statistic optimal K and clusters (PAM clustering)", col=isomds_nbclust_clust2$clustering, pch=16)

writeLines("Done!\n")

writeLines("Doing isoMDS hierarchical clustering...\n\n")
# determine optimal k of RF via hierarchical clustering with BIC
# adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
isomdsplot_clust <- Mclust(isomdsplot$points, G = 1:9)
mclust_grps_isomdsplot2 <- as.numeric(isomdsplot_clust$classification)
max(mclust_grps_isomdsplot2)
# isoMDS with optimal k and clusters of RF via hierarchical clustering
plot(isomdsplot$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="isoMDS RF optimal K and clusters (hierarchical clustering)", col=mclust_grps_isomdsplot2, pch=16)

s.label(isomdsplot$points, clabel=0.5, grid=0)

writeLines("Done!\n")

# Write groups to file for RF proximity scores, cmds, gap statistic and hierarchical clustering optimal K

writeLines("Writing isoMDS groups to file...")

# isoMDS
write.table(file = paste0(opt$prefix, "_isomds_grps_pam_K", isomdsK[[2]], "_run", opt$run, ".txt"), x = isomdsK[[1]]$clustering, quote = F, sep = "\t", col.names = F)

# isoMDS Gap statistic
write.table(file = paste0(opt$prefix, "_isomds_grps_pam_gapstatisticK_run", opt$run, ".txt"), x = isomds_nbclust_clust2$clustering, quote = F, sep = "\t", col.names = F)


# isoMDS Hierarchical
write.table(file = paste0(opt$prefix, "_isomds_grps_hierarchicalK_run", opt$run, ".txt"), x = mclust_grps_isomdsplot2, quote = F, sep = "\t", col.names = F)

writeLines("Done!")

dev.off()

rm(isomdsplot)
rm(mclust_grps_isomdsplot2)
rm(isomdsK)
rm(isomds_nbclust)
rm(isomds_nbclust_clust2)
rm(isomdsplot_clust)
rm(rftest)
gc()

###############################################
###############################################
# t-SNE
###############################################
###############################################

writeLines("\nDoing t-SNE Grid Search...\n")

##### perplexity=50 #####

# prepare plot labels and such
# this assumes you have a population assignment (a priori species) column in the .str file.
colors = rainbow(length(unique(data$pop)))
names(colors) = unique(data$pop)
ecb = function(x,y){plot(x,t='n'); text(x, labels=data$pop, col=colors[data$pop])}
# OR
# this makes it so it is grouped by DAPC clusters
#colors = rainbow(length(unique(results$grp)))
#names(colors) = unique(results$grp)
#ecb = function(x,y){plot(x,t='n'); text(x, labels=results$grp, col=colors[results$grp])}

# t-SNE on principal components of scaled data
# adjust perplexity, initial_dims
# can do k=3 for 3D plot
# should do only <50 variables
# can do it on pca$li (if you reduce the number of components), or on cmdsplot2$points

for (i in seq(5, 50, 5)){
  gc()
  pdf(file = paste0(opt$prefix, "tsne_plots_p", i, "_run", opt$run, ".pdf"), onefile = T)
  writeLines(paste0("\nDoing t-SNE with Perplexity ", i, "...\n"))
  tsne_p50 <- tsne(pca1$li, 
                epoch_callback=ecb, 
                max_iter=20000, 
                perplexity=i, 
                initial_dims=5)
  
  writeLines("\nDone!\n")

  writeLines(paste0("\nDoing PAM clustering for t-SNE Perplexity", i, "\n"))
  # pam clustering with optimal k from DAPC
  tsnep50K <- pamK(tsne_p50, 9)
  
  writeLines("\nDone!\n")
  
  # Plot PAM clustering results
  plot(tsne_p50, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main=paste0("t-SNE p50 K=", tsnep50K[[2]], " optimal K and clusters (PAM clustering)"), col=tsnep50K[[1]]$clustering, pch=16)
  
  # mean silhouette width of k
  # can test multiple k values and select k with highest value
  write.table(file = paste0(opt$prefix, "_tsne_p", i, "_meansilhouette_K", tsnep50K[[2]], "_run", opt$run, ".txt"), x = mean(silhouette(tsnep50K[[1]])[, "sil_width"]), quote = F, sep = "\t", col.names = F, row.names = F)
  
  writeLines(paste0("\nDoing Gap Statistic for t-SNE Perplexity ", i, "...\n"))
  # clustering for perplexity=50
  # determine optimal k using gap statistic
  # can adjust k.max
  tsne_p50_nbclust <- fviz_nbclust(tsne_p50, kmeans, nstart=25, iter.max=50, method = "gap_stat", k.max = 9, nboot = 1000) + labs(subtitle = paste0("Gap statistic method p", i))
  tsne_p50_nbclust
  tsne_p50_nbclust_k <- tsne_p50_nbclust[["layers"]][[4]][["data"]][["xintercept"]]
  
  writeLines("\nDone!\n")
  
  # pam clustering with optimal k from gap statistic
  tsne_p50_nbclust_clust <- pam(tsne_p50, tsne_p50_nbclust_k)
  
  # t-SNE with optimal k of RF via gap statistic and clusters via PAM
  plot(tsne_p50, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main=paste0("t-SNE p", i, " RF gap statistic optimal K and clusters (PAM clustering)"), col=tsne_p50_nbclust_clust$clustering, pch=16)
  
  writeLines(paste0("Doing Hierarchical clustering for t-SNE Perlexity ", i, "...\n"))
  
  # determine optimal k of RF via hierarchical clustering with BIC
  # adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
  tsne_p50_clust <- Mclust(tsne_p50, G=1:9)
  mclust_grps_tsne_p50 <- as.numeric(tsne_p50_clust$classification)
  max(mclust_grps_tsne_p50)
  
  # t-SNE p5 with optimal k and clusters of RF via hierarchical clustering
  plot(tsne_p50, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main=paste0("t-SNE p", i, " RF optimal K and clusters (hierarchical clustering)"), col=mclust_grps_tsne_p50, pch=16)
  
  writeLines("\nDone!\n")
  
  writeLines(paste0("\nWriting t-SNE Perplexity ", i, " groups to file...\n"))
  
  # Save groupings to files
  # t-sne perplexity = 50 PAM clusetering
  write.table(file = paste0(opt$prefix, "_tsnep", i, "_grps_pam_K", tsnep50K[[2]], "_run", opt$run, ".txt"), x = tsnep50K[[1]]$clustering, quote = F, sep = "\t", col.names = F)
  
  # t-sne P5 Gap statistic groups
  write.table(file = paste0(opt$prefix, "_tsnep", i, "_grps_gapstatisticK_pam_run", opt$run, ".txt"), x = tsne_p50_nbclust_clust$clustering, quote = F, sep = "\t", col.names = F)
  
  write.table(x = mclust_grps_tsne_p50, file = paste0(opt$prefix, "_tsnep", i, "_grps_hierarchcialK_run", opt$run, ".txt"))
  
  writeLines("\nGroups written to file!\n")
  
  writeLines(paste0("\nDone with t-SNE perplexity ", i, "!\n"))
  dev.off()
  gc()
}

writeLines("Done with all t-SNE perplexities!\n")

dir.create(path = "./Robjects_psearch", showWarnings = FALSE)

writeLines("\nSaving R objects...\n")
save.image(file = file.path("./Robjects_psearch", paste0(opt$prefix, "psearch_run", opt$run, ".RData")))
writeLines("\nSaved R objects to RData files!\n")

gc()

writeLines("Done with run!")
