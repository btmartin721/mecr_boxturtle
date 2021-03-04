########################################################
## Script by Bradley T. Martin
## This script makes a structure-like barplot for each
## UML method and for each dataset and t-SNE perplexity.
## It then uses gridExtra::grid.arrange to put them all 
## onto one page.
########################################################

########################################################
# Tyler K. Chafin wrote a separate script, UML_alignK.R, 
# to use first. Or in this case for varyingly filtered
# datasets use the UML_alignK_missDataRuns_maf.R script
# (slightly modified by Bradley T. Martin).
# UML_alignK.R uses the CLUMPP algorithm
# (a heuristic search) to permute all groups and 
# standardize the group IDs. 
# It loads the data into a list of two 
# data frames. The first has the group assignments.
# The second has the assignment probabilities for each K.
########################################################


#########################################################
### A) Installing and loading required packages
#########################################################

# For using mixedsort with files
if (!require("gtools")) {
  install.packages("gtools", dependencies = TRUE)
}
library("gtools")

# For making colorful plots
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
}
library("RColorBrewer")

# For doing full join and manipulating data.
if (!require("dplyr")) {
  install.packages("dplyr")
}
library("dplyr")

# For making the plots
if (!require("ggplot2")) {
  install.packages("ggplot2")
}
library("ggplot2")

# Required for using melt method
# to make a long-format data.frame.
if (!require("reshape2")) {
  install.packages("reshape2")
}
library("reshape2")

# Required for using melt method
if (!require("plotrix")) {
  install.packages("plotrix")
}
library("plotrix")

# Required for using count method
if (!require("plyr")) {
  install.packages("plyr")
}
library("plyr")

# To manipulate/ change the order of factors.
if (!require("forcats")) {
  install.packages("forcats")
}
library("forcats")

# To plot the barplots onto one page.
if (!require("gridExtra")){
  install.packages("gridExtra")
}
library("gridExtra")

# For making the grid.arrange title larger.
if (!require("grid")){
  install.packages("grid")
}
# library("grid") # Don't have to load. Called with grid::

# beepr just makes victory or error noises when the 
# code execution completes. Can be removed if you don't want
# that (if so, also remove calls beepr::beep(3) calls below.
if (!require("beepr")){
  install.packages("beepr")
}
library("beepr")

########################################################
### B) Main functions to make plots
########################################################

makeplots <- function(joined, 
                      myoutfile.prefix, 
                      uml_method, 
                      sample.ids,
                      myCols=NULL,
                      plotDIR="./plots") {
  ## Function to make heatmaps of each UML run.
  ## Different from the barplots.
  ## joined: 
      # First list element from object output
      # from UML_alignK_missDataRuns_maf.R
  ## myoutfile.prefix: 
      # prefix for output file (Character string)
  ## uml_method:
      # Dimension reduction method (e.g. cMDS)
      # (Character string)
  ## sample.ids: 
      # data.frame with two columns (V1 and V2)
      # 1. V1: factor of sample ids (character strings).
      # 2. V2: sample order number (integers).
  ## plotDIR:
      # Directory to save plots to (Character string).
  
  # Change column names to sequential integers
  # Did this because full_join makes really long colnames
  runs <- c("sampleID", seq(1, length(joined) - 1))
  colnames(joined) <- runs

  # Join sample ids to get order that I want.
  joined <- merge(joined, 
                  sample.ids, 
                  by.x = "sampleID", 
                  by.y = "V1", 
                  all.x = TRUE)
  
  # Reorder based on column V2.
  joined <- joined[order(joined$V2),]

  # Did this to get the sample ids in the order I want.
  # Reorders factor levels by sorting along another variable.
  joined$sampleID <- forcats::fct_reorder(joined$sampleID, 
                                          joined$V2, 
                                          min)

  # Remove last column.
  joined <- joined[,1:(length(joined)-1)]

  #########################################################
  ### Customizing and plotting the heat map
  #########################################################

  # Create output directory for plots if it doesn't exist.
  dir.create(plotDIR, showWarnings = FALSE)
  
  # Melt the dataframe. For making categorical heatmap
  mymelt <- reshape2::melt(joined, id.var = "sampleID")
  
  if (!is.null(myCols)){
  
  # Plot it
  p <-
    ggplot(mymelt, aes(variable, sampleID)) + 
    geom_tile(aes(fill = as.factor(value)),
      colour = "white") + scale_fill_manual(
        values = myCols
    ) +
    theme(
      axis.text.y = element_text(size = 8, colour = "black"),
      axis.text.x = element_blank(),
      axis.title = element_text(size = 14, colour = "black"),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) + guides(fill = guide_legend(title = "Species")) +
    labs(x = "Run", y = "Sample ID", title = uml_method)

  } else {
    
    # Plot it
    p <-
      ggplot(mymelt, aes(variable, sampleID)) + 
      geom_tile(aes(fill = as.factor(value)),
                colour = "white") + scale_fill_brewer(
                  palette = "Set1"
                ) +
      theme(
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 14, colour = "black"),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) + guides(fill = guide_legend(title = "Species")) +
      labs(x = "Run", y = "Sample ID", title = uml_method)
    
  }

  # Save it to PDF file
  ggsave(filename = paste0(myoutfile.prefix, ".pdf"),
         plot = p,
         device = "pdf", dpi = 300, path = plotDIR, height = 7, width = 7)

  # Save it to PNG file
  ggsave(filename = paste0(myoutfile.prefix, ".png"),
         plot = p,
         device = "png", dpi = 300, path = plotDIR, height = 7, width = 7)
  
  # Convert all but first column to numeric.
  joined[,-1] <- lapply(joined[,-1], function(x) as.numeric(as.character(x)))
  
  return(joined)
}

sum_plots <- function(propsdf, 
                      sample.ids, 
                      pops, 
                      title, 
                      subtitle, 
                      notfirst=FALSE, 
                      myCols=NULL){
  ## Function to summarize K-aligned UML output
  ## and plot it as a structure-like barplot.
  ## propsdf:
      # List element 2 from .rds object output from 
      # UML_alignK_missDataRuns_maf.R
  ## sample.ids
      # data.frame with two columns (V1 and V2)
      # 1. V1: Character factor of sample ids.
      # 2. V2: sample order number (integers).
  ## pops:
      # Two-column tab-separated population map file.
      # First column should be sample ids. colname="Sample"
          # Should be factor.
      # Second column is population ids. colname="Pop"
          # Should be factor.
  ## title:
      # Title of barplot (character string).
      # For one barplot. E.g., cMDS
  ## subtitle
      # subtitle of barplot (character string).
      # For only one clustering algorithm (e.g., GS)
  ## notfirst: 
      # Logical.
      # FALSE if first plot
      # TRUE if not first plot.
      # Default = FALSE
      # Purpose: to plot barplots on a grid.
  
  # Join sample ids to get order that I want.
  joined_tmp <- merge(propsdf, 
                      sample.ids, 
                      by.x = "Sample", 
                      by.y = "V1", 
                      all.x = TRUE)
  
  # Join popmap with propsdf and sample.ids
  joined <- merge(joined_tmp, pops, by="Sample", all.x = TRUE)
  
  # Reorder them based on column V2.
  joined <- joined[order(joined$V2),]
  
  # Did this to get the sample ids in the order I want.
  # Reorders factor levels by sorting along another column.
  joined$Sample <- forcats::fct_reorder(joined$Sample,
                                        joined$V2, 
                                        min)
  
  joined$Pop <- forcats::fct_reorder(joined$Pop, 
                                     joined$V2, 
                                     min)
  
  # Create dummy variable to facet on: 
  # this name will appear in the strip
  joined$tempvar <- title
  
  # Melt data.frame into long format.
  mymelt <- reshape2::melt(joined, 
                           id.var = c("Sample", 
                                      "Pop", 
                                      "V2", 
                                      "tempvar"))
  
  # If first plot.
  if (!notfirst){
    
    if (!is.null(myCols)){
    
    p <- ggplot(mymelt, aes(value, 
                            Sample, 
                            fill = factor(variable))) +
      geom_col(width=1) +
      facet_grid(Pop~., 
                 scales = "free", 
                 space = "free", 
                 switch="y") +
      theme_minimal() + 
      labs(y = "Population", 
           size = 18, 
           title = title, 
           subtitle = subtitle, 
           colour = "black") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_discrete(expand = expand_scale(add = 1)) +
      theme(
        panel.spacing.y = unit(0.3, "lines"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", 
                                    size = 24), 
        plot.title = element_text(colour = "black", 
                                  size = 18),
        plot.subtitle = element_text(colour = "black", 
                                     size = 18),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill="#20b2aa"),
        strip.text = element_text(colour = "black", 
                                  size = 14),
        strip.placement = "outside"
      ) +
      scale_fill_manual(values = myCols)
    
    }else {
      
      p <- ggplot(mymelt, aes(value, 
                              Sample, 
                              fill = factor(variable))) +
        geom_col(width=1) +
        facet_grid(Pop~., 
                   scales = "free", 
                   space = "free", 
                   switch="y") +
        theme_minimal() + 
        labs(y = "Population", 
             size = 18, 
             title = title, 
             subtitle = subtitle, 
             colour = "black") +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_discrete(expand = expand_scale(add = 1)) +
        theme(
          panel.spacing.y = unit(0.3, "lines"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(colour = "black", 
                                      size = 24), 
          plot.title = element_text(colour = "black", 
                                    size = 18),
          plot.subtitle = element_text(colour = "black", 
                                       size = 18),
          panel.grid = element_blank(),
          legend.position = "none",
          strip.background = element_rect(fill="#20b2aa"),
          strip.text = element_text(colour = "black", 
                                    size = 14),
          strip.placement = "outside"
        ) +
        scale_fill_brewer(palette = "Set1")
    }
  
  } else {
    
    if (!is.null(myCols)){
      
      p <- ggplot(mymelt, aes(value, 
                              Sample, 
                              fill = factor(variable))) +
        geom_col(width=1) +
        facet_grid(Pop~., 
                   scales = "free", 
                   space = "free", 
                   switch="y") +
        theme_minimal() + 
        labs(size = 20, 
             title = title, 
             subtitle = subtitle) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_discrete(expand = expand_scale(add = 1)) +
        theme(
          panel.spacing.y = unit(0.3, "lines"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(colour = "black", 
                                    size = 18), 
          plot.subtitle = element_text(colour = "black", 
                                       size = 18),
          panel.grid = element_blank(),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank()
        ) +
        scale_fill_manual(values = myCols)
      
    } else {
    
    p <- ggplot(mymelt, aes(value, 
                            Sample, 
                            fill = factor(variable))) +
      geom_col(width=1) +
      facet_grid(Pop~., 
                 scales = "free", 
                 space = "free", 
                 switch="y") +
      theme_minimal() + 
      labs(size = 20, 
           title = title, 
           subtitle = subtitle) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_discrete(expand = expand_scale(add = 1)) +
      theme(
        panel.spacing.y = unit(0.3, "lines"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(colour = "black", 
                                  size = 18), 
        plot.subtitle = element_text(colour = "black", 
                                     size = 18),
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()
      ) +
      scale_fill_brewer(palette="Set1")
    }
  }
  
  return(p)
}

#########################################################
### C) Reading in data and transform it into matrix format
#########################################################
sample.ids <- 
  read.table(file = "sampleids_numbered.txt", 
                         header = F)

# Directory to save plots to.
plotDIR <- "plots_maf/"

# Create directory to save plots to if it doesn't
# already exist.
dir.create(plotDIR, 
           showWarnings = FALSE)

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

# Read in population map file.
# Two-column tab-separated file.
# SampleIDs = column 1
# Population IDs = column 2.
popmap <-
  read.table(
    file = "maf01.popmap",
    header = FALSE,
    col.names = c("Sample", "Pop")
  )

# Save all plots to one PDF file.
# If you want to save them all to separate
# files, use pdf() just before doing the 
# for loop for t-SNE perplexities.

###################################################
## Make the plots
###################################################

# TRUE if saving all barplots to one PDF.
onefile <- FALSE

if (onefile){
  pdf(
    file = file.path(plotDIR, 
                     "uml_summary_all.pdf"),
    width = 20,
    height = 16,
    onefile = TRUE
  )
}

i <- 50
j <- 25
p <- 15
m <- "0.05" # minor allele frequencies.

writeLines(paste0("\nDoing missInd ",
                  i,
                  " Pop ",
                  j,
                  " MAF ",
                  m,
                  "..."))

# Load aligned VAE objects from
# UML_alignK_missDataRuns_maf.R
vae <-
  readRDS(file = paste0("Robjects_maf/vae",
                        ".rds"))

cmds_gapstat <-
  readRDS(file = paste0(
    "Robjects_maf/cmds_gapstat",
    ".rds"
  ))

# Read in cMDS and isoMDS objects.
cmds_hier <-
  readRDS(file = paste0(
    "Robjects_maf/cmds_hier",
    ".rds"
  ))

cmds_pam <-
  readRDS(file = paste0(
    "Robjects_maf/cmds_pam",
    ".rds"
  ))

cmds_pam_prox <-
  readRDS(file = paste0(
    "Robjects_maf/cmds_pam_prox",
    ".rds"
  ))

isomds_gapstat <-
  readRDS(file = paste0(
    "Robjects_maf/isomds_gapstat",
    ".rds"
  ))

isomds_hier <-
  readRDS(file = paste0(
    "Robjects_maf/isomds_hier",
    ".rds"
  ))

isomds_pam <-
  readRDS(file = paste0(
    "Robjects_maf/isomds_pam",
    ".rds"
  ))

#*********************************************************
### Set plot titles and prefixes here.
#*********************************************************
vae.prefix <- "vae"

vae.title <- "VAE DBSCAN"

cmds.gapstat.prefix <- "cmds_gapstat"

cmds.gapstat.title <- "cMDS Gap Statistic"

# cMDS Hierarchical Clustering
cmds.hier.prefix <- "cmds_hierarchical"

cmds.hier.title <- "cMDS Hierarchical Clustering"

# cMDS PAM Clustering (cMDS Groups)
cmds.pam.cmdsgroups.prefix <- "cmds_pam"

cmds.pam.cmdsgroups.title <- 
  "cMDS PAM Clustering (cMDS Groups)"

# cMDS PAM Clustering (Proximity Scores)
cmds.pam.prox.prefix <- "cmds_pam_prox"

cmds.pam.prox.title <- "cMDS PAM Clustering (Proximity Scores)"

# isoMDS Gap Statistic (PAM Clustering)
isomds.gapstat.prefix <- "isomds_gapstat"

isomds.gapstat.title <- "isoMDS Gap Statistic"

# isoMDS Hierarchical Clustering
isomds.hier.prefix <- "isomds_hierarchical"

isomds.hier.title <- "isoMDS Hierarchical Clustering"

# isoMDS PAM Clustering
isomds.pam.prefix <- "isomds_pam"

isomds.pam.title <- "isoMDS PAM Clustering"

# Make heatmaps of all UML runs.
vae.df <-
  makeplots(vae[[1]],
            vae.prefix,
            vae.title,
            sample.ids,
            myCols = c25)

cmds.gapstat.df <-
  makeplots(cmds_gapstat[[1]],
            cmds.gapstat.prefix,
            cmds.gapstat.title,
            sample.ids,
            myCols = c25)

cmds.hier.df <-
  makeplots(cmds_hier[[1]],
            cmds.hier.prefix,
            cmds.hier.title,
            sample.ids, 
            myCols = c25)

cmds.pam.cmdsgroups.df <-
  makeplots(cmds_pam[[1]],
            cmds.pam.cmdsgroups.prefix,
            cmds.pam.cmdsgroups.title,
            sample.ids, 
            myCols = c25)

cmds.pam.prox.df <-
  makeplots(cmds_pam_prox[[1]],
            cmds.pam.prox.prefix,
            cmds.pam.prox.title,
            sample.ids,
            myCols = c25)

isomds.gapstat.df <-
  makeplots(isomds_gapstat[[1]],
            isomds.gapstat.prefix,
            isomds.gapstat.title,
            sample.ids,
            myCols = c25)

isomds.hier.df <-
  makeplots(isomds_hier[[1]],
            isomds.hier.prefix,
            isomds.hier.title,
            sample.ids,
            myCols = c25)

isomds.pam.df <-
  makeplots(isomds_pam[[1]],
            isomds.pam.prefix,
            isomds.pam.title,
            sample.ids,
            myCols = c25)

###############################################
## Make the barplots for cMDS, isoMDS, and VAE
###############################################
p1 <-
  sum_plots(cmds_gapstat[[2]],
            sample.ids,
            popmap,
            "cMDS",
            "GS",
            myCols = c25)

p2 <-
  sum_plots(cmds_hier[[2]],
            sample.ids,
            popmap,
            "cMDS",
            "H",
            TRUE,
            myCols = c25)

p3 <-
  sum_plots(cmds_pam[[2]],
            sample.ids,
            popmap,
            "cMDS",
            "PAM",
            TRUE, 
            myCols = c25)

p4 <-
  sum_plots(cmds_pam_prox[[2]],
            sample.ids,
            popmap,
            "Prox",
            "PAM",
            TRUE, 
            myCols = c25)

p5 <-
  sum_plots(isomds_gapstat[[2]],
            sample.ids,
            popmap,
            "isoMDS",
            "GS",
            TRUE,
            myCols = c25)

p6 <-
  sum_plots(isomds_hier[[2]],
            sample.ids,
            popmap,
            "isoMDS",
            "H",
            TRUE, 
            myCols = c25)
p7 <-
  sum_plots(isomds_pam[[2]],
            sample.ids,
            popmap,
            "isoMDS",
            "PAM",
            TRUE, 
            myCols = c25)

p.vae <-
  sum_plots(vae[[2]],
            sample.ids,
            popmap,
            "VAE",
            "DBSCAN",
            TRUE, 
            myCols = c25)

# To save each dataset to separate files.
# Only if onefile == FALSE
if (!onefile) {
  pdf(
    file = file.path(
      plotDIR,
      "uml_summary.pdf"
    ),
    width = 20,
    height = 16,
    onefile = TRUE
  )
}

# t-SNE perplexity grid search.
title <- paste0("tSNE",
                  " P",
                  p)
  
writeLines(paste0("Doing t-SNE perplexity ",
                    p,
                    "...\n"))
  
  tsne_gapstat <-
    readRDS(
      file = paste0(
        "Robjects_maf/tsne_gapstat",
        "_P",
        p,
        ".rds"
      )
    )
  
  tsne_hier <-
    readRDS(
      file = paste0(
        "Robjects_maf/tsne_hier",
        "_P",
        p,
        ".rds"
      )
    )
  
  tsne_pam <-
    readRDS(file = paste0(
      "Robjects_maf/tsne_pam",
      "_P",
      p,
      ".rds"
    ))
  
  #*********************************************************
  ### Set plot titles and prefixes here.
  #*********************************************************
  tsnep5.gapstat.prefix <-
    paste0("tsneP15_gapstat",
           "_P",
           p)
  
  tsnep5.gapstat.title <-
    paste0("t-SNE Gap Statistic",
           "(Perplexity = ",
           p,
           ")")
  
  # tsne perplexity = 5 Hierarchical Clustering
  tsnep5.hier.prefix <-
    paste0("tsneP15_hierarchical",
           "_P",
           p)
  
  tsnep5.hier.title <-
    paste0(
      "t-SNE Hierarchical Clustering",
      "(Perplexity = ",
      p,
      ")"
    )
  
  tsnep5.pam.prefix <- "tsneP15_pam"
  
  tsnep5.pam.title <-
    paste0("t-SNE PAM Clustering",
           "(Perplexity = ",
           p,
           ")")
  
  ##################################
  ### Make the t-SNE plots
  ##################################
  
  # Heatmaps for all runs per method.
  tsnep50.gapstat.df <-
    makeplots(tsne_gapstat[[1]],
              tsnep5.gapstat.prefix,
              tsnep5.gapstat.title,
              sample.ids, 
              myCols = c25)
  
  tsnep50.hier.df <-
    makeplots(tsne_hier[[1]],
              tsnep5.hier.prefix,
              tsnep5.hier.title,
              sample.ids,
              myCols = c25)
  
  tsnep50.pam.df <-
    makeplots(tsne_pam[[1]],
              tsnep5.pam.prefix,
              tsnep5.pam.title,
              sample.ids, 
              myCols = c25)
  
  # Make the structure-like barplots.
  p8 <-
    sum_plots(tsne_gapstat[[2]],
              sample.ids,
              popmap,
              paste0("t-SNE P", p),
              "GS",
              TRUE, 
              myCols = c25)
  
  p9 <-
    sum_plots(tsne_hier[[2]],
              sample.ids,
              popmap,
              paste0("t-SNE P", p),
              "H",
              TRUE,
              myCols = c25)
  
  p10 <-
    sum_plots(tsne_pam[[2]],
              sample.ids,
              popmap,
              paste0("t-SNE P", p),
              "PAM",
              TRUE,
              myCols = c25)
  
  # Save onto one page (a grid)
  g <-
    grid.arrange(
      p1,
      p2,
      p3,
      p4,
      p5,
      p6,
      p7,
      p8,
      p9,
      p10,
      p.vae,
      nrow = 1,
      ncol = 11,
      top = grid::textGrob(title,
                           gp = grid::gpar(fontsize = 15,
                                           font = 8))
    )
  print(g)

if (!onefile) {
  dev.off()
}

if (onefile) {
  dev.off()
}
beepr::beep(3)
