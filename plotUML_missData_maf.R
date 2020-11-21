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
  
  # Plot it
  p <-
    ggplot(mymelt, aes(variable, sampleID)) + 
    geom_tile(aes(fill = as.factor(value)),
      colour = "white") + scale_fill_brewer(
      palette="Set1"
    ) +
    theme(
      axis.text.y = element_text(size = 8, colour = "black"),
      axis.text.x = element_blank(),
      axis.title = element_text(size = 14, colour = "black"),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) + guides(fill = guide_legend(title = "Species")) +
    labs(x = "Run", y = "Sample ID", title = uml_method)

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
                      notfirst=FALSE){
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
  joined$Sample <- forcats::fct_reorder(joined$Sample, joined$V2, min)
  joined$Pop <- forcats::fct_reorder(joined$Pop, joined$V2, min)
  
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
      scale_fill_brewer(palette="Set1")
  }
  
  else {
    
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
  
  return(p)
}

#########################################################
### C) Reading in data and transform it into matrix format
#########################################################
sample.ids <- 
  read.table(file = file.path("popmaps", 
                              "sampleids_numbered.txt"), 
                         header = F)

# Directory to save plots to.
plotDIR <- "plots_maf_final2"

# Create directory to save plots to if it doesn't
# already exist.
dir.create(file.path("missDataRuns", 
                     plotDIR), 
           showWarnings = FALSE)

# Make color vector
colors <- brewer.pal(n = 9, name = "Set1")

# Read in population map file.
# Two-column tab-separated file.
# SampleIDs = column 1
# Population IDs = column 2.
popmap <-
  read.table(
    file = "popmaps/phylogen.popmap.final.txt",
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
    file = file.path("missDataRuns", 
                     plotDIR, 
                     "uml_summary_all.pdf"),
    width = 20,
    height = 16,
    onefile = TRUE
  )
}
mafs <- c("0.0", "0.01", "0.03", "0.05") # minor allele frequencies.
for (i in seq(25, 100, 25)) { # individual missing data %
  for (j in seq(25, 100, 25)) { # population missing data %
    for (m in seq_along(mafs)) { # minor allele frequencies
      
      writeLines(paste0("\nDoing missInd ", 
                        i, 
                        " Pop ", 
                        j, 
                        " MAF ", 
                        mafs[m], 
                        "..."))
      
      # This dataset only had 9 SNPs after using the MAF filters.
      # The UML runs failed accordingly.
      # Thus, it is skipped here.
      if (i == 100 & j == 25){
        next
      }
      
      # Load aligned VAE objects from 
      # UML_alignK_missDataRuns_maf.R
      vae <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/vae_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            mafs[m],
            ".rds"
          )
        )
      
      cmds_gapstat <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/cmds_gapstat_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            mafs[m],
            ".rds"
          )
        )
      
      # Read in cMDS and isoMDS objects.
      cmds_hier <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/cmds_hier_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            mafs[m],
            ".rds"
          )
        )
      
      cmds_pam <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/cmds_pam_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            mafs[m],
            ".rds"
          )
        )
      
      cmds_pam_prox <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/cmds_pam_prox_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            mafs[m],
            ".rds"
          )
        )
      
      isomds_gapstat <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/isomds_gapstat_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            mafs[m],
            ".rds"
          )
        )
      
      isomds_hier <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/isomds_hier_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            mafs[m],
            ".rds"
          )
        )
      
      isomds_pam <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/isomds_pam_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            mafs[m],
            ".rds"
          )
        )
      
      #*********************************************************
      ### Set plot titles and prefixes here.
      #*********************************************************
      vae.prefix <-
        paste0("vae_missInd", 
               i, 
               "_Pop",
               j, 
               "_maf", 
               mafs[m])
      
      vae.title <-
        paste0("VAE DBSCAN missInd ", 
               i, 
               " Pop ", 
               j, 
               "_maf", 
               mafs[m])
      
      cmds.gapstat.prefix <-
        paste0("cmds_gapstat_missInd", 
               i, 
               "_Pop", 
               j, 
               "_maf", 
               mafs[m])
      
      cmds.gapstat.title <-
        paste0("cMDS Gap Statistic missInd ", 
               i, 
               " Pop ", 
               j, 
               "_maf", 
               mafs[m])
      
      # cMDS Hierarchical Clustering
      cmds.hier.prefix <-
        paste0("cmds_hierarchical_missInd", 
               i, 
               "_Pop", 
               j, 
               "_maf", 
               mafs[m])
      
      cmds.hier.title <-
        paste0("cMDS Hierarchical Clustering missInd ",
               i,
               " Pop ",
               j,
               "_maf",
               mafs[m])
      
      # cMDS PAM Clustering (cMDS Groups)
      cmds.pam.cmdsgroups.prefix <-
        paste0("cmds_pam_missInd", 
               i, 
               "_Pop", 
               j, 
               "_maf", 
               mafs[m])
      
      cmds.pam.cmdsgroups.title <-
        paste0("cMDS PAM Clustering (cMDS Groups) missInd ",
               i,
               " Pop ",
               j,
               "_maf",
               mafs[m])
      
      # cMDS PAM Clustering (Proximity Scores)
      cmds.pam.prox.prefix <-
        paste0("cmds_pam_prox_missInd", 
               i, 
               "_Pop", 
               j, 
               "_maf", 
               mafs[m])
      
      cmds.pam.prox.title <-
        paste0("cMDS PAM Clustering (Proximity Scores) missInd ",
               i,
               " Pop ",
               j,
               "_maf",
               mafs[m])
      
      # isoMDS Gap Statistic (PAM Clustering)
      isomds.gapstat.prefix <-
        paste0("isomds_gapstat_missInd", 
               i, 
               "_Pop", 
               j, 
               "_maf", 
               mafs[m])
      
      isomds.gapstat.title <-
        paste0("isoMDS Gap Statistic missInd ", 
               i, 
               " Pop ", 
               j, 
               "_maf", 
               mafs[m])
      
      # isoMDS Hierarchical Clustering
      isomds.hier.prefix <-
        paste0("isomds_hierarchical_missInd", 
               i, 
               "_Pop", 
               j, 
               "_maf", 
               mafs[m])
      
      isomds.hier.title <-
        paste0("isoMDS Hierarchical Clustering missInd ",
               i,
               " Pop ",
               j,
               "_maf",
               mafs[m])
      
      # isoMDS PAM Clustering
      isomds.pam.prefix <-
        paste0("isomds_pam_missInd", 
               i, 
               "_Pop", 
               j, 
               "_maf", 
               mafs[m])
      
      isomds.pam.title <-
        paste0("isoMDS PAM Clustering missInd ", 
               i, 
               " Pop ", 
               j, 
               "_maf", 
               mafs[m])
      
      # Make heatmaps of all UML runs.
      vae.df <-
        makeplots(vae[[1]],
                  vae.prefix,
                  vae.title,
                  sample.ids)
      
      cmds.gapstat.df <-
        makeplots(cmds_gapstat[[1]],
                  cmds.gapstat.prefix,
                  cmds.gapstat.title,
                  sample.ids)
      
      cmds.hier.df <-
        makeplots(cmds_hier[[1]],
                  cmds.hier.prefix,
                  cmds.hier.title,
                  sample.ids)
      
      cmds.pam.cmdsgroups.df <-
        makeplots(
          cmds_pam[[1]],
          cmds.pam.cmdsgroups.prefix,
          cmds.pam.cmdsgroups.title,
          sample.ids
        )
      
      cmds.pam.prox.df <-
        makeplots(cmds_pam_prox[[1]],
                  cmds.pam.prox.prefix,
                  cmds.pam.prox.title,
                  sample.ids)
      
      isomds.gapstat.df <-
        makeplots(isomds_gapstat[[1]],
                  isomds.gapstat.prefix,
                  isomds.gapstat.title,
                  sample.ids)
      
      isomds.hier.df <-
        makeplots(isomds_hier[[1]],
                  isomds.hier.prefix,
                  isomds.hier.title,
                  sample.ids)
      
      isomds.pam.df <-
        makeplots(isomds_pam[[1]],
                  isomds.pam.prefix,
                  isomds.pam.title,
                  sample.ids)
      
      ###############################################
      ## Make the barplots for cMDS, isoMDS, and VAE
      ###############################################
      p1 <-
        sum_plots(cmds_gapstat[[2]], 
                  sample.ids, 
                  popmap, 
                  "cMDS", 
                  "GS")
      
      p2 <-
        sum_plots(cmds_hier[[2]], 
                  sample.ids, 
                  popmap, 
                  "cMDS", 
                  "H", 
                  TRUE)
      
      p3 <-
        sum_plots(cmds_pam[[2]], 
                  sample.ids, 
                  popmap, 
                  "cMDS", 
                  "PAM", 
                  TRUE)
      
      p4 <-
        sum_plots(cmds_pam_prox[[2]],
                  sample.ids,
                  popmap,
                  "Prox",
                  "PAM",
                  TRUE)
      
      p5 <-
        sum_plots(isomds_gapstat[[2]],
                  sample.ids,
                  popmap,
                  "isoMDS",
                  "GS",
                  TRUE)
      
      p6 <-
        sum_plots(isomds_hier[[2]],
                  sample.ids,
                  popmap,
                  "isoMDS",
                  "H",
                  TRUE)
      p7 <-
        sum_plots(isomds_pam[[2]],
                  sample.ids,
                  popmap,
                  "isoMDS",
                  "PAM",
                  TRUE)
      
      p.vae <-
        sum_plots(vae[[2]], 
                  sample.ids, 
                  popmap, 
                  "VAE", 
                  "DBSCAN", 
                  TRUE)
      
      # To save each dataset to separate files.
      # Only if onefile == FALSE
      if (!onefile) {
        pdf(
          file = file.path("missDataRuns",
                           plotDIR,
                           paste0("uml_summary_missInd",
                           i,
                           "_Pop", 
                           j,
                           "_maf",
                           mafs[m],
                           ".pdf")),
          width = 20,
          height = 16,
          onefile = TRUE
        )
      }
      
      # t-SNE perplexity grid search.
      for (p in seq(5, 50, 5)) {
        
        title <- paste0("missInd", 
                        i, 
                        " Pop", 
                        j, 
                        " MAF", 
                        mafs[m], 
                        " P", 
                        p)
        
        writeLines(paste0("Doing t-SNE perplexity ", 
                          p, 
                          "...\n"))
        
        tsne_gapstat <-
          readRDS(
            file = paste0(
              "missDataRuns/Robjects_maf/tsne_gapstat_missInd",
              i,
              "_Pop",
              j,
              "_maf",
              mafs[m],
              "_P",
              p,
              ".rds"
            )
          )
        
        tsne_hier <-
          readRDS(
            file = paste0(
              "missDataRuns/Robjects_maf/tsne_hier_missInd",
              i,
              "_Pop",
              j,
              "_maf",
              mafs[m],
              "_P",
              p,
              ".rds"
            )
          )

        tsne_pam <-
          readRDS(
            file = paste0(
              "missDataRuns/Robjects_maf/tsne_pam_missInd",
              i,
              "_Pop",
              j,
              "_maf",
              mafs[m],
              "_P",
              p,
              ".rds"
            )
          )
        
        #*********************************************************
        ### Set plot titles and prefixes here.
        #*********************************************************
        tsnep5.gapstat.prefix <-
          paste0("tsneP50_gapstat_missInd",
                 i,
                 "_Pop",
                 j,
                 "_maf",
                 mafs[m],
                 "_P",
                 p)
        
        tsnep5.gapstat.title <-
          paste0("t-SNE Gap Statistic missInd ",
                 i,
                 " Pop",
                 j,
                 "_maf",
                 mafs[m],
                 "(Perplexity = ",
                 p,
                 ")")
        
        # tsne perplexity = 5 Hierarchical Clustering
        tsnep5.hier.prefix <-
          paste0("tsneP50_hierarchical_missInd",
                 i,
                 "_Pop",
                 j,
                 "_maf",
                 mafs[m],
                 "_P",
                 p)
        
        tsnep5.hier.title <-
          paste0(
            "t-SNE Hierarchical Clustering missInd ",
            i,
            " Pop",
            j,
            "_maf",
            mafs[m],
            "(Perplexity = ",
            p,
            ")"
          )
        
        tsnep5.pam.prefix <-
          paste0("tsneP50_pam_missInd", 
                 i, 
                 "_Pop", 
                 j, 
                 "_maf",
                 mafs[m])
        
        tsnep5.pam.title <-
          paste0("t-SNE PAM Clustering missInd ",
                 i,
                 " Pop",
                 j,
                 "_maf",
                 mafs[m],
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
                    sample.ids)
        
        tsnep50.hier.df <-
          makeplots(tsne_hier[[1]],
                    tsnep5.hier.prefix,
                    tsnep5.hier.title,
                    sample.ids)
        
        tsnep50.pam.df <-
          makeplots(tsne_pam[[1]],
                    tsnep5.pam.prefix,
                    tsnep5.pam.title,
                    sample.ids)
        
        # Make the structure-like barplots.
        p8 <-
          sum_plots(tsne_gapstat[[2]],
                    sample.ids,
                    popmap,
                    paste0("t-SNE P", p),
                    "GS",
                    "notnull")
        
        p9 <-
          sum_plots(tsne_hier[[2]],
                    sample.ids,
                    popmap,
                    paste0("t-SNE P", p),
                    "H",
                    "notnull")
        
        p10 <-
          sum_plots(tsne_pam[[2]],
                    sample.ids,
                    popmap,
                    paste0("t-SNE P", p),
                    "PAM",
                    "notnull")
        
        # Save onto one page (a grid)
        g <-
          grid.arrange(p1,
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
                                            gp=grid::gpar(fontsize=15,
                                                          font=8)))
        print(g)
      }
      if (!onefile){
        dev.off()
      }
    }
  }
}
if (onefile){
  dev.off()
}
beepr::beep(3)

