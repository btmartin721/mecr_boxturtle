##################################################################
## Script by Bradley T. Martin
##
## Plot Terrapene iqtree phylogeny with divergence times alongside 
## UML barplots
## 
## Individuals on the barplots are aligned with the phylogeny tips
##################################################################

library("ggplot2")
library("ggtree")
library("treeio")
library("ape")
library("RColorBrewer")
library("reshape2")
library("grid")

sum_plots <- function(propsdf, sample.ids, pops, title, subtitle){
  # Function to plot UML results as a barplot.
  # propsdf: data.frame with assignment proportions.
  # sample.ids: Two-column tab-separated file with 
    # sampleIDs as column 1 and plotting order as column 2.
  # pops: Two-column tab-separated file with 
    # population map file with sampleIDs as column 1
    # and population ID as column 2.
  # title: Title of barplot (i.e. dimension reduction algorithm).
  # subtitle: Subtitle of barplot (i.e. clustering algorithm).
  
  # Join sample ids to get order that I want.
  joined_tmp <- merge(propsdf, 
                      sample.ids, 
                      by.x = "Sample", 
                      by.y = "V1", 
                      all = TRUE)
  
  joined <- merge(joined_tmp, pops, by="Sample", all.x = TRUE)
  
  joined <- joined[order(joined$V2),]
  
  # Did this to get the sample ids in the order I want.
  joined$Sample <- forcats::fct_reorder(joined$Sample, joined$V2, min)
  joined$Pop <- forcats::fct_reorder(joined$Pop, joined$V2, min)
  
  # Create dummy variable to facet on: this name will appear in the strip
  joined$tempvar <- title
  
  mymelt <- melt(joined, id.var = c("Sample", "Pop", "V2", "tempvar"))
  
  return(mymelt)
}

plot_tree <-
  function(p,
           clades = FALSE,
           OUT = NULL,
           DS = NULL,
           ON = NULL,
           FL = NULL,
           MX = NULL,
           TT = NULL,
           EA = NULL,
           GUFL = NULL,
           GUMS_CH = NULL,
           ONDS = NULL,
           terrapene = NULL,
           east = NULL) {
    # Function to plot a time-calibrated Terrapene phylogeny.
    # p: ggtree object.
    # clades: Logical indicating whether the below arguments are specified.
    # Remaining argument: Specific clades I want labelled. Only plotted
      # if clades == TRUE.
  
  p <- revts(p)
  
  if (!clades) {
    
    p.final <- p +
      # outgroup
      geom_cladelabel(
        node = 423,
        label = "",
        color = c("black", "black"),
        offset = 2.0,
        offset.text = 0.8,
        barsize = 4
      ) +
      
      # DS
      geom_cladelabel(
        node = 388,
        label = "",
        color = c("#08ffff", "black"),
        offset = 2.0,
        offset.text = 0.8,
        barsize = 4
      ) +
      
      # ON
      geom_cladelabel(
        node = 395,
        label = "",
        color = c("#0033cc", "black"),
        offset = 2.0,
        offset.text = 0.8,
        barsize = 4
      ) +
      
      # FL
      geom_cladelabel(
        node = 383,
        label = "",
        color = c("#ffff00", "black"),
        offset = 2.0,
        offset.text = 0.8,
        barsize = 4
      ) +
      
      # MX
      geom_cladelabel(
        node = 381,
        label = "",
        color = c("#ff3333", "black"),
        offset = 2.0,
        offset.text = 0.8,
        barsize = 4
      ) +
      
      # TT
      geom_cladelabel(
        node = 355,
        label = "",
        offset = 2.0,
        offset.text = 0.8,
        color = c("#006600", "black"),
        barsize = 4
      ) +
      
      # EA
      geom_cladelabel(
        node = 288,
        label = "",
        offset = 2.0,
        offset.text = 0.8,
        color = c("#ff6600", "black"),
        barsize = 4
      ) +
      
      # GUFL
      geom_cladelabel(
        node = 259,
        label = "",
        offset = 2.0,
        offset.text = 0.8,
        color = c("#cc00cc", "black"),
        barsize = 4
      ) +
      
      # CH/GUMS
      geom_cladelabel(
        node = 221,
        label = "",
        offset = 2.0,
        offset.text = 0.8,
        color = c("#666666", "black"),
        barsize = 4
      ) +
      
      geom_nodepoint(
        aes(
          label = label,
          subset = as.numeric(sub("/.*", "",
                                  label)) >= 95 &
            !as.numeric(sub(".*/", "",
                            label)) >= 50 & 
            node != east
        ),
        fill = "dodgerblue2",
        size = 2,
        shape = 21,
        na.rm = TRUE
      ) +
      
      geom_nodepoint(
        aes(
          label = label,
          subset = !as.numeric(sub("/.*", "",
                                   label)) >= 95 &
            as.numeric(sub(".*/", "",
                           label)) >= 50
        ),
        fill = "firebrick",
        size = 2,
        shape = 21,
        na.rm = TRUE
      ) +
      
      geom_nodepoint(
        aes(
          label = label,
          subset = as.numeric(sub("/.*", "",
                                  label)) >= 95 &
            as.numeric(sub(".*/", "",
                           label)) >= 50
        ),
        fill = "darkorchid",
        size = 2,
        shape = 21,
        na.rm = TRUE
      ) +
      
      geom_nodepoint(
        aes(label = label,
            subset = node == terrapene),
        fill = "darkorchid",
        size = 2,
        shape = 22,
        na.rm = TRUE
      ) +
      
      geom_nodepoint(
        aes(label = label,
            subset = node == east),
        fill = "dodgerblue2",
        size = 2,
        shape = 22,
        na.rm = TRUE
      ) +
      
      geom_nodepoint(
        aes(label = label,
            subset = node == ONDS),
        fill = "darkorchid",
        size = 2,
        shape = 22,
        na.rm = TRUE
      )
    
  } else {
    
    p.final <- p +
      # outgroup
      geom_cladelabel(
        node = OUT,
        label = "",
        color = c("black", "black"),
        offset = 2.0,
        offset.text = 0.8,
        barsize = 4
      ) +
      
      # DS
      geom_cladelabel(
        node = DS,
        label = "",
        color = c("#08ffff", "black"),
        offset = 2.0,
        offset.text = 0.8,
        barsize = 4
      ) +
      
      # ON
      geom_cladelabel(
        node = ON,
        label = "",
        color = c("#0033cc", "black"),
        offset = 2.0,
        offset.text = 0.8,
        barsize = 4
      ) +
      
      # FL
      geom_cladelabel(
        node = FL,
        label = "",
        color = c("#ffff00", "black"),
        offset = 2.0,
        offset.text = 0.8,
        barsize = 4
      ) +
      
      # MX
      geom_cladelabel(
        node = MX,
        label = "",
        color = c("#ff3333", "black"),
        offset = 2.0,
        offset.text = 0.8,
        barsize = 4
      ) +
      
      # TT
      geom_cladelabel(
        node = TT,
        label = "",
        offset = 2.0,
        offset.text = 0.8,
        color = c("#006600", "black"),
        barsize = 4
      ) +
      
      # EA
      geom_cladelabel(
        node = EA,
        label = "",
        offset = 2.0,
        offset.text = 0.8,
        color = c("#ff6600", "black"),
        barsize = 4
      ) +
      
      # GUFL
      geom_cladelabel(
        node = GUFL,
        label = "",
        offset = 2.0,
        offset.text = 0.8,
        color = c("#cc00cc", "black"),
        barsize = 4
      ) +
      
      # CH/GUMS
      geom_cladelabel(
        node = GUMS_CH,
        label = "",
        offset = 2.0,
        offset.text = 0.8,
        color = c("#666666", "black"),
        barsize = 4
      ) +
      
      geom_nodepoint(
        aes(
          label = label,
          subset = as.numeric(sub("/.*", "",
                                  label)) >= 95 &
            !as.numeric(sub(".*/", "",
                            label)) >= 50 & 
            node != east & 
            node != terrapene & 
            node != ONDS
        ),
        fill = "dodgerblue2",
        size = 2,
        shape = 21,
        na.rm = TRUE
      ) +
      
      geom_nodepoint(
        aes(
          label = label,
          subset = !as.numeric(sub("/.*", "",
                                   label)) >= 95 &
            as.numeric(sub(".*/", "",
                           label)) >= 50 & 
            node != east & 
            node != terrapene & 
            node != ONDS
        ),
        fill = "firebrick",
        size = 2,
        shape = 21,
        na.rm = TRUE
      ) +
      
      geom_nodepoint(
        aes(
          label = label,
          subset = as.numeric(sub("/.*", "",
                                  label)) >= 95 &
            as.numeric(sub(".*/", "",
                           label)) >= 50 & 
            node != east & 
            node != terrapene & 
            node != ONDS
        ),
        fill = "darkorchid",
        size = 2,
        shape = 21,
        na.rm = TRUE
      ) +
      
      geom_nodepoint(
        aes(label = label,
            subset = node == terrapene),
        fill = "darkorchid",
        size = 2,
        shape = 22,
        na.rm = TRUE
      ) +
      
      geom_nodepoint(
        aes(label = label,
            subset = node == east),
        fill = "dodgerblue2",
        size = 2,
        shape = 22,
        na.rm = TRUE
      ) +
      
      geom_nodepoint(
        aes(label = label,
            subset = node == ONDS),
        fill = "darkorchid",
        size = 2,
        shape = 22,
        na.rm = TRUE
      ) +
      
      geom_rootpoint(fill="black", size=2, shape=22, na.rm=TRUE)
    
  }
  
  return(p.final)
}

# Read sampleIDs from file.
# Two columns: sampleID\tPlotting order
sample.ids <- read.table(file = "popmaps/sampleids_numbered.txt", 
                         header = F)

# Create output directory if it doesn't already exist.
dir.create("missDataRuns/plots_maf_final", 
           showWarnings = FALSE)

# Make color vector
colors <- brewer.pal(n = 9, 
                     name = "Set1")

popmap <-
  read.table(
    file = "popmaps/phylogen.popmap.final.txt",
    header = FALSE,
    col.names = c("Sample", "Pop")
  )

RobjectDIR <- "missDataRuns/Robjects_best_maf"

cmds_gapstat <-
  readRDS(
    file = paste0(
      RobjectDIR, 
      "/cmds_gapstat_missInd",
      25,
      "_Pop",
      25,
      "_maf",
      0.05,
      ".rds"
    )
  )
cmds_hier <-
  readRDS(
    file = paste0(
      RobjectDIR, 
      "/cmds_hier_missInd",
      25,
      "_Pop",
      25,
      "_maf",
      0.05,
      ".rds"
    )
  )
cmds_pam <-
  readRDS(
    file = paste0(
      RobjectDIR, 
      "/cmds_pam_missInd",
      25,
      "_Pop",
      25,
      "_maf",
      0.05,
      ".rds"
    )
  )
isomds_gapstat <-
  readRDS(
    file = paste0(
      RobjectDIR, 
      "/isomds_gapstat_missInd",
      25,
      "_Pop",
      25,
      "_maf",
      0.01,
      ".rds"
    )
  )
isomds_hier <-
  readRDS(
    file = paste0(
      RobjectDIR, 
      "/isomds_hier_missInd",
      25,
      "_Pop",
      25,
      "_maf",
      0.01,
      ".rds"
    )
  )
isomds_pam <-
  readRDS(
    file = paste0(
      RobjectDIR, 
      "/isomds_pam_missInd",
      25,
      "_Pop",
      25,
      "_maf",
      0.01,
      ".rds"
    )
  )
tsne_gapstat <-
  readRDS(
    file = paste0(
      RobjectDIR, 
      "/tsne_gapstat_missInd",
      25,
      "_Pop",
      25,
      "_maf",
      0.05,
      "_P",
      15,
      ".rds"
    )
  )
tsne_hier <-
  readRDS(
    file = paste0(
      RobjectDIR, 
      "/tsne_hier_missInd",
      25,
      "_Pop",
      25,
      "_maf",
      0.05,
      "_P",
      15,
      ".rds"
    )
  )
tsne_pam <-
  readRDS(
    file = paste0(
      RobjectDIR, 
      "/tsne_pam_missInd",
      25,
      "_Pop",
      25,
      "_maf",
      0.05,
      "_P",
      15,
      ".rds"
    )
  )
vae <-
  readRDS(
    file = paste0(
      RobjectDIR, 
      "/vae_missInd",
      25,
      "_Pop",
      25,
      "_maf",
      0.05,
      ".rds"
    )
  )

##################################################
## Summarize UML into a barplot
## X-axis = assignment proportion
## Y-axis = samples
##################################################

p1 <-
  sum_plots(cmds_gapstat[[2]], sample.ids, popmap, "cMDS", "GS")

p2 <-
  sum_plots(cmds_hier[[2]], sample.ids, popmap, "cMDS", "H")

p3 <-
  sum_plots(cmds_pam[[2]], sample.ids, popmap, "cMDS", "PAM")

p4 <-
  sum_plots(isomds_gapstat[[2]],
            sample.ids,
            popmap,
            "isoMDS",
            "GS")

p5 <-
  sum_plots(isomds_hier[[2]],
            sample.ids,
            popmap,
            "isoMDS",
            "H")
p6 <-
  sum_plots(isomds_pam[[2]],
            sample.ids,
            popmap,
            "isoMDS",
            "PAM")

p7 <-
  sum_plots(tsne_gapstat[[2]],
            sample.ids,
            popmap,
            paste0("t-SNE P", 15),
            "GS")

p8 <-
  sum_plots(tsne_hier[[2]],
            sample.ids,
            popmap,
            paste0("t-SNE P", 15),
            "H")

p9 <-
  sum_plots(tsne_pam[[2]],
            sample.ids,
            popmap,
            paste0("t-SNE P", 15),
            "PAM")

p10 <-
  sum_plots(vae[[2]],
            sample.ids,
            popmap,
            "VAE",
            "DBSCAN")

###############################################
## Make the tree using ggtree.
###############################################

tfile <- 
  file.path("missDataRuns", "output_dating_part_iqtree",
            "box_dating_full.timetree.nex")

# Read IQ-TREE time-calibrated tree file.
mytree <- treeio::read.beast(file = tfile)

# Make ggtree object from tree..
p.gg <-
  ggtree::ggtree(mytree, right = TRUE) + 
  theme_tree2() + 
  scale_x_continuous(labels = abs) + 
  geom_rootedge(0.4)

# Plot the tree.
p.final <- plot_tree(p.gg)

# Get data from ggtree object.
df.t <- p.final$data

# Drop samples not used in UML analyses.
samples2drop <- subset(p1, is.na(value))
samples2drop <- as.character(unique(samples2drop$Sample))
samples2drop <- append(samples2drop, c("TCAL_BXTC111_U52", 
                                       "GUFL_BXGU38_502", 
                                       "YUYU_BX1183", 
                                       "SSMX_BX1620_x2", 
                                       "TCAL_BX329"))

# drop tips from tree.
mytree2 <- ape::drop.tip(as.phylo(mytree), 
                    samples2drop, 
                    subtree = FALSE, 
                    trim.internal = TRUE)

# Re-make ggtree object from tree with dropped tips.
# Did this again because using drop.tip removes the data
# from the tree object. E.g., divergence times get removed.
# So I made an initial ggtree object, saved the data from it,
# dropped the tips, made a new ggtree object, then re-added the data.
p.gg <-
  ggtree::ggtree(mytree2, right = TRUE)

# Re-add the data to the ggtree object.
p.gg <- p.gg %<+% df.t +
  theme_tree2() + 
  scale_x_continuous(labels = abs) + 
  geom_rootedge(0.4)

# Get most recent common ancestors for specific clades.
OUT <- ape::getMRCA(mytree2, grep("OGCG|KPFP", p.gg$data$label))
DS <- ape::getMRCA(mytree2, grep("DS", p.gg$data$label))
ON <- ape::getMRCA(mytree2, grep("ON", p.gg$data$label))
FL <- ape::getMRCA(mytree2, grep("FLFL", p.gg$data$label))
MX <- ape::getMRCA(mytree2, grep("MXMX", p.gg$data$label))
TT <- ape::getMRCA(mytree2, grep("TT", p.gg$data$label))
EA <- ape::getMRCA(mytree2, grep("EA", p.gg$data$label))
GUFL <-
  ape::getMRCA(
    mytree2,
    grep(
      "GUFL_BX626|GUFL_BXGU35_AA13|GUFL_BXGU33|GUFL_BXGU36_AA14|GUFL_BX503|GUFL_BX684|GUFL_BXGU63_AA37|GUFL_BX504| GUFL_BXGU61_U57|GUFL_BXGU62_AA36|GUFL_BXGU32|GUFL_BXGU65_AA39|GUFL_BX685",
      p.gg$data$label
    )
  )
GUMS_CH <- ape::getMRCA(mytree2, grep("GUMS|CH", p.gg$data$label))

# Fossil calibration nodes.
ONDS <- ape::getMRCA(mytree2, grep("ON|DS", p.gg$data$label))
east <- ape::getMRCA(mytree2, grep("EANC|FLFL|TTTX|MXMX|GUFL|CHCH|GUMS", 
                                   p.gg$data$label))
terrapene <- ape::getMRCA(mytree2, grep("EANC|TTTX|MXMX|FLFL|CHCH|GUFL|GUMS|ONTX|DSNM|", 
                                        p.gg$data$label))

# Plot the tree with again with the labelled clades.
p.final <- plot_tree(p.gg, 
                     clades=TRUE, 
                     OUT, 
                     DS, 
                     ON, 
                     FL, 
                     MX, 
                     TT, 
                     EA, 
                     GUFL, 
                     GUMS_CH, 
                     ONDS=ONDS, 
                     east=east, 
                     terrapene=terrapene)

##############################################
## Plot tree alongside barplots.
##############################################

pa <-
  ggtree::facet_plot(
    p = p.final,
    panel = "cMDS (GS)",
    data = p1,
    geom = ggstance::geom_barh,
    mapping = aes(x = value, fill = variable),
    stat = "identity",
    width = 1
  ) + 
  scale_fill_brewer(palette = "Set1")

pb <-
  ggtree::facet_plot(
    p = pa,
    panel = "cMDS (H)",
    data = p2,
    geom = ggstance::geom_barh,
    mapping = aes(x = value, fill = variable),
    stat = "identity",
    width = 1
  )

pc <-
  ggtree::facet_plot(
    p = pb,
    panel = "cMDS (PAM)",
    data = p3,
    geom = ggstance::geom_barh,
    mapping = aes(x = value, fill = variable),
    stat = "identity",
    width = 1
  )

pd <-
  ggtree::facet_plot(
    p = pc,
    panel = "isoMDS (GS)",
    data = p4,
    geom = ggstance::geom_barh,
    mapping = aes(x = value, fill = variable),
    stat = "identity",
    width = 1
  )

pe <-
  ggtree::facet_plot(
    p = pd,
    panel = "isoMDS (H)",
    data = p5,
    geom = ggstance::geom_barh,
    mapping = aes(x = value, fill = variable),
    stat = "identity",
    width = 1
  )

pf <-
  ggtree::facet_plot(
    p = pe,
    panel = "isoMDS (PAM)",
    data = p6,
    geom = ggstance::geom_barh,
    mapping = aes(x = value, fill = variable),
    stat = "identity",
    width = 1
  )

pg <-
  ggtree::facet_plot(
    p = pf,
    panel = "t-SNE P15 (GS)",
    data = p7,
    geom = ggstance::geom_barh,
    mapping = aes(x = value, fill = variable),
    stat = "identity",
    width = 1
  )

ph <-
  ggtree::facet_plot(
    p = pg,
    panel = "t-SNE P15 (H)",
    data = p8,
    geom = ggstance::geom_barh,
    mapping = aes(x = value, fill = variable),
    stat = "identity",
    width = 1
  )

pi <-
  ggtree::facet_plot(
    p = ph,
    panel = "t-SNE P15 (PAM)",
    data = p9,
    geom = ggstance::geom_barh,
    mapping = aes(x = value, fill = variable),
    stat = "identity",
    width = 1
  )

pj <-
  ggtree::facet_plot(
    p = pi,
    panel = "VAE (DBSCAN)",
    data = p10,
    geom = ggstance::geom_barh,
    mapping = aes(x = value, fill = variable),
    stat = "identity",
    width = 1
  ) + 
  theme(legend.position = "none")

# Change tree panel width (makes it 3 times larger).
gt <- ggplot_gtable(ggplot_build(pj))
gt$layout$l[grep('panel-1', gt$layout$name)]
gt$widths[5] = 3*gt$widths[5]
grid::grid.draw(gt)

# Save the final plot.
ggtree::ggsave(filename = file.path("missDataRuns",
                                    "plots_maf_final2_best",
                                    "missInd25_Pop25_best_withTree.pdf"), 
       plot = gt, 
       device="pdf", 
       width=16, 
       height = 9, 
       units = "in", 
       dpi=300)
