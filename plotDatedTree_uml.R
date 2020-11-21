
# Plot iqtree phylogeny with divergence times.
# For Terrapene

library("ggplot2")
library("ggtree")
library("treeio")
library("gridExtra")
library("ape")
library("RColorBrewer")
library("aplot")
library("reshape2")

sum_plots <- function(propsdf, sample.ids, pops, title, subtitle){
  
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

sample.ids <- read.table(file = "popmaps/sampleids_numbered.txt", header = F)

dir.create("missDataRuns/plots_maf_final", showWarnings = FALSE)

# Make color vector
colors <- brewer.pal(n = 9, name = "Set1")

popmap <-
  read.table(
    file = "popmaps/phylogen.popmap.final.txt",
    header = FALSE,
    col.names = c("Sample", "Pop")
  )

cmds_gapstat <-
  readRDS(
    file = paste0(
      "missDataRuns/Robjects_maf/cmds_gapstat_missInd",
      50,
      "_Pop",
      50,
      "_maf",
      0.03,
      ".rds"
    )
  )
cmds_hier <-
  readRDS(
    file = paste0(
      "missDataRuns/Robjects_maf/cmds_hier_missInd",
      50,
      "_Pop",
      50,
      "_maf",
      0.03,
      ".rds"
    )
  )
cmds_pam <-
  readRDS(
    file = paste0(
      "missDataRuns/Robjects_maf/cmds_pam_missInd",
      50,
      "_Pop",
      50,
      "_maf",
      0.03,
      ".rds"
    )
  )
cmds_pam_prox <-
  readRDS(
    file = paste0(
      "missDataRuns/Robjects_maf/cmds_pam_prox_missInd",
      50,
      "_Pop",
      50,
      "_maf",
      0.03,
      ".rds"
    )
  )
isomds_gapstat <-
  readRDS(
    file = paste0(
      "missDataRuns/Robjects_maf/isomds_gapstat_missInd",
      50,
      "_Pop",
      50,
      "_maf",
      0.03,
      ".rds"
    )
  )
isomds_hier <-
  readRDS(
    file = paste0(
      "missDataRuns/Robjects_maf/isomds_hier_missInd",
      50,
      "_Pop",
      50,
      "_maf",
      0.03,
      ".rds"
    )
  )
isomds_pam <-
  readRDS(
    file = paste0(
      "missDataRuns/Robjects_maf/isomds_pam_missInd",
      50,
      "_Pop",
      50,
      "_maf",
      0.03,
      ".rds"
    )
  )


#*********************************************************
### D) Set function values here
#*********************************************************
cmds.gapstat.prefix <-
  paste0("cmds_gapstat_missInd", 50, "_Pop", 50, "_maf", 0.03)
cmds.gapstat.title <-
  paste0("cMDS Gap Statistic missInd ", 50, " Pop ", 50, "_maf", 0.03)

# cMDS Hierarchical Clustering
cmds.hier.prefix <-
  paste0("cmds_hierarchical_missInd", 50, "_Pop", 50, "_maf", 0.03)
cmds.hier.title <-
  paste0("cMDS Hierarchical Clustering missInd ",
         50,
         " Pop ",
         50,
         "_maf",
         0.03)

# cMDS PAM Clustering (cMDS Groups)
cmds.pam.cmdsgroups.prefix <-
  paste0("cmds_pam_missInd", 50, "_Pop", 50, "_maf", 0.03)
cmds.pam.cmdsgroups.title <-
  paste0("cMDS PAM Clustering (cMDS Groups) missInd ",
         50,
         " Pop ",
         50,
         "_maf",
         0.03)

# cMDS PAM Clustering (Proximity Scores)
cmds.pam.prox.prefix <-
  paste0("cmds_pam_prox_missInd", 50, "_Pop", 50, "_maf", 0.03)
cmds.pam.prox.title <-
  paste0("cMDS PAM Clustering (Proximity Scores) missInd ",
         50,
         " Pop ",
         50,
         "_maf",
         0.03)

# isoMDS Gap Statistic (PAM Clustering)
isomds.gapstat.prefix <-
  paste0("isomds_gapstat_missInd", 50, "_Pop", 50, "_maf", 0.03)
isomds.gapstat.title <-
  paste0("isoMDS Gap Statistic missInd ", 50, " Pop ", 50, "_maf", 0.03)

# isoMDS Hierarchical Clustering
isomds.hier.prefix <-
  paste0("isomds_hierarchical_missInd", 50, "_Pop", 50, "_maf", 0.03)
isomds.hier.title <-
  paste0("isoMDS Hierarchical Clustering missInd ",
         50,
         " Pop ",
         50,
         "_maf",
         0.03)

# isoMDS PAM Clustering
isomds.pam.prefix <-
  paste0("isomds_pam_missInd", 50, "_Pop", 50, "_maf", 0.03)
isomds.pam.title <-
  paste0("isoMDS PAM Clustering missInd ", 50, " Pop ", 50, "_maf", 0.03)

# Make all the plots
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

p1 <-
  sum_plots(cmds_gapstat[[2]], sample.ids, popmap, "cMDS", "GS")

p2 <-
  sum_plots(cmds_hier[[2]], sample.ids, popmap, "cMDS", "H")

p3 <-
  sum_plots(cmds_pam[[2]], sample.ids, popmap, "cMDS", "PAM")

p4 <-
  sum_plots(cmds_pam_prox[[2]],
            sample.ids,
            popmap,
            "Prox",
            "PAM")

p5 <-
  sum_plots(isomds_gapstat[[2]],
            sample.ids,
            popmap,
            "isoMDS",
            "GS")

p6 <-
  sum_plots(isomds_hier[[2]],
            sample.ids,
            popmap,
            "isoMDS",
            "H")
p7 <-
  sum_plots(isomds_pam[[2]],
            sample.ids,
            popmap,
            "isoMDS",
            "PAM")

writeLines(paste0("Doing t-SNE perplexity ", 15, "...\n"))
tsne_gapstat <-
  readRDS(
    file = paste0(
      "missDataRuns/Robjects_maf/tsne_gapstat_missInd",
      50,
      "_Pop",
      50,
      "_maf",
      0.03,
      "_P",
      15,
      ".rds"
    )
  )
tsne_hier <-
  readRDS(
    file = paste0(
      "missDataRuns/Robjects_maf/tsne_hier_missInd",
      50,
      "_Pop",
      50,
      "_maf",
      0.03,
      "_P",
      15,
      ".rds"
    )
  )
tsne_pam <-
  readRDS(
    file = paste0(
      "missDataRuns/Robjects_maf/tsne_pam_missInd",
      50,
      "_Pop",
      50,
      "_maf",
      0.03,
      "_P",
      15,
      ".rds"
    )
  )

# tsne perplexity = 5 Gap Statistic
tsnep5.gapstat.prefix <-
  paste0("tsneP50_gapstat_missInd",
         50,
         "_Pop",
         50,
         "_maf",
         0.03,
         "_P",
         15)
tsnep5.gapstat.title <-
  paste0("t-SNE Gap Statistic missInd ",
         50,
         " Pop",
         50,
         "_maf",
         0.03,
         "(Perplexity = ",
         15,
         ")")

# tsne perplexity = 5 Hierarchical Clustering
tsnep5.hier.prefix <-
  paste0("tsneP50_hierarchical_missInd",
         50,
         "_Pop",
         50,
         "_maf",
         0.03,
         "_P",
         15)

tsnep5.hier.title <-
  paste0(
    "t-SNE Hierarchical Clustering missInd ",
    50,
    " Pop",
    50,
    "_maf",
    0.03,
    "(Perplexity = ",
    15,
    ")"
  )

# tsne perplexity = 5 PAM Clustering
tsnep5.pam.prefix <-
  paste0("tsneP50_pam_missInd", 50, "_Pop", 50, "_maf", 0.03)

tsnep5.pam.title <-
  paste0("t-SNE PAM Clustering missInd ",
         50,
         " Pop",
         50,
         "_maf",
         0.03,
         "(Perplexity = ",
         15,
         ")")


#*********************************************************

##########################################################
### E) Run the functions for each uml method
##########################################################

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

###########################################################
### F) Summarize Runs (Rows)
###########################################################

p8 <-
  sum_plots(tsne_gapstat[[2]],
            sample.ids,
            popmap,
            paste0("t-SNE P", 15),
            "GS")

p9 <-
  sum_plots(tsne_hier[[2]],
            sample.ids,
            popmap,
            paste0("t-SNE P", 15),
            "H")

p10 <-
  sum_plots(tsne_pam[[2]],
            sample.ids,
            popmap,
            paste0("t-SNE P", 15),
            "PAM")

plot_tree <- function(p, clades=NULL, OUT=NULL, DS=NULL, ON=NULL, FL=NULL, MX=NULL, TT=NULL, EA=NULL, GUFL=NULL, GUMS_CH=NULL, ONDS=NULL, terrapene=NULL, east=NULL) {
  
  p <- revts(p)
  
  if (is.null(clades)) {

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
                            label)) >= 50 & node != east
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
                            label)) >= 50 & node != east & node != terrapene & node != ONDS
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
                           label)) >= 50 & node != east & node != terrapene & node != ONDS
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
                           label)) >= 50 & node != east & node != terrapene & node != ONDS
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

tfile <- 
  file.path("output_dating_part_iqtree",
            "box_dating_full.timetree.nex")

mytree <- read.beast(file = tfile)

p.gg <-
  ggtree(mytree, right = TRUE) + 
  theme_tree2() + 
  scale_x_continuous(labels = abs) + 
  geom_rootedge(0.4)

p.final <- plot_tree(p.gg)

df.t <- p.final$data

samples2drop <- subset(p1, is.na(value))
samples2drop <- as.character(unique(samples2drop$Sample))
samples2drop <- append(samples2drop, c("TCAL_BXTC111_U52", "GUFL_BXGU38_502", "YUYU_BX1183", "SSMX_BX1620_x2", "TCAL_BX329"))

mytree2 <- drop.tip(as.phylo(mytree), samples2drop, subtree = FALSE, trim.internal = TRUE)

p.gg <-
  ggtree(mytree2, right = TRUE)

p.gg <- p.gg %<+% df.t +
  theme_tree2() + 
  scale_x_continuous(labels = abs) + 
  geom_rootedge(0.4)

p.final <- plot_tree(p.gg)

OUT <- getMRCA(mytree2, grep("OGCG|KPFP", p.final$data$label))
DS <- getMRCA(mytree2, grep("DS", p.final$data$label))
ON <- getMRCA(mytree2, grep("ON", p.final$data$label))
FL <- getMRCA(mytree2, grep("FLFL", p.final$data$label))
MX <- getMRCA(mytree2, grep("MXMX", p.final$data$label))
TT <- getMRCA(mytree2, grep("TT", p.final$data$label))
EA <- getMRCA(mytree2, grep("EA", p.final$data$label))
GUFL <-
  getMRCA(
    mytree2,
    grep(
      "GUFL_BXGU27|GUFL_BXGU65_AA39|GUFL_BX504|GUFL_BX503|GUFL_BX684|GUFL_BXGU32|GUFL_BX626|GUFL_BX627",
      p.final$data$label
    )
  )

GUMS_CH <- getMRCA(mytree2, grep("GUMS|CH", p.final$data$label))

ONDS <- getMRCA(mytree2, grep("ON|DS", p.final$data$label))
east <- getMRCA(mytree2, c("EANC_BX316", "FLFL_BX683"))
terrapene <- getMRCA(mytree2, c("EANC_BX316", "ONTX_BX765"))

p.final <- plot_tree(p.gg, clades="use", OUT, DS, ON, FL, MX, TT, EA, GUFL, GUMS_CH, ONDS=ONDS, east=east, terrapene=terrapene)

pa <- facet_plot(p = p.final, panel = "cMDS (GS)", data = p1, geom = ggstance::geom_barh, mapping = aes(x = value, fill = variable), stat = "identity", width = 1) + scale_fill_brewer(palette = "Set1")

pb <- facet_plot(p = pa, panel = "cMDS (H)", data = p2, geom=ggstance::geom_barh, mapping = aes(x = value, fill = variable), stat = "identity", width = 1)

pc <- facet_plot(p = pb, panel = "cMDS (PAM)", data = p3, geom=ggstance::geom_barh, mapping = aes(x = value, fill = variable), stat = "identity", width = 1)

pd <- facet_plot(p = pc, panel = "cMDS (PAM Prox)", data = p4, geom=ggstance::geom_barh, mapping = aes(x = value, fill = variable), stat = "identity", width = 1)

pe <- facet_plot(p = pd, panel = "isoMDS (GS)", data = p5, geom=ggstance::geom_barh, mapping = aes(x = value, fill = variable), stat = "identity", width = 1)

pf <- facet_plot(p = pe, panel = "isoMDS (H)", data = p6, geom=ggstance::geom_barh, mapping = aes(x = value, fill = variable), stat = "identity", width = 1)

pg <- facet_plot(p = pf, panel = "isoMDS (PAM)", data = p7, geom=ggstance::geom_barh, mapping = aes(x = value, fill = variable), stat = "identity", width = 1)

ph <- facet_plot(p = pg, panel = "t-SNE P15 (GS)", data = p8, geom=ggstance::geom_barh, mapping = aes(x = value, fill = variable), stat = "identity", width = 1)

pi <- facet_plot(p = ph, panel = "t-SNE P15 (H)", data = p9, geom=ggstance::geom_barh, mapping = aes(x = value, fill = variable), stat = "identity", width = 1)

pj <- facet_plot(p = pi, panel = "t-SNE P15 (PAM)", data = p10, geom=ggstance::geom_barh, mapping = aes(x = value, fill = variable), stat = "identity", width = 1) + theme(legend.position = "none")

gt <- ggplot_gtable(ggplot_build(pj))
gt$layout$l[grep('panel-1', gt$layout$name)]
gt$widths[5] = 3*gt$widths[5]
grid::grid.draw(gt)

ggsave(filename = "missInd50_Pop50_maf0.03.pdf", plot = gt, device="pdf", width=16, height = 9, units = "in", dpi=300)