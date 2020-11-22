library("ggplot2")
library("dplyr")
library("tibble")
library("reshape2")
library("RColorBrewer")
library("gridExtra")
library("ggpubr")
library("wesanderson")

plotUMLsummary <- function(runs, ind, pop, K){
  p <- ggplot(data = runs, aes(x = ind, y = pop, color = K)) + geom_line()
}

pad_columns <- function(df){
  
  if (ncol(df) < 101){
    while(ncol(df) < 101){
      pad <- ncol(df)
      df <- df %>%
        mutate("X{pad}" := NA)
    }
    return(df)
  } else{
    return(df)
  }
}

pad_rows <- function(df, max_inds){
  
  if (nrow(df) < max_inds){
    while(nrow(df) < max_inds){
      pad <- nrow(df)
      df <- df %>%
        mutate("X{pad}" := NA)
    }
    return(df)
  } else{
    return(df)
  }
}

my_theme <- function(...){
  list(
    geom_tile(color = "white") ,
    theme_bw() +
      theme(
        axis.title = element_text(colour = "black", size = 22),
        axis.text = element_text(colour = "black", size = 22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        strip.text = element_text(colour = "black", size = 22),
        legend.title = element_text(size=22, colour="black"),
        legend.text = element_text(size=22, colour="black"),
        ...
      )
  )
}

maf_vals <- c("0.0", "0.01", "0.03", "0.05")

uml_list <- list()
avgK <- list()
counter <- 1
for (i in seq(25, 100, 25)) {
  for (j in seq(25, 100, 25)) {
    for (m in 1:length(maf_vals)) {
      
      if (i == 100 & j == 25){
        next
      }
      
      writeLines(paste0("\nDoing missInd ", i, " Pop ", j, "_maf", maf_vals[m], "..."))
      # Load aligned objects from UML_alignK.R
      cmds_gapstat <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/cmds_gapstat_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            maf_vals[m],
            ".rds"
          )
        )
      
      cmds_gapstat <- cmds_gapstat[[1]]
      
      cmds_gapstat <- pad_columns(cmds_gapstat)
      
      cmds_gapstat$method <- "cmds"
      cmds_gapstat$clust <- "gs"
      cmds_gapstat$perp <- NA
      cmds_gapstat$ind <- i
      cmds_gapstat$pop <- j
      cmds_gapstat$maf <- maf_vals[m]
      
      uml_list[[counter]] <- cmds_gapstat
      
      counter <- counter + 1
      
      cmds_hier <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/cmds_hier_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            maf_vals[m],
            ".rds"
          )
        )
      
      cmds_hier <- cmds_hier[[1]]
      
      cmds_hier <- pad_columns(cmds_hier)
      
      cmds_hier$method <- "cmds"
      cmds_hier$clust <- "hier"
      cmds_hier$perp <- NA
      cmds_hier$ind <- i
      cmds_hier$pop <- j
      cmds_hier$maf <- maf_vals[m]
      
      uml_list[[counter]] <- cmds_hier
      counter <- counter + 1
      
      cmds_pam <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/cmds_pam_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            maf_vals[m],
            ".rds"
          )
        )
      
      cmds_pam <- cmds_pam[[1]]
      
      cmds_pam <- pad_columns(cmds_pam)
      
      cmds_pam$method <- "cmds"
      cmds_pam$clust <- "pam"
      cmds_pam$perp <- NA
      cmds_pam$ind <- i
      cmds_pam$pop <- j
      cmds_pam$maf <- maf_vals[m]
      
      uml_list[[counter]] <- cmds_pam
      
      counter <- counter + 1
      
      isomds_gapstat <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/isomds_gapstat_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            maf_vals[m],
            ".rds"
          )
        )
      
      isomds_gapstat <- isomds_gapstat[[1]]
      
      isomds_gapstat <- pad_columns(isomds_gapstat)
      
      isomds_gapstat$method <- "isomds"
      isomds_gapstat$clust <- "gs"
      isomds_gapstat$perp <- NA
      isomds_gapstat$ind <- i
      isomds_gapstat$pop <- j
      isomds_gapstat$maf <- maf_vals[m]
      
      uml_list[[counter]] <- isomds_gapstat
      counter <- counter + 1
      
      isomds_hier <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/isomds_hier_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            maf_vals[m],
            ".rds"
          )
        )
      
      isomds_hier <- isomds_hier[[1]]
      
      isomds_hier <- pad_columns(isomds_hier)
      
      isomds_hier$method <- "isomds"
      isomds_hier$clust <- "hier"
      isomds_hier$perp <- NA
      isomds_hier$ind <- i
      isomds_hier$pop <- j
      isomds_hier$maf <- maf_vals[m]
      
      uml_list[[counter]] <- isomds_hier
      counter <- counter + 1
      
      isomds_pam <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/isomds_pam_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            maf_vals[m],
            ".rds"
          )
        )
      
      isomds_pam <- isomds_pam[[1]]
      
      isomds_pam <- pad_columns(isomds_pam)
      
      isomds_pam$method <- "isomds"
      isomds_pam$clust <- "pam"
      isomds_pam$perp <- NA
      isomds_pam$ind <- i
      isomds_pam$pop <- j
      isomds_pam$maf <- maf_vals[m]
      
      uml_list[[counter]] <- isomds_pam
      counter <- counter + 1
      
      vae <-
        readRDS(
          file = paste0(
            "missDataRuns/Robjects_maf/vae_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            maf_vals[m],
            ".rds"
          )
        )
      
      vae <- vae[[1]]
      
      vae <- pad_columns(vae)
      
      vae$method <- "vae"
      vae$clust <- "DBSCAN"
      vae$perp <- NA
      vae$ind <- i
      vae$pop <- j
      vae$maf <- maf_vals[m]
      
      uml_list[[counter]] <- vae
      counter <- counter + 1
      
      for (p in seq(5, 50, 5)) {
        tsne_gapstat <-
          readRDS(
            file = paste0(
              "missDataRuns/Robjects_maf/tsne_gapstat_missInd",
              i,
              "_Pop",
              j,
              "_maf",
              maf_vals[m],
              "_P",
              p,
              ".rds"
            )
          )
        
        tsne_gapstat <- tsne_gapstat[[1]]
        
        tsne_gapstat <- pad_columns(tsne_gapstat)
        
        tsne_gapstat$method <- "tsne"
        tsne_gapstat$clust <- "gs"
        tsne_gapstat$perp <- p
        tsne_gapstat$ind <- i
        tsne_gapstat$pop <- j
        tsne_gapstat$maf <- maf_vals[m]
        
        uml_list[[counter]] <- tsne_gapstat
        counter <- counter + 1
        
        tsne_hier <-
          readRDS(
            file = paste0(
              "missDataRuns/Robjects_maf/tsne_hier_missInd",
              i,
              "_Pop",
              j,
              "_maf",
              maf_vals[m],
              "_P",
              p,
              ".rds"
            )
          )
        
        tsne_hier <- tsne_hier[[1]]
        
        tsne_hier <- pad_columns(tsne_hier)
        
        tsne_hier$method <- "tsne"
        tsne_hier$clust <- "hier"
        tsne_hier$perp <- p
        tsne_hier$ind <- i
        tsne_hier$pop <- j
        tsne_hier$maf <- maf_vals[m]
        
        uml_list[[counter]] <- tsne_hier
        counter <- counter + 1
        
        tsne_pam <-
          readRDS(
            file = paste0(
              "missDataRuns/Robjects_maf/tsne_pam_missInd",
              i,
              "_Pop",
              j,
              "_maf",
              maf_vals[m],
              "_P",
              p,
              ".rds"
            )
          )
        
        tsne_pam <- tsne_pam[[1]]
        
        tsne_pam <- pad_columns(tsne_pam)
        
        tsne_pam$method <- "tsne"
        tsne_pam$clust <- "pam"
        tsne_pam$perp <- p
        tsne_pam$ind <- i
        tsne_pam$pop <- j
        tsne_pam$maf <- maf_vals[m]
        
        uml_list[[counter]] <- tsne_pam
        counter <- counter + 1
      }
    }
  }
}

uml_list_working <- lapply(uml_list, function(x) {
  cbind(lapply(x[,c(2:101)], as.integer), x[,c(102:ncol(x))])
})

# Get optimal K for each run in each data.frame of a list of data.frames.
maxK <- lapply(uml_list_working, function(x){
  cbind(lapply(x[,c(1:100)], function(y) length(unique(y))), x[, c(101:ncol(x))])
})

maxK_2 <- lapply(maxK, function(x){
  x$runid <- paste0(x$ind, "_", x$pop)
  x$runid_perp <- paste0(x$ind, "_", x$pop, "_", x$perp)
  x$runid_maf <- paste0(x$ind, "_", x$pop, "_", x$maf)
  x$runid_maf_perp <- paste0(x$ind, "_", x$pop, "_", x$maf, "_", x$perp)
  return(x)
})

max_inds <- max(sapply(maxK_2, function(x) nrow(x)))

combined <- maxK_2 %>% bind_rows(.id="runid")

mymelt <-
  melt(
    combined,
    na.rm = TRUE,
    id.vars = c(
      "method",
      "clust",
      "perp",
      "ind",
      "pop",
      "maf",
      "runid",
      "runid_perp",
      "runid_maf",
      "runid_maf_perp",
      "X101"
    ),
    variable.name = "runNUM",
    value.name = "K"
  )

se <- function(x) {
  qt(0.975,df=99)*sd(x)/sqrt(100)
}

hm_data <- mymelt %>% group_by(runid, method, clust, ind, pop, maf, perp) %>% summarise(meanK=mean(K), n=n(), sd=sd(K), var=var(K), seK=qt(0.975,df=99)*sd/sqrt(100),ci_low = mean(K)-qt(0.975, df=99)*seK, ci_high=mean(K)+qt(0.975, df=99)*seK)

hm_data$meanK <- ifelse(hm_data$meanK < 2, round(hm_data$meanK), hm_data$meanK)

hm_tsne <- subset(hm_data, hm_data$method == "tsne")

################################################
## Plot it!
################################################

clust_names <- as_labeller(c(`cmds` = "cMDS", `isomds` = "isoMDS", `tsne` = "t-SNE", `gs` = "GS", `hier` = "HC", `pam` = "PAM", `0.0` = "0%", `0.01` = "1%", `0.03` = "3%", `0.05` = "5%"))

clust_names_vae <- as_labeller(c(`0.0` = "0%", `0.01` = "1%", `0.03` = "3%", `0.05` = "5%"))

perp_names <- as_labeller(c(`5` = "P5", `10` = "P10", `15` = "P15", `20` = "P20", `25` = "P25", `30` = "P30", `35` = "P35", `40` = "P40", `45` = "P45", `50` = "P50", `gs` = "GS", `hier` = "HC", `pam` = "PAM", `0.0` = "0%", `0.01` = "1%", `0.03` = "3%", `0.05` = "5%"))

psd.avg <-
  ggplot(data = hm_data, aes(
    x = factor(ind),
    y = factor(pop),
    fill = sd
  )) + my_theme() + scale_fill_gradientn(
    name = "SD",
    limits = c(0, 4),
    breaks = c(0, 1, 2, 3, 4),
    colours = c("blue", "#FDFD96", "red")
  ) +
  scale_y_discrete(
    labels = c("25", "", "75", ""),
    breaks = c("25", "", "75", "")
  ) +
  xlab("Per-individual Filtering (%)") +
  ylab("Per-population Filtering (%)") + 
  facet_grid(clust+factor(maf)~method, labeller = clust_names)

pK.avg <-
  ggplot(data = hm_data, aes(
    x = factor(ind),
    y = factor(pop),
    fill = meanK
  )) + my_theme() + scale_fill_gradientn(
    name = "Optimal K",
    limits = c(2, 9),
    breaks = c(2, 4, 6, 8),
    colours = c("blue", "#FDFD96", "red")
  ) +
  scale_y_discrete(
    labels = c("25", "", "75", ""),
    breaks = c("25", "", "75", "")
  ) +
  xlab("Per-individual Filtering (%)") +
  ylab("Per-population Filtering (%)") + 
  facet_grid(clust+factor(maf)~method, labeller = clust_names)

psd.avg.vae <-
  ggplot(data = subset(hm_data, method == "vae"), 
  aes(
    x = factor(ind),
    y = factor(pop),
    fill = sd
  )) + 
  my_theme() + 
  scale_fill_gradientn(
    name = "SD",
    limits = c(0, 4),
    breaks = c(0, 1, 2, 3, 4),
    colours = c("blue", "#FDFD96", "red")
  ) +
  scale_y_discrete(
    labels = c("25", "50", "75", "100"),
    breaks = c("25", "50", "75", "100")
  ) +
  xlab("Per-individual Filtering (%)") +
  ylab("Per-population Filtering (%)") + 
  facet_wrap(~factor(maf), 
             labeller = clust_names_vae, 
             nrow = 4, 
             ncol = 1) +
  coord_fixed()

pK.avg.vae <-
  ggplot(data = subset(hm_data, method == "vae"), 
    aes(
      x = factor(ind),
      y = factor(pop),
      fill = meanK
  )) + 
  my_theme() + 
  scale_fill_gradientn(
    name = "Optimal K",
    limits = c(2, 10),
    breaks = c(2, 4, 6, 8, 10),
    colours = c("blue", "#FDFD96", "red")
  ) +
  scale_y_discrete(
    labels = c("25", "50", "75", "100"),
    breaks = c("25", "50", "75", "100")
  ) +
  xlab("Per-individual Filtering (%)") +
  ylab("Per-population Filtering (%)") + 
  facet_wrap(~factor(maf), 
             labeller = clust_names_vae, 
             nrow = 4, 
             ncol = 1) +
  coord_fixed()


#####################################################
## t-SNE perplexities
#####################################################

tsne.perp.avgK <-
  ggplot(data = hm_tsne, aes(
    x = factor(ind),
    y = factor(pop),
    fill = meanK
  )) + my_theme() +
  scale_fill_gradientn(
    name = "Optimal K",
    limits = c(1, 9),
    breaks = c(2,4,6,8),
    colours = c("blue", "#FDFD96", "red")
  ) +
  scale_y_discrete(
    labels = c("25", "", "75", ""),
    breaks = c("25", "", "75", "")
    ) +
  xlab("Per-individual Filtering (%)") +
  ylab("Per-population Filtering (%)") + 
  facet_grid(perp~clust+maf, labeller = perp_names)

tsne.perp.sd <-
  ggplot(data = hm_tsne, aes(
    x = factor(ind),
    y = factor(pop),
    fill = sd
  )) + my_theme() +
  scale_fill_gradientn(
    name = "SD",
    limits = c(0, 4),
    breaks = c(0,1,2,3,4),
    colours = c("blue", "#FDFD96", "red")
  ) + 
  scale_y_discrete(
    labels = c("25", "", "75", ""),
    breaks = c("25", "", "75", "")
  ) +
  xlab("Per-individual Filtering (%)") +
  ylab("Per-population Filtering (%)") + 
  facet_grid(perp~clust+maf, labeller = perp_names)

ggsave(filename = "missDataRuns/plots_maf_final/avgK_mean.pdf", plot = pK.avg, device = "pdf", width = 10, height = 10, units = "in", dpi = 300)

ggsave(filename = "missDataRuns/plots_maf_final/avgK_sd.pdf", plot = psd.avg, device = "pdf", width = 10, height = 10, units = "in", dpi = 300)

ggsave(filename = "missDataRuns/plots_maf_final2/avgK_vae_sd.pdf", plot = psd.avg.vae, device = "pdf", width = 10, height = 10, units = "in", dpi = 300)

ggsave(filename = "missDataRuns/plots_maf_final2/avgK_vae_mean.pdf", plot = pK.avg.vae, device = "pdf", width = 10, height = 10, units = "in", dpi = 300)

ggsave(filename = "missDataRuns/plots_maf_final2/avgK_tsne_mean.pdf", plot = tsne.perp.avgK, device = "pdf", width = 11, height = 11, units = "in", dpi = 300)

ggsave(filename = "missDataRuns/plots_maf_final2/avgK_tsne_sd.pdf", plot = tsne.perp.sd, device = "pdf", width = 11, height = 11, units = "in", dpi = 300)

###################################################
## Line plots
###################################################

my_theme2 <- function(...){
  list(
    theme_bw() +
      theme(
        axis.title = element_text(colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 14),
        panel.grid = element_blank(), 
        strip.text = element_text(size=14, color = "black"),
        legend.text = element_text(size=14, color = "black"),
        legend.title = element_text(size=14, color = "black"),
        ...
      ) ,
    scale_fill_manual(values=wes_palette(n=3, name="Royal1"), name="Clustering Algorithm", labels = c("GS", "HC", "PAM"))
  )
}

my_theme3 <- function(...){
  list(
    theme_bw() +
      theme(
        axis.title = element_text(colour = "black", size = 22),
        axis.text = element_text(colour = "black", size = 22),
        panel.grid = element_blank(), 
        strip.text = element_text(size=22, color = "black"),
        legend.text = element_text(size=14, color = "black"),
        legend.title = element_text(size=14, color = "black"),
        ...
      )
  )
}

my_theme4 <- function(...){
  list(
    theme_bw() +
      theme(
        axis.title = element_text(colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 14),
        panel.grid = element_blank(), 
        strip.text = element_text(size=14, color = "black"),
        legend.text = element_text(size=14, color = "black"),
        legend.title = element_text(size=14, color = "black"),
        ...
      ) ,
    scale_fill_brewer(
      palette = "Set3",
      name = "Clustering Alg. (MAF)",
      labels = c(
        "GS (0%)",
        "H (0%)",
        "PAM (0%)",
        "GS (1%)",
        "H (1%)",
        "PAM (1%)",
        "GS (3%)",
        "H (3%)",
        "PAM (3%)",
        "GS (5%)",
        "H (5%)",
        "PAM (5%)"
      )
    )
  )
}

plot_grpbox <- function(df){
  
  p.maf1 <-
    ggplot(subset(df, method != "tsne"), 
           aes(x = method, y = K, 
               fill = interaction(clust, maf), 
               group = runid)) +
    geom_boxplot(
      outlier.size = NA,
      outlier.colour = NA,
      outlier.fill = NA,
      outlier.stroke = NA,
      outlier.alpha = NA,
      outlier.shape = NA,
      notch = FALSE
    ) +
    facet_grid(forcats::fct_rev(factor(pop)) ~ factor(ind)) +
    my_theme4() +
    xlab("Dimension Reduction Method") +
    ylab(expression(italic("K"))) +
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10))
  
  p.maf2 <-
    ggplot(subset(df, method == "tsne"),
           aes(x = factor(perp),
               y = K, 
               fill = factor(maf))) +
    geom_boxplot(
      outlier.size = NA,
      outlier.colour = NA,
      outlier.fill = NA,
      outlier.stroke = NA,
      outlier.alpha = NA,
      outlier.shape = NA,
      notch = FALSE
    ) +
    facet_grid(forcats::fct_rev(factor(pop)) ~ factor(ind)) +
    my_theme3() +
    xlab("t-SNE Perplexity") +
    ylab(expression(italic("K"))) +
    scale_y_continuous(breaks = c(2, 4, 6, 8, 10)) +
    scale_fill_brewer(palette = "Set3", name = "MAF Filter (%)") + 
    scale_x_discrete(labels = c("", "10", "", "20", "", "30", "", "40", "", "50"))
  
  gc()
  return(list(p.maf1, p.maf2))
}

myp <- plot_grpbox(mymelt)

ggsave(filename = "missDataRuns/plots_maf_final/boxplots_maf.png", myp[[1]], device="png", width = 11, height = 7, units = "in", dpi=300)

ggsave(filename = "missDataRuns/plots_maf_final/boxplots_perp_maf.png", myp[[2]], device="png", width = 16, height = 9, units = "in", dpi=300)

ggsave(filename = "missDataRuns/plots_maf_final/boxplots_maf.pdf", myp[[1]], device="pdf", width = 11, height = 7, units = "in", dpi=300)

ggsave(filename = "missDataRuns/plots_maf_final/boxplots_perp_maf.pdf", myp[[2]], device="pdf", width = 16, height = 9, units = "in", dpi=300)




