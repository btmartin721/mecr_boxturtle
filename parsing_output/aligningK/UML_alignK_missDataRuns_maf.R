library(pophelper)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(combinat)
library(stringr)

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

runs2table <- function(dir, sampleNamesFromFile=F){
  file.names <- dir(dir)
  first=TRUE
  for (f in seq(1,length(file.names))){
    fpath<-file.path(dir, file.names[[f]])
    fdf<-read.table(fpath, header=F, sep="\t", stringsAsFactors=FALSE)
    if (first==TRUE){
      if (sampleNamesFromFile==TRUE){
        t <- data.frame("Sample"=(fdf[1]))
      }else{
        t<-data.frame(Sample=(seq(1,nrow(fdf))))
      }
      colnames(t) <- c("Sample")
      first=FALSE
    }
    t[paste0("X",f)] <- fdf[,-1]
  }
  return(t)
}

tab2Qlist <- function(x){
  ret <- list()
  count=1
  for (i in seq(2, ncol(x))){
    sub <- x[,c(1,i)]
    qdf <- data.frame(Sample=sub$Sample)
    for (k in sort(unique(sub[,2]))){
      #if (k == 0){next}
      qdf[as.character(k)] <- 0.0
      qdf[is.element(qdf$Sample, sub[sub[2]==k,"Sample"]),as.character(k)] <- 1.0
    }
    ret[[count]] <- as.data.frame(qdf[,-c(1)])
    count <- count+1
  }
  return(ret)
}

#pass a Q list (from alignK in pophelper)
alignQ2tab <- function(x, sampleNames, zeroIndexing=F){
  df <- data.frame(Sample=sampleNames)
  for (run in seq(1, length(x))){
    sub<-x[[run]]
    run_label <- paste0("X", as.character(run))
    df[run_label] <- NA
    for (c in seq(1, ncol(sub))){
      label <- as.character(c)
      if (zeroIndexing == TRUE){label <- as.character(c-1)}
      df[sub[c]==1,run_label]<-label
    }
  }
  return(df)
}

getMeanBest <- function(props){
  best <- apply(props[,2:ncol(props)], 1, max)
  return(mean(best))
}


exhaustiveAlignK <- function(t){
  current_best<-getMeanBest(propsFromTable(t))
  current_t<-t
  for (col in seq(2, ncol(t))){
    print(paste0("Column: ", col))
    ord_t<-unique(t[,col])
    if(length(ord_t) == 1){next}
    perm_ords <- combinat::permn(ord_t)
    perm_ords <- perm_ords[2:length(perm_ords)]
    perm_bests<-vector()
    perm_tables<-list()
    for (perm in seq(1, length(perm_ords))){
      new_col <- t[,col]
      for (value in seq(1,length(ord_t))){
        new_col[new_col==ord_t[value]] <- paste0("r",perm_ords[[perm]][value])
      }
      new_col<-str_replace_all(new_col, "r", "")
      temp<-current_t
      temp[,col]<-as.vector(new_col)
      perm_bests[[perm]]<-getMeanBest(propsFromTable(temp))
      perm_tables[[perm]]<-temp
    }
    best<-which(perm_bests==max(perm_bests))
    if(length(best)>1){
      best<-best[[1]]
    }
    if(perm_bests[best] > current_best){
      current_best <- perm_bests[[best]]
      current_t<-perm_tables[[best]]
      print("Found better")
    }
  }
  return(current_t)
}

#returns table of assignment proportions, from a table of cluster labels across runs
propsFromTable <- function(x){
  t<-melt(x, id="Sample")
  tt <- t %>% dplyr::group_by(Sample) %>% 
    select(-variable) %>% 
    add_count(value) %>% 
    group_by(Sample, value) %>% 
    summarise(n=first(n)) %>% 
    arrange(Sample) %>%
    mutate(freq = n/sum(n)) %>%
    dcast(Sample~value, value.var="freq") %>%
    replace(is.na(.), 0)
  return(tt)
}


##############################################################################

align_uml_heuristic <- function(rundir, 
                                 method, 
                                 alg, 
                                 plotDIR){
  #read in table (alternatively, build table from set of runs)
  #t<-read.table("clumpak_adapt/example_data.txt", sep="\t", header=T, stringsAsFactors=F)
  
  #convert run outputs (directory of text files) 
  # to table format
  t <- runs2table(dir=rundir, sampleNamesFromFile=T)
  
  #convert to qlist format
  qlist <- as.qlist(tab2Qlist(t))
  
  #align across K using clumpak-like algorithm
  aligned_Q <- alignK(qlist, type="auto")
  
  #convert back to table format for all runs
  aligned_T <- alignQ2tab(aligned_Q, 
                          t$Sample, zeroIndexing=T)

  # convert to assignment proportions
  assign_props <- propsFromTable(aligned_T)
  
  # Create directory for plots.
  dir.create(plotDIR, showWarnings = FALSE)
  
  mymelt <- melt(aligned_T, id.var = "Sample")
  
  # Plot it
  p <-
    ggplot(mymelt, aes(variable, Sample)) +
    geom_tile(aes(fill = as.factor(value)),
              colour = "white") +
    theme(
      axis.text.y = element_text(size = 8,
                                 colour = "black"),
      axis.text.x = element_blank(),
      axis.title = element_text(size = 14,
                                colour = "black"),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(fill = guide_legend(title = "Species")) +
    labs(x = "Run", y = "Sample ID",
         title = "aligned clusters")
  
  ggsave(
    filename = paste0(method, "_", alg, ".pdf"),
    plot = p,
    device = "pdf",
    path = plotDIR,
    width = 8,
    height = 11,
    units = "in",
    dpi = 300
  )
  
  return(list(aligned_T, assign_props))
}


align_uml_exhaustive <- function(rundir, 
                           method, 
                           alg, 
                           plotDIR){
  #read in table (alternatively, build table from set of runs)
  #t<-read.table("clumpak_adapt/example_data.txt", sep="\t", header=T, stringsAsFactors=F)
  
  #convert run outputs (directory of text files) 
  # to table format
  t <- runs2table(dir=rundir, sampleNamesFromFile=T)

  #convert to qlist format
  qlist <- as.qlist(tab2Qlist(t))
  
  #align across K using clumpak-like algorithm
  aligned_Q <- alignK(qlist, type="auto")
  
  #convert back to table format for all runs
  aligned_T <- alignQ2tab(aligned_Q, 
                          t$Sample, zeroIndexing=T)
  # exhaustive 2nd pass. 
  # This one goes across runs and exhaustively checks
  # if switching any two labels maximizes assignment
  # proportions across samples
  # might take a long time at high K values though!
  aligned2_T <- exhaustiveAlignK(aligned_T)

  # convert to assignment proportions
  assign_props <- propsFromTable(aligned2_T)
  
  #assign_props <- propsFromTable(aligned_T)
  
  
  # Create directory for plots.
  dir.create(plotDIR, showWarnings = FALSE)
  
  mymelt <- melt(aligned2_T, id.var = "Sample")
  
  # Plot it
  p <-
    ggplot(mymelt, aes(variable, Sample)) +
    geom_tile(aes(fill = as.factor(value)),
              colour = "white") +
    theme(
      axis.text.y = element_text(size = 8,
                                 colour = "black"),
      axis.text.x = element_blank(),
      axis.title = element_text(size = 14,
                                colour = "black"),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(fill = guide_legend(title = "Species")) +
    labs(x = "Run", y = "Sample ID",
         title = "aligned clusters")
  
  ggsave(
    filename = paste0(method, "_", alg, ".pdf"),
    plot = p,
    device = "pdf",
    path = plotDIR,
    width = 8,
    height = 11,
    units = "in",
    dpi = 300
  )
  
  return(list(aligned2_T, assign_props))
  
}

dir.create("missDataRuns/plots_maf", showWarnings = FALSE)
dir.create("missDataRuns/Robjects_maf", showWarnings = FALSE)


plotDIR <- "missDataRuns/plots_maf"

for (i in seq(25, 100, 25)) {
  for (j in seq(25, 100, 25)) {
    for (m in seq(0.01, 0.05, 0.02)) {
      
      if (i == 100 & j == 25){
        next
      }
      
      writeLines(paste0("\nDoing missInd ",
                        i,
                        " Pop ",
                        j,
                        " MAF ",
                        m,
                        "...\n"))
      
      writeLines("Doing cmds gapstat...")
      
      cmds_gapstat <-
        align_uml_heuristic(
          rundir = paste0(
            "missDataRuns/final_output/maf/uml/runs_used/missInd",
            i,
            "_Pop",
            j,
            "_maf",
            m,
            "/cmds/gapstat/"
          ),
          method = "cmds",
          alg = "gapstat",
          plotDIR = plotDIR
        )
      
      saveRDS(
        object = cmds_gapstat,
        file = paste0(
          "missDataRuns/Robjects_maf/cmds_gapstat_missInd",
          i,
          "_Pop",
          j,
          "_maf",
          m,
          ".rds"
        )
      )
      
      writeLines("Done!\n")
      
      writeLines("Doing cmds Hierarchical...")
      
      cmds_hier <-
        align_uml_heuristic(
          rundir = paste0(
            "missDataRuns/final_output/maf/uml/runs_used/missInd",
            i,
            "_Pop",
            j,
            "_maf",
            m,
            "/cmds/hierarchical/"
          ),
          method = "cmds",
          alg = "hier",
          plotDIR = plotDIR
        )
      
      saveRDS(
        object = cmds_hier,
        file = paste0(
          "missDataRuns/Robjects_maf/cmds_hier_missInd",
          i,
          "_Pop",
          j,
          "_maf",
          m,
          ".rds"
        )
      )
      
      writeLines("Done with cmds hier!\n")
      
      writeLines("Doing cmds pam (cmds groups)...")
      
      cmds_pam <-
        align_uml_heuristic(
          rundir = paste0(
            "missDataRuns/final_output/maf/uml/runs_used/missInd",
            i,
            "_Pop",
            j,
            "_maf",
            m,
            "/cmds/pam/cmds_groups/"
          ),
          method = "cmds",
          alg = "pam",
          plotDIR = plotDIR
        )
      
      saveRDS(
        object = cmds_pam,
        file = paste0(
          "missDataRuns/Robjects_maf/cmds_pam_missInd",
          i,
          "_Pop",
          j,
          "_maf",
          m,
          ".rds"
        )
      )
      
      writeLines("Done with cmds pam (cmds groups!\n")
      
      writeLines("Doing cmds PAM (proximity scores...")
      
      cmds_pam_prox <-
        align_uml_heuristic(
          rundir = paste0(
            "missDataRuns/final_output/maf/uml/runs_used/missInd",
            i,
            "_Pop",
            j,
            "_maf",
            m,
            "/cmds/pam/prox/"
          ),
          method = "cmds",
          alg = "pam_prox",
          plotDIR = plotDIR
        )
      
      saveRDS(
        object = cmds_pam_prox,
        file = paste0(
          "missDataRuns/Robjects_maf/cmds_pam_prox_missInd",
          i,
          "_Pop",
          j,
          "_maf",
          m,
          ".rds"
        )
      )
      
      writeLines("Done with cmds PAM (Proximity scores)!\n")
      
      writeLines("Doing isomds gapstat...")
      
      isomds_gapstat <-
        align_uml_heuristic(
          rundir = paste0(
            "missDataRuns/final_output/maf/uml/runs_used/missInd",
            i,
            "_Pop",
            j,
            "_maf",
            m,
            "/isomds/gapstat/"
          ),
          method = "isomds",
          alg = "gapstat",
          plotDIR = plotDIR
        )
      
      saveRDS(
        object = isomds_gapstat,
        file = paste0(
          "missDataRuns/Robjects_maf/isomds_gapstat_missInd",
          i,
          "_Pop",
          j,
          "_maf",
          m,
          ".rds"
        )
      )
      
      writeLines("Done with isomds gapstat!\n")
      
      writeLines("Doing isomds hierarchical...")
      
      isomds_hier <-
        align_uml_heuristic(
          rundir = paste0(
            "missDataRuns/final_output/maf/uml/runs_used/missInd",
            i,
            "_Pop",
            j,
            "_maf",
            m,
            "/isomds/hierarchical/"
          ),
          method = "isomds",
          alg = "hier",
          plotDIR = plotDIR
        )
      
      saveRDS(
        object = isomds_hier,
        file = paste0(
          "missDataRuns/Robjects_maf/isomds_hier_missInd",
          i,
          "_Pop",
          j,
          "_maf",
          m,
          ".rds"
        )
      )
      
      writeLines("Done with isomds hierarchical!\n")
      
      writeLines("Doing isomds PAM...")
      
      isomds_pam <-
        align_uml_heuristic(
          rundir = paste0(
            "missDataRuns/final_output/maf/uml/runs_used/missInd",
            i,
            "_Pop",
            j,
            "_maf",
            m,
            "/isomds/pam/"
          ),
          method = "isomds",
          alg = "pam",
          plotDIR = plotDIR
        )
      
      saveRDS(
        object = isomds_pam,
        file = paste0(
          "missDataRuns/Robjects_maf/isomds_pam_missInd",
          i,
          "_Pop",
          j,
          "_maf",
          m,
          ".rds"
        )
      )
      
      writeLines("Done with isomds PAM!\n")
      
      
      ###################################################
      ###################################################
      ## t-SNE
      ###################################################
      ###################################################
      
      for (p in seq(5, 50, 5)) {
        writeLines(paste0(
          "Doing tsne missInd",
          i,
          " Pop",
          j,
          "_maf",
          m,
          " P",
          p,
          "...\n"
        ))
        
        writeLines("Doing tsne gapstat...\n")
        
        tsne_gapstat <-
          align_uml_heuristic(
            rundir = paste0(
              "missDataRuns/final_output/maf/uml/runs_used/missInd",
              i,
              "_Pop",
              j,
              "_maf",
              m,
              "/tsne/gapstat/p",
              p,
              "/"
            ),
            method = paste0("tsneP", p),
            alg = "gapstat",
            plotDIR = plotDIR
          )
        
        saveRDS(
          object = tsne_gapstat,
          file = paste0(
            "missDataRuns/Robjects_maf/tsne_gapstat_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            m,
            "_P",
            p,
            ".rds"
          )
        )
        
        writeLines("Doing tsne hierarchical...\n")
        
        tsne_hier <-
          align_uml_heuristic(
            rundir = paste0(
              "missDataRuns/final_output/maf/uml/runs_used/missInd",
              i,
              "_Pop",
              j,
              "_maf",
              m,
              "/tsne/hierarchical/p",
              p,
              "/"
            ),
            method = paste0("tsneP", p),
            alg = "hier",
            plotDIR = plotDIR
          )
        
        saveRDS(
          object = tsne_hier,
          file = paste0(
            "missDataRuns/Robjects_maf/tsne_hier_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            m,
            "_P",
            p,
            ".rds"
          )
        )
        
        writeLines("Doing tsne PAM...\n")
        
        tsne_pam <-
          align_uml_heuristic(
            rundir = paste0(
              "missDataRuns/final_output/maf/uml/runs_used/missInd",
              i,
              "_Pop",
              j,
              "_maf",
              m,
              "/tsne/pam/p",
              p,
              "/"
            ),
            method = paste0("tsneP", p),
            alg = "pam",
            plotDIR = plotDIR
          )
        
        saveRDS(
          object = tsne_pam,
          file = paste0(
            "missDataRuns/Robjects_maf/tsne_pam_missInd",
            i,
            "_Pop",
            j,
            "_maf",
            m,
            "_P",
            p,
            ".rds"
          )
        )
        
        writeLines(paste0("Done with tsne",
                          i,
                          " Pop",
                          j, "_maf", m,
                          " P",
                          p,
                          "!\n"))
        
        
      }
    }
  }
}
beepr::beep(3)