library(pophelper)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(combinat)
library(stringr)

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
    #print(qdf)
    ret[[count]] <- qdf[,-c(1)]
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

aligned2_T <- exhaustiveAlignK(aligned_T)

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


setwd ("~/Dropbox/Academic/Manuscripts/Ongoing_Collabs/Martin_etal_BOX2_phy/")

#read in table (alternatively, build table from set of runs)
#t<-read.table("clumpak_adapt/example_data.txt", sep="\t", header=T, stringsAsFactors=F)

#convert run outputs (directory of text files) to table format
t <- runs2table(dir="clumpak_adapt/isomds_hier/", sampleNamesFromFile=T)

#convert to qlist format
qlist <- as.qlist(tab2Qlist(t))

#align across K using clumpak-like algorithm
aligned_Q <- alignK(qlist, type="auto", maxiter=10000, threshold=1e-20)

#convert back to table format for all runs
aligned_T <-alignQ2tab(aligned_Q, t$Sample, zeroIndexing=T)

#exhaustive 2nd pass. This one goes across runs and exhaustively checks 
#if switching any two labels maximizes assignment proportions across samples
#might take a long time at high K values though!
aligned2_T <- exhaustiveAlignK(aligned_T)

#convert to assignment proportions
assign_props <- propsFromTable(aligned2_T)




mymelt <- melt(aligned_T, id.var = "Sample")
# Plot it
p <-
  ggplot(mymelt, aes(variable, Sample)) + geom_tile(aes(fill = as.factor(value)),
                                                      colour = "white") +
  theme(
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.text.x = element_blank(),
    axis.title = element_text(size = 14, colour = "black"),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + guides(fill = guide_legend(title = "Species")) + 
  labs(x = "Run", y = "Sample ID", title = "aligned clusters")



