library(stringr)
# for use, assign this to "glib" for the next function to work.
get_glib <- function(){
  # Read in JT's primer csv
  lib <- read_csv("~/Projects/phix_host/PrelimAssay/data/vcpe_oligo_pool.csv", col_names = F) # dir corrected 2/6/23
  tmp <- lib$X1
  lib$gene <- str_split(tmp,"_", simplify = T)[,1]
  lib$codon <- str_split(tmp, "_", simplify = T)[,3]
  lib$mut <- str_split(lib$codon, "-", simplify = T)[,1]
  lib$cod.seq <- str_split(lib$codon, "-", simplify = T)[,2]
  lib$pos <- as.numeric(str_extract(lib$mut, "[0-9]+"))
  glib <- filter(lib, gene == "G1" | gene == "G2") # subset out just G1 and G2 entries
  glib[glib$gene == "G2",]$pos <- glib[glib$gene == "G2",]$pos + 134 #correct G2 pos to the pos in relation to whole gene
  glib$pos <- glib$pos + 1 # make the numbering to include start codon as 1
  aa <- strsplit(glib$mut, split = "[0-9]+") # split out AAs to correct residue numbering
  glib$mut.corr <- NA  # create a new column for AA subs with corrected residue numbering
  for (i in 1:nrow(glib)){
    glib$mut.corr[i] <- paste(aa[[i]][1], glib$pos[i], aa[[i]][2], sep = "")
  }
  return(glib)
}

# # "filepath" to sam file as arg, produces a list of 2 dataframes with count/mut info
analyze_dadasam <- function(sampath){
  sam <- read.table(sampath)
  ct <- str_split(sam[,1], "=", simplify = TRUE)[,2]
  sam$counts <- as.numeric(str_replace(ct, ";", ""))
  sam$mismatches <- as.numeric((str_split(sam[,12], ":", simplify = TRUE)[,3]))
  
  # Save singletons (counts == 1) as separate df
  sam_singletons <- sam[(sam$counts == 1),]
  # delete singleton rows from main df
  sam <- sam[!(sam$counts == 1),]
  # Loop through each read, compare codons
  # first, get ref codons
  wt.seq <- sam$V10[1] # first one, with largest counts, is the wt seq.
  start <- seq(1, nchar(as.character(wt.seq)), 3)
  stop <- start+2
  wt.cod <- str_sub(wt.seq, start, stop)
  
  # create placeholder dataframes (each row will be unique variant seq)
  # and collect them in a list 
  out <- list()
  pos_vector <- c("start", paste("pos", 2:175, sep = ""))
  out$pos_df <- as.data.frame(matrix(NA, nrow = nrow(sam), ncol = 176)) # name mut at each pos
  colnames(out$pos_df) <- pos_vector
  out$sum_df <- as.data.frame(matrix(NA, nrow=nrow(sam), ncol = 6)) # summary df
  colnames(out$sum_df) <- c("var","counts","num_mut","muts","in_lib","not_lib")
  # parse codons and enter muts to pos_df
  # at end of each nth var, counts & mut info into sum_df
  for (n in 1:nrow(sam)){
    out$sum_df$var[n] <- str_split(sam$V1[n], ";", simplify = T)[1] #sq1, sq2, etc
    out$sum_df$counts[n] <- sam$counts[n]
    x <- as.character(sam$V10[n]) # sequence
    x.cod <- str_sub(x, start, stop) # sequence broken into codons
    
    for (i in 1:175){
      if (x.cod[i] == wt.cod[i]){
        out$pos_df[n,i] <- "wt"
      } else {
        if (i %in% glib$pos) {  # start codon and residues 134-135 not in the primer lib
          tmp <- filter(glib, glib$pos == i-1)
          if (x.cod[i] %in% tmp$cod.seq){
            out$pos_df[n,i] <- tmp[tmp$cod.seq == x.cod[i],]$mut.corr
          } else {
            out$pos_df[n,i] <- paste(wt.cod[i], ">", x.cod[i], "@", i, sep = "") #for muts not in lib, indicate codon change
          }
        } else {   #for two residues between G1, G2 fragments
          out$pos_df[n,i] <- paste(wt.cod[i], ">", x.cod[i], "@", i, sep = "")
        }
      }
    }
    if (n==1){
      out$sum_df$num_mut[n] <- 0
      out$sum_df$muts[n] <- 0
      out$sum_df$in_lib[n] <- 0
      out$sum_df$not_lib[n] <- 0
    } else {
      var_n <- out$pos_df[n,] #turn nth row of pos_df into a vector
      muts_all <- var_n[which(var_n!="wt")] #get all non-wt entries from all positions
      muts_notlib <- var_n[which(nchar(var_n)>7)] #muts not from library will be 7 characters
      muts_lib <- muts_all[which(nchar(muts_all)<7)] #muts from lib, indicated by AA sub
      out$sum_df$num_mut[n] <- length(muts_all) #number of muts in each var
      out$sum_df$muts[n] <- paste(muts_all, collapse=";") #all muts
      out$sum_df$in_lib[n] <- length(muts_lib) #number of muts that are in primer lib
      out$sum_df$not_lib[n] <- length(muts_notlib) #number of muts not in primer lib
    }
  }
  return(out)
}


