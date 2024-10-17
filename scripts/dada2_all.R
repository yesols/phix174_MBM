
#######################################################
# G Mutants Growth Assay Seq Analysis
#######################################################

#### command-line:
# 1. fastp to QC
# 2. FLASh to merge forwqrd and reverse reads
# 3. samtools to extract mapped reads with CIGAR 528M (ensure identical length)
#     samtools view -h in.sam | awk '$6 == "528M"' > out.sam
# 4. Reattach header:
#     samtools view -H file.sam > sam_head.txt
#     cat sam_head.txt headerless.sam > out.sam
# 5. Convert to fastq: 
#     picard SamToFastq I=input.sam F=output.fastq
#    Make sure no Ns in fastq (dada will throw an error)
# 
#### R
# 6. dada2 workflow from the fastq
# 7. export to fasta (uniquesToFasta)
# 8. map fasta to sam: bwa mem ../../../phix_anc.fasta G1G2_t35.fasta > G1G2t35_dada.sam
# 9. Remove header from sam files and read it into R as table again, for matching muts and counts



# starting with step 5, *extractedwH.sam files in Projects/phix_host/PrelimAssay/seqdata/filtered/sam/


#setwd("~/Projects/phix_host/PrelimAssay/dada2_all/")
#setwd("~/phix_prelim/")  # Run on server
setwd("~/phix_zoe/")
library(dada2)
library(stringr)

#path <- "~/Projects/phix_host/PrelimAssay/dada2_all/"
#path <- "~/phix_prelim/"
path <- "~/phix_zoe/"
#files <- list.files(paste(path, "input_fastq/", sep = ""))
files <- list.files(paste(path, "converted_fastq/", sep = ""))
#files.in <- paste("input_fastq/", files, sep = "")
files.in <- paste("converted_fastq/", files, sep = "")
#samples <- str_split(files, ".extracted.", simplify = T)[,1]
samples <- str_split(files, "_converted.", simplify = T)[,1]

for (i in 3:length(samples)){
  err.i <- learnErrors(files.in[i], multithread=TRUE)
  dds.i <- dada(files.in[i], err=err.i, multithread=TRUE, OMEGA_A = 1e-10, OMEGA_C = 1e-10, 
                   DETECT_SINGLETONS = TRUE)
  uniquesToFasta(getUniques(dds.i), paste(path, "uniquesToFasta_out/", samples[i], ".fasta", sep = ""))
}
  


# Map fasta from dada2:
#     bwa mem ../../../phix_anc.fasta G1G2_t35.fasta
# then remove first 2 lines (header) so that sam can be read in as a table in R:
#     sed -i '1,2d' file.sam 
#     (-i deletes it from the source file. Otherwise sed sends to stout)
# For multiple files:
# for i in *.sam; do
#      sed -i '1,2d' $i
# done
# -i throwing error. do redirect (wrote a script "sed_sam.sh")


### Edited final sam in sed_sam dir
### R to match mutations and count freq

#### Separate script fn_dadasam.R contains 2 functions:
#### get_glib() to read in JT's primers file to get prot G muts
#### analyze_dadasam() to loop through sam files and get counts, muts info

glib <- get_glib()

dir <- "~/Projects/phix_host/PrelimAssay/dada2_all/sed_sam/"
files <- list.files(dir)
samples <- str_split(files, "\\.", simplify = TRUE)[,1]
outlist <- list()
for (k in 1:length(samples)){
  outlist[[k]] <- analyze_dadasam(paste(dir, files[k], sep = ""))
}
names(outlist) <- samples
yesol_df <- data.frame()
for (i in 1:length(outlist)){
  outlist[[i]]$sum_df$sample <- samples[i]
  yesol_df <- rbind(yesol_df, outlist[[i]]$sum_df[,c("sample","var","counts","num_mut","muts","in_lib","not_lib")])
}
#write.csv(yesol_df, "~/Projects/phix_host/PrelimAssay/dada2_all//R_outputs/dada2result_yesol.csv")
write.csv(yesol_df, "~/Projects/phix_host/PrelimAssay/dada2_all//R_outputs/dada2result2_yesol.csv")



# Zoe's samples
dir <- "~/Projects/phix_host/phix_zoe/sed_sam/"
files <- list.files(dir)
samples <- str_split(files, "\\.", simplify = TRUE)[,1]
zoelist <- list()
for (i in 1:length(samples)){
  zoelist[[i]] <- analyze_dadasam(paste(dir, files[i], sep = ""))
}
names(zoelist) <- samples

zoe_df <- data.frame()
for (i in 1:length(zoelist)){
  zoelist[[i]]$sum_df$sample <- samples[i]
  zoe_df <- rbind(zoe_df, zoelist[[i]]$sum_df[,c("sample","var","counts","num_mut","muts","in_lib","not_lib")])
}
#write.csv(zoe_df, "~/Projects/phix_host/phix_zoe/R_outputs/dada2result_zoe.csv")
write.csv(zoe_df, "~/Projects/phix_host/phix_zoe/R_outputs/dada2result2_zoe.csv")



#################################################
#  New Analysis with revised script/fn   6/3/20 #
#  Create summary df for counts of reads/vars over time
#################################################

#df <- yesol_df
df <- zoe_df

samples <- unique(df$sample)
mut_counts <- data.frame(
  matrix(ncol=11, nrow=length(samples), 
         dimnames=list(samples, c("total_vars","lib_vars","vars_sing","vars_doub","vars_more",
                                  "total_ct","lib_ct","wt_pc","sing_pc","doub_pc","more_pc")))
)
for (i in 1:length(samples)){
  df_i <- df[which(df$sample==samples[i]),]
  df_lib <- df_i[which(df_i$not_lib==0),]
  mut_counts$total_vars[i] <- nrow(df_i)
  mut_counts$total_ct[i] <- sum(df_i$counts)
  mut_counts$lib_vars[i] <- nrow(df_lib)
  mut_counts$lib_ct[i] <- sum(df_lib$counts)
  mut_counts$wt_pc[i] <- 100*df_lib[which(df_lib$var=="sq1"),]$counts/sum(df_lib$counts)
  mut_counts$vars_sing[i] <- nrow(df_lib[which(df_lib$in_lib==1),])
  mut_counts$sing_pc[i] <- 100*sum(df_lib[which(df_lib$in_lib==1),]$counts)/sum(df_lib$counts)
  mut_counts$vars_doub[i] <- nrow(df_lib[which(df_lib$in_lib==2),])
  mut_counts$doub_pc[i] <- 100*sum(df_lib[which(df_lib$in_lib==2),]$counts)/sum(df_lib$counts)
  mut_counts$vars_more[i] <- nrow(df_lib[which(df_lib$in_lib>2),])
  mut_counts$more_pc[i] <- 100*sum(df_lib[which(df_lib$in_lib>2),]$counts)/sum(df_lib$counts)
}

#yesol_summary <- mut_counts
zoe_summary <- mut_counts


### 

