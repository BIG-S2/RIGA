## D_S (project x centre cells as strata, all subjects; 2026-09-22): copy of the pair33 driver with the D_S null model,
## QC_label7 (the corrected per-variant label, as the coding chain) and its own output directory.
#####################################################################
# Gene-centric analysis for noncoding rare variants of ncRNA
# genes using STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 02/17/2024
#####################################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
## aGDS directory
agds_dir <- get(load("${ADSP_ROOT}/output/agds_dir.Rdata"))
## Null model
# obj_nullmodel <- get(load("${ADSP_ROOT}/output/obj.STAAR.ADSP.AD.20240908.Rdata"))
obj_nullmodel <- get(load("${ADSP_ROOT}/output/matched/obj_null_stratum_all_matched_lpa_pc.pair33_projcell.RData"))
obj_nullmodel$sparse_kins <- TRUE
obj_nullmodel$use_SPA <- FALSE
obj_nullmodel$relatedness <- TRUE
## QC_label
QC_label <- "annotation/info/QC_label7"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("${ADSP_ROOT}/output/Annotation_name_catalog.Rdata"))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## output path
output_path <- "${ADSP_ROOT}/output/matched_lpa_pc_pair33_projcell/ncRNA/"
# output_path <- "${ADSP_ROOT}/output/matched_lpa_pc_nofilter_pair33/SPA/ncRNA/"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

## output file name
output_file_name <- "ADSP_AD_ncRNA"
## input array id from batch file
arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
## gene number in job
gene_num_in_array <- 100
group.num.allchr <- ceiling(table(ncRNA_gene[,1])/gene_num_in_array)
sum(group.num.allchr)

chr <- which.max(arrayid <= cumsum(group.num.allchr))
group.num <- group.num.allchr[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(group.num.allchr)[chr-1]
}

ncRNA_gene_chr <- ncRNA_gene[ncRNA_gene[,1]==chr,]
sub_seq_num <- dim(ncRNA_gene_chr)[1]

if(groupid < group.num)
{
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):(groupid*gene_num_in_array)
}else
{
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):sub_seq_num
}

## exclude large ncRNA masks
if(arrayid==117)
{
  sub_seq_id <- setdiff(sub_seq_id,53)
}

if(arrayid==218)
{
  sub_seq_id <- setdiff(sub_seq_id,19)
}

if(arrayid==220)
{
  sub_seq_id <- setdiff(sub_seq_id,c(208,274))
}

if(arrayid==221)
{
  sub_seq_id <- setdiff(sub_seq_id,311)
}

if(arrayid==156)
{
  sub_seq_id <- setdiff(sub_seq_id,41)
}

if(arrayid==219)
{
  sub_seq_id <- setdiff(sub_seq_id,103)
}

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

results_ncRNA <- c()
for(kk in sub_seq_id)
{
  print(kk)
  gene_name <- ncRNA_gene_chr[kk,2]
  results <- c()
  results <- try(ncRNA(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                       rare_maf_cutoff=0.01,rv_num_cutoff=2,
                       QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                       Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                       Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
  results_ncRNA <- rbind(results_ncRNA,results)
}

save(results_ncRNA,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)

