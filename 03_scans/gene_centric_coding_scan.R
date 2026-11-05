#####################################################################
# Gene-centric analysis for coding rare variants using STAARpipeline
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
library(MultiSTAAR)
###########################################################
#           User Input
###########################################################
## aGDS directory
jobs_num <- get(load("${ADSP_ROOT}/output/jobs_num.Rdata"))
## aGDS directory
agds_dir <- get(load("${ADSP_ROOT}/output/agds_dir.Rdata"))
## Null model
# obj_nullmodel <- get(load("${ADSP_ROOT}/output/obj.STAAR.ADSP.AD.20240908.Rdata"))
# obj_nullmodel <- get(load(paste0("${ADSP_ROOT}/output/obj.STAAR.ADSP.AD.",race,".20241011.Rdata")))
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
output_path <- "${ADSP_ROOT}/output/matched_lpa_pc_pair33_projcell/coding/"
# output_path <- "${ADSP_ROOT}/output/matched_lpa_pc/SPA/coding_qc5/"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

## output file name
output_file_name <- "ADSP_AD_Coding"
## input array id from batch file
arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
## gene number in job
gene_num_in_array <- 50
group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)
sum(group.num.allchr)

chr <- which.max(arrayid <= cumsum(group.num.allchr))
print(chr)
group.num <- group.num.allchr[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(group.num.allchr)[chr-1]
}

genes_info_chr <- genes_info[genes_info[,2]==chr,]
sub_seq_num <- dim(genes_info_chr)[1]

if(groupid < group.num)
{
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):(groupid*gene_num_in_array)
}else
{
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):sub_seq_num
}

## exclude large coding masks
if(arrayid==57)
{
  sub_seq_id <- setdiff(sub_seq_id,840)
}

if(arrayid==112)
{
  sub_seq_id <- setdiff(sub_seq_id,c(543,544))
}

if(arrayid==113)
{
  sub_seq_id <- setdiff(sub_seq_id,c(575,576,577,578,579,580,582))
}

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

genes <- genes_info

results_coding <- c()
for(kk in sub_seq_id)
{
  print(kk)
  gene_name <- genes_info_chr[kk,1]
  results <- Gene_Centric_Coding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                 rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                 QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  results_coding <- append(results_coding,results)
}

save(results_coding,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)

