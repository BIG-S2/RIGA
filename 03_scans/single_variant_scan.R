###########################################################
# Individual analysis using STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 02/16/2024
###########################################################
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
## Number of jobs for each chromosome
jobs_num <- get(load("${ADSP_ROOT}/output/jobs_num.Rdata"))
## aGDS directory
agds_dir <- get(load("${ADSP_ROOT}/output/agds_dir.Rdata"))
## input race group of interest
# race <- commandArgs(TRUE)[1]
## Null model
# obj_nullmodel <- get(load(paste0("${ADSP_ROOT}/output/obj.STAAR.ADSP.AD.",race,".20241011.Rdata")))
# obj_nullmodel <- get(load(paste0("${ADSP_ROOT}/output/obj.STAAR.ADSP.AD.20240908.Rdata")))
obj_nullmodel <- get(load("${ADSP_ROOT}/output/matched/obj_null_stratum_all_matched_lpa_pc.pair33_projcell.RData"))
obj_nullmodel$sparse_kins <- TRUE
obj_nullmodel$use_SPA <- FALSE
obj_nullmodel$relatedness <- TRUE
## QC_label
QC_label <- "annotation/info/QC_label7"   ## 2026-10-13 scan on the per-record label directly
## variant_type
variant_type <- "variant"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## output path (refanc variant — separate dir so existing matched/ results stay intact)
output_path <- paste0("${ADSP_ROOT}/output/matched_lpa_pc_pair33_projcell/individual_analysis_qc7/")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
## output file name
output_file_name <- "ADSP_AD_Individual_Analysis"
## input array id from batch file
arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
chr <- which.max(arrayid <= cumsum(jobs_num$individual_analysis_num))
group.num <- jobs_num$individual_analysis_num[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(jobs_num$individual_analysis_num)[chr-1]
}

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)
start_loc <- (groupid-1)*10e6 + jobs_num$start_loc[chr]
end_loc <- start_loc + 10e6 - 1
end_loc <- min(end_loc,jobs_num$end_loc[chr])

a <- Sys.time()
results_individual_analysis <- c()
if(start_loc <= end_loc)
{
  results_individual_analysis <- Individual_Analysis(chr=chr,start_loc=start_loc,end_loc=end_loc,
                                                     genofile=genofile,obj_nullmodel=obj_nullmodel,mac_cutoff=50,
                                                     QC_label=QC_label,variant_type=variant_type,
                                                     geno_missing_imputation=geno_missing_imputation)
}
b <- Sys.time()
b - a

save(results_individual_analysis,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)

