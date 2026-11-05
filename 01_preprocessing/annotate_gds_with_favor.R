###########################################################
# The code was modified by the one composed
#       by Xihao Li, Zilin Li.
# https://github.com/xihaoli/STAARpipeline-Tutorial
###########################################################

rm(list=ls())

##### Parse parameters #####
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
        cat("\n")
        cat("Usage: Rscript --vanilla annotate_gds_with_favor.R <CHR> \n")
	cat("\tNote:chr <CHR> gds file will be modified!!!\n")
	cat("\tInput:\n")
	cat("\t\t<CHR_GDS> is a GDS file for a chromosome..\n")
	cat("\t\t<CHR> is a number.\n")
	cat("\t\tThe folder contains Anno_chr_chr<CHR>_*.csv.\n")
        q(save="no")
}

gc()

##########################################################################
#           Input
##########################################################################

### gds file
chr <- as.numeric(args[1])
# gds <- paste0("${ADSP_ROOT}/processed/r5_chr",chr,".split.gds")
gds <- paste0("~/overflow/adsp/r5_chr",chr,".gds")
### annotation file (output of Annotate.R)
FAVOR_Anno_file <- paste0("${ADSP_ROOT}/output/chr",chr,"/Anno_chr",chr,"_STAARpipeline.csv")

###########################################################################
#           Main Function 
###########################################################################

### load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(readr)


## open GDS
genofile <- seqOpen(gds, readonly = FALSE)

### read annotation data
FunctionalAnnotation <- read_csv(FAVOR_Anno_file,
col_types=list(col_character(),col_double(),col_double(),col_double(),col_double(),
col_double(),col_double(),col_double(),col_double(),col_double(),
col_character(),col_character(),col_character(),col_double(),col_character(),
col_character(),col_character(),col_character(),col_character(),col_double(),
col_double(),col_character()))

dim(FunctionalAnnotation)

## rename colnames
colnames(FunctionalAnnotation)[2] <- "apc_conservation"
colnames(FunctionalAnnotation)[7] <- "apc_local_nucleotide_diversity"
colnames(FunctionalAnnotation)[9] <- "apc_protein_function"

recursive_delete <- function(node) {
  # Get the list of sub-nodes
  subnodes <- ls.gdsn(node)

  # Recursively delete each sub-node
  for(subnode_name in subnodes) {
    subnode <- index.gdsn(node, subnode_name)
    recursive_delete(subnode)  # Recursive call
  }

  # After all sub-nodes are deleted, delete the node itself
  delete.gdsn(node)
}

Anno.folder <- index.gdsn(genofile, "annotation/info")

# if ("FunctionalAnnotation" %in% ls.gdsn(Anno.folder)) {
  # functional_annotation_node <- index.gdsn(genofile, "annotation/info/FunctionalAnnotation")
  # recursive_delete(functional_annotation_node)
# }
add.gdsn(Anno.folder, "FunctionalAnnotation", val=FunctionalAnnotation, compress="LZMA_RA", replace=TRUE)



seqClose(genofile)
