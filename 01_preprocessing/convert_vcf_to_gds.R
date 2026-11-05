#!/usr/bin/env Rscript

# LOAD REQUIRED LIBRARIES
library(SeqVarTools)
library(gdsfmt)
library(SeqArray)
library(GenomicRanges)

# Get user arguments
args = commandArgs(trailingOnly=TRUE)
vcf=args[1]
gds=args[2]
pN=as.integer(args[3])

## Perform the conversion
seqVCF2GDS(vcf, gds, verbose=FALSE, parallel=pN, storage.option="LZMA_RA")

## Now we need to inject "PASS" into every entry for variants to be analysed from here:
bedFile<- args[4]               # e.g. "processed/r5_chr1.bed"

# Open the GDS for writing
genofile<-seqOpen(gds, readonly = FALSE)
genofile
# extract and assign annotation data to Anno.folder (gdsfmt)
Anno.folder <- index.gdsn(genofile, "annotation/info")

# Add a new GDS node to the GDS file. (gdsfmt)
# generate vector we can use to count variants to apply function through rep
position <- as.integer(seqGetData(genofile, "position"))
add.gdsn(Anno.folder, "QC_label", val=factor(rep("PASS", length(position))), compress="LZMA_ra", closezip=TRUE)

# Read BED file
#    (Assumes BED has no header, columns: chrom, start, end)
bed <- read.table(bedFile, header=FALSE)
colnames(bed) <- c("chrom", "start", "end")

# Get existing variant annotation from GDS
chrom <- seqGetData(genofile, "chromosome")
pos   <- seqGetData(genofile, "position")
varid <- seqGetData(genofile, "variant.id")

# Convert to GenomicRanges for easier subsetting
varGR <- GRanges(seqnames = chrom,
                 ranges   = IRanges(start=pos, width=1))

# NOTE: If your BED is zero-based, you might want bed$start+1. 
# E.g.:
bedGR <- GRanges(seqnames = sub("^chr", "", bed$chrom),
                 ranges   = IRanges(start=bed$start+1,
                                    end  = bed$end))

# Identify overlap
ov <- findOverlaps(varGR, bedGR)
keepIndices <- unique(queryHits(ov))  # row indices for variants in the intervals
keepVarIDs  <- varid[keepIndices]     # actual variant IDs

# Extract the *existing* QC_label
oldQC  <- seqGetData(genofile, "annotation/info/QC_label")  
# This should be a factor with levels "PASS" (and possibly "FAIL"?

# Combine conditions
# Condition A: oldQC == "PASS"
# Condition B: variant is in bedOverlapIndices
# Both must be TRUE => final "PASS". Otherwise => "FAIL".
isPassOldQC  <- (oldQC == "PASS")
isInBedRange <- logical(length(varid))
isInBedRange[keepIndices] <- TRUE

qualifiedIndices <- which(isPassOldQC & isInBedRange)

# Build QC_label2
# Initialize everything to "FAIL"
qcVec2 <- rep("FAIL", length(varid))

# Mark the qualified variants as "PASS"
qcVec2[qualifiedIndices] <- "PASS"

# Make it a factor with known levels
qcVec2 <- factor(qcVec2, levels = c("FAIL","PASS"))

# Add as new GDS node in annotation/info
add.gdsn(
  node = Anno.folder,
  name = "QC_label2",
  val  = qcVec2,
  compress = "LZMA_ra",
  closezip = TRUE
)

delete.gdsn(index.gdsn(Anno.folder, "QC_label"))

#10) Close the GDS
seqClose(genofile)
  

