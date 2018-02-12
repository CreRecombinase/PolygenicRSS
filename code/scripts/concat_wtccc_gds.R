#library(SNPRelate)
library(SeqArray)
library(tidyverse)
## inf_seqgdsf <- dir("~/Desktop/scratch/polyg_scratch/gds",full.names=T,pattern="*seq_wtcc_geno.gds")
## outf <- "/media/nwknoblauch/Data/wtcc_input/combined.gds"
## gn <- gsub(".+snpgds/(.+)_seq_wtcc.+","\\1",inf_seqgdsf)
inf_seqgdsf <- snakemake@input[["inf_seqgds"]]
gn <- snakemake@params[["gn"]]
outf <- snakemake@output[["outf"]]

stopifnot(length(gn)==length(inf_seqgdsf))


gdsl <- map(inf_seqgdsf,seqOpen)

for(i in 1:length(gdsl)){
  x <- gdsl[[i]]
chrom_pos <- seqGetData(x,"$chrom_pos")
# tfiles <- map2(gdsl,chrom_pos,function(x,y){
  seqSetFilter(x,variant.sel=which(y %in% unique_cp))
  tf <- tempfile()
  tf <- seqExport(x,tf)
  refv[i] <- tf
}
## ninf <- map(refv,seqOpen)
## walk(refv,cleanup.gds)
## seqMerge(refv,outf)
## chrom_pos <- map(gdsl,function(tgds){
##     # tgds <- seqOpen(x,readonly=F)
##     sample_n <- index.gdsn(tgds, "sample.id")
##     sample_d <- read.gdsn(sample_n)
##     if(all(!is.na(as.integer(sample_d)))){
##         mgn <- as.character(seqGetData(tgds,"sample.annotation/trait"))
##         n_sample_d <- paste0(sample_d,"_",mgn)
##         add.gdsn(tgds, "sample.id", n_sample_d, replace=TRUE)
##     }
##     #seqClose(tgds)
## })

seqMerge(inf_seqgdsf,outf,storage.option = "LZ4_RA.fast")

