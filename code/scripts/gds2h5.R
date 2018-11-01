library(SeqArray)
library(EigenH5)
library(SeqSupport)
library(tidyverse)
#data("break_df")
#inf <- "~/Desktop/scratch/polyg_scratch/gds/scombined_19.gds"
inf <- snakemake@input[["input_gds"]]
outf <- snakemake@output[["outf"]]
#mapf <- snakemake@input[["mapf"]]
## save.image()
## stop()

makeSeq <- function(inf,readonly=TRUE){
  mgds <- safely(seqOpen,otherwise=NA)(inf)
  gdsfmt::showfile.gds(T)
  if(is.na(mgds$result)){
    tempf <- tempfile()
    inf <- seqSNP2GDS(inf,tempf,storage.option ="LZ4_RA.fast")
  }
  return(seqOpen(inf,readonly = readonly))
}

tgds <- makeSeq(inf)
cat("Reading SNP info\n")
si_df <- read_SNPinfo_gds(tgds,alleles=T,MAF=T,region_id=F,map = F,info=T,more=list(rs="annotation/id")) %>% mutate(chr=as.integer(chr)) %>% mutate(snp_id=1:n(),nsnp_id=snp_id)
cat("Beginning HDF5 Conversion\n")
SeqSupport::gds2hdf5(tgds,outf,si_df)
pl <- snakemake@wildcards
pl <- as_data_frame(pl[names(pl)!=""])
write_df_h5(pl,outf, "Wildcards")
