library(tidyverse)
library(SNPRelate)
library(gdsfmt)
library(SeqSupport)
library(SeqArray)
library(EigenH5)



#dgn_covarbf <- "

DGN_geno_gds <- "/run/media/nwknoblauch/Data/polyg_scratch/gds/dgn_geno.gds"
dgn_expf <-  "/run/media/nwknoblauch/Data/DGN/Lev/data_used_for_eqtl_study/cis_data.txt"
dgn_h5 <- "/run/media/nwknoblauch/Data/polyg_scratch/h5/ALL_dgn_geno.h5"
DGN_bedf <- snakemake@input[["dgn_bed"]]
DGN_famf  <- snakemake@input[["dgn_fam"]]
DGN_bimf <-  snakemake@input[["dgn_bim"]]
dgn_expf <- snakemake@input[["dgn_expf"]]
DGN_geno_gds <- snakemake@output[["gdsf"]]
dgn_h5 <- snakemake@output[["h5f"]]
SNPRelate::snpgdsBED2GDS(bed.fn = DGN_bedf,fam.fn = DGN_famf,bim.fn = DGN_bimf,out.gdsfn = DGN_geno_gds,family = F)


makeSeq <- function(inf,readonly=TRUE){
  mgds <- safely(seqOpen,otherwise=NA)(inf)
  gdsfmt::showfile.gds(T)
  if(is.na(mgds$result)){
    tempf <- tempfile()
    inf <- seqSNP2GDS(inf,tempf,storage.option ="LZ4_RA.fast")
  }
  return(seqOpen(inf,readonly = readonly))
}

tgds <- makeSeq(DGN_geno_gds)
si_df <- read_SNPinfo_gds(tgds) %>% mutate(chr=as.integer(chr))
si_df <- mutate(si_df,t_snp_id=1:n())
sample_info <- tibble::data_frame(sample_id=seqGetData(tgds,"sample.id"))

#
#
# dgn_gds <- gdsfmt::openfn.gds(DGN_geno_gds,readonly = F)
# dgn_leg <- data_frame(rsid=read.gdsn(index.gdsn(dgn_gds,"snp.id")),
#                       chrom=as.integer(read.gdsn(index.gdsn(dgn_gds,"snp.chromosome"))),
#                       pos=as.integer(read.gdsn(index.gdsn(dgn_gds,"snp.position"))),
#                       allele=as.character(read.gdsn(index.gdsn(dgn_gds,"snp.allele"))))%>% mutate(ind=1:n())
# dgn_leg <- separate(dgn_leg,allele,into = c("ref","alt"),convert = T)
# ndgn_leg <- distinct(dgn_leg,chrom,pos,.keep_all=T) %>% filter(pos>0)
# dgn_ids <- read.gdsn(index.gdsn(dgn_gds,"sample.id"))
# #dgn_snp <- read.gdsn(index.gdsn(dgn_gds,"genotype"))
# dgn_snp <- dgn_snp[,ndgn_leg$ind]
# if(file.exists(DGN_geno_h5)){
#   file.remove(DGN_geno_h5)
# }
# ndgn_leg <- select(ndgn_leg,-ind)
# write_df_h5(df=ndgn_leg,groupname = "SNPinfo",outfile = DGN_geno_h5)
#



sample_info <- mutate(sample_info,new_sample_id=gsub("WG[0-9]+-DNA.+_[0-9]+_(LD[0-9]+)-.+","\\1",x = sample_id),index=1:n())

dgn_exp <- read_delim(dgn_expf,col_names = T,delim = "\t")
dgn_exp <- rename(dgn_exp,Id=X1) %>% select(Id)
n_sample_info <-inner_join(sample_info,dgn_exp,by=c("new_sample_id"="Id"))
snp_df <- distinct(si_df,chr,pos,.keep_all=T) %>% filter(pos>0)

seqSetFilter(tgds,variant.sel = snp_df$t_snp_id,sample.sel=n_sample_info$index)
gds2hdf5(tgds,dgn_h5)
