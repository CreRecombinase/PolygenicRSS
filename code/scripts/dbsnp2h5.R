library(tidyverse)
dbsnpf <- "/media/nwknoblauch/Data/snp150.txt.gz"
234104110

nvcff <- "/media/nwknoblauch/Data/wtcc_input/bd_19.vcf.gz"
nvc_gdsf <- "/media/nwknoblauch/Data/wtcc_input/bd_19.vcf.gds"
snpgdsVCF2GDS(nvcff,nvc_gdsf,compress.annotation = "LZ4_RA.fast")
gds <- snpgdsOpen(nvc_gdsf)
gdsi <- read_SNPinfo_gds(nvc_gdsf,alleles=T)

nsnp_df <- data_frame(SNP=read.gdsn(index.gdsn(gds,"snp.rs.id")),
           pos=read.gdsn(index.gdsn(gds,"snp.position")),
           chr=read.gdsn(index.gdsn(gds,"snp.chromosome")))



map_files <- dir("/media/nwknoblauch/Data/1kg/CEU",full.names = T)
tm_df <- read_delim(map_files[1],delim="\t",trim_ws = T)
