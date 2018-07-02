library(tidyverse)
#library(EigenH5)
library(SeqArray)
library(SeqSupport)
library(LDshrink)
# setwd("~/Dropbox/PolygenicRSS/code/snakemake_files/")
# load("ssSNP.RData")


## save.image()
## stop()
inf_a <- snakemake@input[["gdsf_a"]]
inf_b <- snakemake@input[["gdsf_b"]]
ldetectf <- snakemake@input[["ldetectf"]]
outf_a <- snakemake@output[["outf_a"]]
outf_b <- snakemake@output[["outf_b"]]

outf_gwas <- snakemake@output[["gwaso"]]
outf_ld <- snakemake@output[["ldo"]]

chromosome <- snakemake@params[["chrom"]]
AF_cutoff <- as.numeric(snakemake@params[["AF"]])
SNPCT <- as.integer(as.numeric(snakemake@params[["SNPCT"]]))
N  <- as.integer(as.numeric(snakemake@params[["N"]]))




stopifnot(file.exists(ldetectf),!is.null(ldetectf))
break_df <- read_delim(ldetectf,delim="\t",trim_ws = T) %>% mutate(chr=as.integer(gsub("chr","",chr))) %>% filter(!is.na(chr)) %>% filter(chr %in% 1:22) %>% mutate(region_id=1:n())

stopifnot(!is.na(AF_cutoff),!is.na(SNPCT),!is.na(N),length(AF_cutoff)>0,
          length(SNPCT)>0,length(N)>0)
                                        #gwas_inds <- sort(sample(1:N,ceiling(N/2),replace = F))
tgds_a <- seqOpen(inf_a)
if(inf_a==inf_b){
    tgds_b <- tgds_a
}else{
    tgds_b <- seqOpen(inf_b)
}
inds_a <- seqGetData(tgds_a,"sample.id")
inds_b <- seqGetData(tgds_b,"sample.id")
if(N>0){
    sN <- min(c(length(inds_a),length(inds_b)))
    subset_ind_a <- 1:N
    saveRDS(subset_ind_a,outf_gwas)
    saveRDS(subset_ind_a,outf_ld)
}else{
    subset_ind_a <- integer()
    saveRDS(1:length(inds_a),outf_gwas)
    saveRDS(1:length(inds_b),outf_ld)
}
subset_ind_b <- subset_ind_a


subset_snp_df <- function(inf,chromosome,AF_cutoff=0,subset_ind){

    cat("chromosome:",chromosome)
    if(grepl("-",chromosome)){
        chr_bounds <- as.integer(strsplit(chromosome,"-")[[1]])
        tchromosome <- chr_bounds[1]:chr_bounds[2]
    }else{
        if(grepl(",",chromosome)){
            tchromosome <- as.integer(strsplit(chromosome,",")[[1]])
        }else{
            tchromosome <- as.integer(chromosome)
        }
    }
    stopifnot(all(!is.na(tchromosome)))

    seqSetFilterChrom(inf, tchromosome)
    if(length(subset_ind)>0){
        seqSetFilter(inf,sample.id=subset_ind,action="intersect")
    }
    AF <- seqAlleleFreq(inf)
    mAF <- pmin(AF,1-AF)
    seqSetFilter(inf, variant.sel=(mAF>=AF_cutoff), action="intersect")
    subsnp_df <- read_SNPinfo_gds(inf,MAF=T,allele=T) %>% mutate(chr=as.integer(chr))  %>%
        filter(!((chr==6)&between(pos, 27866528, 34775446))) %>%
        distinct(chr, pos, .keep_all=T)
    seqSetFilter(inf,variant.id=subsnp_df$snp_id,action="intersect")

    subsnp_df <-assign_snp_block(subsnp_df,break_df,assign_all = F) %>% filter(!is.na(region_id))
    return(subsnp_df)
}


snp_a <- subset_snp_df(tgds_a,chromosome,AF_cutoff=AF_cutoff,subset_ind_a)
if(inf_a==inf_b){
    snp_b <- snp_a
}else{
    snp_b <- subset_snp_df(tgds_b,chromosome,AF_cutoff=AF_cutoff,subset_ind_b)
}

both_snp <- inner_join(snp_a,snp_b,by=c("chr","pos"),suffix=c("_a","_b"))

subsnp_a <-semi_join(snp_a,both_snp)
subsnp_b <- semi_join(snp_b,both_snp)

## bdf <-flip_allele_exp(allele_a = subsnp_a$allele,allele_b=subsnp_b$allele)
## both_snp <- filter(both_snp,!is.na(bdf))

p <- nrow(both_snp)

stopifnot(p >= SNPCT)
if(SNPCT>0){
    both_snp <- sample_n(both_snp,SNPCT,replace=F) %>% arrange(chr,pos)
}
stopifnot(!is.unsorted(both_snp$snp_id_a),
          !is.unsorted(both_snp$snp_id_b))


trait_snp <- select(both_snp,AF=MAF_a,SNP=SNP_a,allele=allele_a,chr,pos,snp_id=snp_id_a)
LD_snp <-    select(both_snp,AF=MAF_b,SNP=SNP_b,allele=allele_b,chr,pos,snp_id=snp_id_b)


write_delim(trait_snp,outf_a,delim="\t")
write_delim(LD_snp,outf_b,delim="\t")
