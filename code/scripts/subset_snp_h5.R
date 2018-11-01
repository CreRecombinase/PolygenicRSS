
## save.image("subset.RData")
## stop()

library(tidyverse)
library(EigenH5)
library(SeqSupport)
library(LDshrink)


inf_a <- snakemake@input[["gdsf_a"]]
inf_b <- snakemake@input[["gdsf_b"]]
dosagef <- snakemake@input[["dosagef"]]


if(is.null(inf_b)){
    inf_b <- inf_a
}

ldetectf <- snakemake@input[["ldetectf"]]

trait_indf <- snakemake@output[["trait_indf"]]
trait_snpf  <- snakemake@output[["trait_snpf"]]

ld_indf <- snakemake@output[["ld_indf"]]
ld_snpf <- snakemake@output[["ld_snpf"]]


chromosome <- snakemake@params[["chrom"]]
AF_cutoff <- as.numeric(snakemake@params[["AF"]])
SNPCT <- as.integer(as.numeric(snakemake@params[["SNPCT"]]))
N  <- as.integer(as.numeric(snakemake@params[["N"]]))



stopifnot(file.exists(ldetectf),!is.null(ldetectf))
break_df <- read_delim(ldetectf,delim="\t",trim_ws = T) %>%
    mutate(chr=as.integer(gsub("chr","",chr)),start=as.integer(start),stop=as.integer(stop)) %>%
    filter(!is.na(chr)) %>%
    filter(chr %in% 1:22) %>%
    mutate(region_id=1:n())

stopifnot(!is.na(AF_cutoff),
          !is.na(SNPCT),
          !is.na(N),
          length(AF_cutoff)>0,
          length(SNPCT)>0,
          length(N)>0)
#gwas_inds <- sort(sample(1:N,ceiling(N/2),replace = F))
inds_a  <- read_vector_h5(inf_a,"SampleInfo/sample_id")
dosage_dims <- dim_h5(inf_a,"dosage")
SNPfirst  <- length(inds_a)==dosage_dims[2]


inds_b  <- read_vector_h5(inf_b,"SampleInfo/sample_id")

if(N>0){
    sN <- min(c(length(inds_a),length(inds_b)))
    subset_ind_a <- data_frame(sample_id=1:N)
    write_df_h5(subset_ind_a,trait_indf,"SampleInfo")
    write_df_h5(subset_ind_a,ld_indf,"SampleInfo")
}else{
    subset_ind_a <- integer()
    subset_df_a <- data_frame(sample_id=1:length(inds_a))
    subset_df_b <- data_frame(sample_id=1:length(inds_b))
    write_df_h5(subset_df_a,trait_indf,"SampleInfo")
    write_df_h5(subset_df_b,ld_indf,"SampleInfo")
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
    all_alleles  <- outer(c("A","C","T","G"),c("A","C","T","G"),function(x,y)paste(x,y,sep=","))
    all_alleles <- data_frame(allele=c(all_alleles[upper.tri(all_alleles)], all_alleles[lower.tri(all_alleles)]))

    subsnp_df <- read_df_h5(inf, "SNPinfo") %>%
        filter(chr%in%tchromosome)  %>%
        filter(!((chr==6)&between(pos, 27866528, 34775446))) %>%
        distinct(chr, pos, .keep_all=T)  %>%
        semi_join(all_alleles)



    if(is.null(subsnp_df$MAF)||(length(subset_ind)>0)){
        if(length(subset_ind)>0){
            if(SNPfirst){
                snp_idl <- BBmisc::chunk(subsnp_df$snp_id, chunk.size=10000) %>% map(~list(subset_rows=.x,
                                                                                           filename=inf,subset_cols=subset_ind,
                                                                                           datapath="dosage"))
            }else{
                snp_idl <- BBmisc::chunk(subsnp_df$snp_id, chunk.size=10000) %>% map(~list(subset_cols=.x,
                                                                                           filename=inf,subset_rows=subset_ind,
                                                                                           datapath="dosage"))
            }
        }else{
            if(SNPfirst){
                snp_idl <- BBmisc::chunk(subsnp_df$snp_id, chunk.size = 10000) %>% map(~list(subset_rows=.x,
                                                                                             filename=inf,
                                                                                             datapath="dosage"))
            }else{
                snp_idl <- BBmisc::chunk(subsnp_df$snp_id, chunk.size = 10000) %>% map(~list(subset_cols=.x,
                                                                                             filename=inf,
                                                                                             datapath="dosage"))
            }
            subsnp_df <- mutate(subsnp_df, AF=calc_af_h5(snp_idl,list(SNPfirst=SNPfirst)))
        }
    }else{
        subsnp_df <-rename(subsnp_df, AF=MAF)
    }
    subsnp_df <- subsnp_df %>% filter(pmin(AF, 1 - AF)>AF_cutoff)
    subsnp_df <-assign_snp_block(subsnp_df, break_df, assign_all = F) %>% filter(!is.na(region_id))
    return(subsnp_df)
}


snp_a <- subset_snp_df(inf_a, chromosome, AF_cutoff=AF_cutoff, subset_ind_a)
if(all.equal(subset_ind_a, subset_ind_b) && inf_a == inf_b){
    if(SNPCT>0){
        snp_a <- sample_n(snp_a, SNPCT, replace=F) %>% arrange(chr, pos)
    }
    write_df_h5(snp_a, trait_snpf, "SNPinfo")
    write_df_h5(snp_a, ld_snpf, "SNPinfo")
#    snp_b <- snp_a
}else{
    snp_b <- subset_snp_df(inf_b, chromosome, AF_cutoff = AF_cutoff, subset_ind_b)


    both_snp <- inner_join(snp_a, snp_b, by=c("chr", "pos"), suffix=c("_a", "_b"))

    subsnp_a <-semi_join(snp_a, both_snp)
    subsnp_b <- semi_join(snp_b, both_snp)

    bdf <-flip_allele_exp(allele_a = subsnp_a$allele, allele_b=subsnp_b$allele)
    both_snp <- filter(both_snp, !is.na(bdf))

    p <- nrow(both_snp)

    stopifnot(p >= SNPCT)
    if(SNPCT>0){
        both_snp <- sample_n(both_snp, SNPCT, replace=F) %>% arrange(chr, pos)
    }
    stopifnot(!is.unsorted(both_snp$snp_id_a),
              !is.unsorted(both_snp$snp_id_b))


    trait_snp <- select(both_snp,AF=AF_a,allele=allele_a,chr,pos,snp_id=snp_id_a)
    LD_snp <- select(both_snp,AF=AF_b,allele=allele_b,chr,pos,snp_id=snp_id_b)


    write_df_h5(trait_snp,trait_snpf,"SNPinfo")
    write_df_h5(LD_snp,ld_snpf,"SNPinfo")
}
