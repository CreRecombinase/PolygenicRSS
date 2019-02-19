library(SeqArray)
library(EigenH5)
library(SeqSupport)
library(dplyr)
library(ldshrink)
library(rlang)
library(readr)



## input_f_a  <- "/scratch/t.cri.nknoblauch/polyg_scratch/gds/FRAM/grmi/FRAM.chr19.gds"
## input_f_b  <- "/scratch/t.cri.nknoblauch/polyg_scratch/gds/EUR/EUR.chr19.gds"
## threads <- 8
## MAF_filter <- 0.01
## min_snps <- 10000

input_f_a <- snakemake@input[["gds_a"]]
input_f_b <- snakemake@input[["gds_b"]]

output_f_a <- snakemake@output[["gdslist_a"]]
output_f_b <- snakemake@output[["gdslist_b"]]


MAF_filter <- as.numeric(snakemake@params[["MAF"]] %||% 0)
snp_only <- (snakemake@params[["snp_only"]] %||% "T") == "T"
min_snps <- snakemake@params[["min_snps"]]




snp_df_a <- read_df_h5(input_f_a,"SNPinfo")%>% filter(between(MAF,MAF_filter,1-MAF_filter))
snp_df_b <- read_df_h5(input_f_b,"SNPinfo") %>% filter(between(MAF,MAF_filter,1-MAF_filter))

if(snp_only){
    snp_df_a <- filter(snp_df_a,nchar(allele)==3)
    snp_df_b <- filter(snp_df_b,nchar(allele)==3)

}else{
    snp_df_a <- filter(snp_df_a,nchar(allele)>=3)
    snp_df_b <- filter(snp_df_b,nchar(allele)>=3)
}


snp_df_both <- dplyr::inner_join(snp_df_a,snp_df_b,by=c("chr","pos"),suffix=c("_a","_b")) %>%
    filter(abs(flip_alleles(allele_a,allele_b))==1L) %>% distinct(chr,pos,.keep_all=T) %>% arrange(chr,pos)

semi_join(snp_df_a,snp_df_both) %>% distinct(chr,pos,.keep_all=T) %>% arrange(chr,pos) %>%  write_tsv(output_f_a)
semi_join(snp_df_b,snp_df_both) %>% distinct(chr,pos,.keep_all=T) %>% arrange(chr,pos) %>%  write_tsv(output_f_b)
