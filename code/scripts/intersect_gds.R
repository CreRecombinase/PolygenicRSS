library(SeqArray)
library(SeqSupport)
library(dplyr)
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
min_snps <- snakemake@params[["min_snps"]]

threads <- as.integer(snakemake@threads)
seqParallelSetup(cluster=threads, verbose=TRUE)
gds_a  <- seqOpen(input_f_a)

gds_b  <- seqOpen(input_f_b)


snp_df_a <- read_SNPinfo_gds(gds_a,MAF=T,alleles=T)%>% filter(between(MAF,MAF_filter,1-MAF_filter))
snp_df_b <- read_SNPinfo_gds(gds_b,MAF=T,alleles=T) %>% filter(between(MAF,MAF_filter,1-MAF_filter))

snp_df_both <- dplyr::inner_join(snp_df_a,snp_df_b,by=c("chr","pos"),suffix=c("_a","_b")) %>% distinct(chr,pos,.keep_all=T) %>% arrange(chr,pos)

semi_join(snp_df_a,snp_df_both) %>% distinct(chr,pos,.keep_all=T) %>% arrange(chr,pos) %>%  write_tsv(output_f_a)
semi_join(snp_df_b,snp_df_both) %>% distinct(chr,pos,.keep_all=T) %>% arrange(chr,pos) %>%  write_tsv(output_f_b)


seqParallelSetup(FALSE)
