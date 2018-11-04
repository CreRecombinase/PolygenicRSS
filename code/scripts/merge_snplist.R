library(readr)
library(dplyr)
## gwasf <- "/scratch/t.cri.nknoblauch/polyg_scratch/summary_statistics/bc_mchc_summary_statistics.tsv.gz"
## vcff <- "/scratch/t.cri.nknoblauch/polyg_scratch/vcf/EUR/snplist/EUR.chr2.txt"
gwasf <- snakemake@input[["gwasf"]]
vcff <- snakemake@input[["vcff"]]
output_f <- snakemake@output[["outputf"]]


gwas_cols <-readr::cols_only("chrom"="c",
                        "pos"="i",
                        "snp"="c",
                        "ref_allele"="c",
                        "alt_allele"="c",
                        "beta_hat"="d",
                        "se"="d",
                        "sample_size"="d")

vcf_cn  <- c("chrom",
             "pos",
             "snp",
             "ref_allele",
             "alt_allele",
             "score",
             "filter",
             "INFO")

vcf_cols <- readr::cols("chrom"="c",
                        "pos"="i",
                        "snp"="c",
                        "ref_allele"="c",
                        "alt_allele"="c")


gwas_df <-read_delim(gwasf,delim="\t",col_types=gwas_cols)

snp_df <- read_delim(vcff,delim="\t",col_names=vcf_cn,col_types=vcf_cols)

b_df <- inner_join(gwas_df,snp_df)

cat("Fraction of SNPs found :",(nrow(b_df)/nrow(filter(gwas_df,chrom==snp_df$chrom[1]))),"\n")

write_delim(b_df,output_f,delim="\t")
