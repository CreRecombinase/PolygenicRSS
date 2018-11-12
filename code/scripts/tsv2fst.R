library(readr)
library(fst)

gwasif <- snakemake@input[["gwasf"]]
gwasof <- snakemake@output[["gwasf"]]


gwas_cols <-readr::cols_only("chrom"="i",
                        "pos"="i",
                        "snp"="c",
                        "ref_allele"="c",
                        "alt_allele"="c",
                        "beta_hat"="d",
                        "se"="d",
                        "sample_size"="d")


gwas_df <-read_delim(gwasif,delim="\t",col_types=gwas_cols)

write_fst(gwas_df,path=gwasof,compress=100)
