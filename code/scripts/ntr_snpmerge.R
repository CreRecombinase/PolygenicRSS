library(tidyverse)
vcf_inf <- paste0("/scratch/t.cri.nknoblauch/polyg_scratch/vcf_snplist/EUR/",1:22,".txt")
ntr_inf <- "~/Downloads/PolygenicRSS/data/Snakemake_inputs/ntr_snps.txt"
outf <- paste0("~/Downloads/PolygenicRSS/code/snakemake_files/",1:22,"_matched_kg.txt.gz")

outf <- snakemake@output[["match_df"]]
vcf_inf <- snakemake@input[["vcff"]]
ntr_inf <- snakemake@input[["ntrf"]]






allmap <- map_df(vcf_inf,~semi_join(ntr_df,read_delim(.x,delim="\t",col_names=c("chrom","pos","snp","ref","alt"),col_types=cols(chrom="i",pos="i",snp="c",ref="c",alt="c"))))  %>% mutate(ld_snp_id=1:n())

split(allmap,allmap$chrom) %>% walk2(outf,~write_delim(.x,.y,delim="\t"))
##  inner_join(ntr_df, bigsnp_df)  %>% write_delim(snakemake@output[["match_df"]]
## bigsnp_df <-
## split(allmap,allmap$chrom)
