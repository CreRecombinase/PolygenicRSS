## save.image("ws.RData")
## stop()
#load("ws.RData")
library(readr)
library(fst)
library(dplyr)
library(ldshrink)
library(purrr)

## gwasf <- "/scratch/t.cri.nknoblauch/polyg_scratch/summary_statistics/bc_mchc_summary_statistics.tsv.gz"
## vcff <- "/scratch/t.cri.nknoblauch/polyg_scratch/vcf/EUR/snplist/EUR.chr2.txt"
gwasf <- snakemake@input[["gwasf"]]
vcff <- snakemake@input[["vcff"]]
mapff <- snakemake@input[["mapf"]]
output_f <- snakemake@output[["outputf"]]



mapdf <- EigenH5::read_df_h5(mapff,"SNPinfo")  %>% mutate(tmap=round(map,2)) %>% group_by(chr) %>% distinct(tmap,.keep_all=T) %>% ungroup() %>% select(-tmap)


vcf_cn  <- c("chrom",
             "pos",
             "snp",
             "ref_allele",
             "alt_allele",
             "score",
             "filter",
             "INFO")

vcf_cols <- readr::cols_only("chrom"="i",
                        "pos"="i",
                        "snp"="c",
                        "ref_allele"="c",
                        "alt_allele"="c")
## formals(try_interpolate_map) <- alist(map = , map_pos = ,target_pos,progress=TRUE)
gwas_df <- assign_genetic_map(dplyr::rename(read_fst(gwasf),chr=chrom),mapdf,strict=F) %>% distinct(chr,pos,.keep_all=T) %>%  filter(map>0) %>% dplyr::rename(chrom=chr)

checknwrite <- function(x){
    cat("Head!")
    options(tibble.width=Inf)
    ## mutate(x,omap=order(map),odiff=lead(omap)-omap,cdiff=lead(map)-map,tid=1:n())  %>%
    ##     mutate(leadcdiff=lead(cdiff),lagcdiff=lag(cdiff)) %>%
    ##     filter(odiff!=1)  %>%
    ##     head()  %>% select(snp,chrom,pos,map,contains("cdiff"),contains("id"),omap,odiff)  %>% print()
    if(is.unsorted(x$map,strictly=T)){
        saveRDS(x,"ws.RDS")
    }
    stopifnot(!is.unsorted(x$map,strictly=T))
    write_fst(x,output_f[x$chrom[1]],compress=100)
}


snp_df <- map_df(vcff,read_delim,delim="\t",col_names=vcf_cn,col_types=vcf_cols)

b_df <- inner_join(gwas_df,snp_df) %>%
    mutate(ld_snp_id=1:n())

cat("Fraction of SNPs found :",(nrow(b_df)/nrow(filter(gwas_df,chrom==snp_df$chrom[1]))),"\n")


group_by(b_df,chrom) %>% do(checknwrite(.))
