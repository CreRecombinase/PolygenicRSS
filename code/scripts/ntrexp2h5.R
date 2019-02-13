inf <- dir("/scratch/t.cri.nknoblauch/polyg_scratch/Wrights_summary_hapmap/chr1",full.names=T,pattern="Rdata")
chrom <- 1
snp_inf <- paste0("~/Downloads/PolygenicRSS/code/snakemake_files/",1:22,"_matched_kg.txt.gz")
outf <- "/scratch/t.cri.nknoblauch/polyg_scratch/gwas_uh/chr1-22AF0SNP0N0_ntr_ntr_data_0_sim.h5"
nuhf <- "/scratch/t.cri.nknoblauch/polyg_scratch/gwas_uh/chr1-22AF0SNP0N0_ntr_EUR_data_0_sim.h5"

## inf <- snakemake@input[["inf"]]
snp_inf <- snakemake@input[["snp_inf"]]
outf <- snakemake@output[["outf"]]





library(EigenH5)
library(tidyverse)
library(fs)
snp_df <- map_df(snp_inf,read_delim,progress=T, col_names=T,delim="\t",col_types=cols(chrom="i",snp="c",pos="i",map="d",snp_id="i",ld_snp_id="i"))

gene_df <- tibble(filename=inf,probe_id=path_ext_remove(path_file(inf))) %>% mutate(fgeneid=1:n())

write_df_h5(gene_df,outf,"EXPinfo")
write_df_h5(snp_df,outf,"SNPinfo")

load_exp <- function(filename,idx){
    load(filename)
    return(tstat.int[idx]/1000)
}

p <- as.integer(nrow(snp_df))
g <- as.integer(nrow(gene_df))
create_matrix_h5(outf,"uh",numeric(),dim=c(p,g),chunksizes=c(as.integer(p/10),1L))
gt <- progress::progress_bar$new(total=g)
pwalk(gene_df,function(filename, fgeneid, idx, ...){
    td <- matrix(load_exp(filename,idx),p,1)
    write_matrix_h5(td,outf,"uh",offset=c(0L, fgeneid-1L), datasize=c(p, 1L))
    gt$tick()
},idx=snp_df$snp_id)
