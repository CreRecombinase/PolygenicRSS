library(EigenH5)
library(tidyverse)
inf <- snakemake@input[["inf"]]
outf_info <- snakemake@output[["outf_info"]]
outf_mean <- snakemake@output[["outf_mean"]]

snp_df <- read_df_h5(inf,"SNPinfo",subcols = c("SNP","allele","chr","pos")) %>%
    rename(rs=SNP) %>%
    separate(allele,c("A","B")) %>% mutate(rs=paste0("rs",rs),                                                               snp_id=1:n())


p <- nrow(snp_df)
chunks <- 100000
split_v <- rep(1:ceiling(p/chunks), each=chunks, length.out=p)
snp_dfl <- split(snp_df, split_v)


snp_df <- map_df(snp_dfl,function(x){
    cat("Reading dosage\n")
    tdosage <- round(read_matrix_h5(inf,"/","dosage",subset_rows=x$snp_id))

    cat("binding cols\n")
    tsnpdf <- bind_cols(select(x,rs,"A","B"),as_data_frame(tdosage))

    cat("Writing res\n")
    write_delim(tsnpdf,outf_mean,delim=" ",col_names=F,append=T)

    cat("calculating AF\n")
    return(mutate(x,af=apply(tdosage,1,mean)/2))
})
snp_df %>% select(rs,A,B,af,chr,pos) %>% write_delim(outf_info,delim="\t")


