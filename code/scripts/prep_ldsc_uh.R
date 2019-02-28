
library(tidyverse)
library(EigenH5)
# library(SeqSupport)
# library(readr)
# library(purrr)
# library(progress)
#rss_rdsf <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/RSSp_genome_gwas_uh/sscombined_wtcc_NoConfoundSmall_sim.h5"
#tparam_df <- read_df_h5(rss_rdsf,"Traitinfo")
#outf <- paste0("/home/nwknoblauch/Desktop/scratch/polyg_scratch/ldsc_sim_gwas_genome/sscombined_sim_",tparam_df$fgeneid,"_uh_NoConfoundSmall.txt")
#mfgeneid <- as.character(tparam_df$fgeneid)
                                        #names(outf) <- mfgeneid
#
## save.image("ph.RData")
## stop()
rss_rdsf <- snakemake@input[["rdsf"]]
mfgeneid <- as.character(snakemake@params[["fgeneid"]])
outf <- snakemake@output[["ldscf"]]
names(outf) <- mfgeneid



tparam_df <- read_df_h5(rss_rdsf[1],"Traitinfo")

N <- unique(tparam_df$n)

## N <- unique(
nfgeneid <- tparam_df$fgeneid
stopifnot(all(nfgeneid==mfgeneid))


snpinfo_l <- map(rss_rdsf,function(x){
    snp_df <- read_df_h5(x,"SNPinfo") %>%
        select(SNP,allele) %>%
        separate(allele,c("A1","A2")) #%>% mutate(SNP=paste0("rs",SNP))
    return(snp_df)
})


names(snpinfo_l) <- map_chr(rss_rdsf,~as.character(read_vector_h5(.x,"SNPinfo/chr",subset=1L)))
names(rss_rdsf) <- names(snpinfo_l)

## asnp_df <- bind_rows(snpinfo_l

for(i in names(snpinfo_l)){

    uhff <- rss_rdsf[i]
    uhm <- read_matrix_h5(uhff,"uh")
    colnames(uhm) <- mfgeneid
    tsnp_df <- snpinfo_l[[i]]
    for(j in mfgeneid){
        ap <- file.exists(outf[j])
        mutate(tsnp_df,Z=uhm[,j],N=N) %>% select(SNP,N,Z,A1,A2) %>% write_delim(outf[j],delim="\t",append= ap)
    }
}








## uh_df <- as_tibble(read_matrix_h5(rss_rdsf,"uh")) %>%
##   magrittr::set_colnames(as.character(tparam_df$fgeneid)) %>%
##   mutate(snp_id=1:n()) %>%
##   gather("fgeneid","uh",-snp_id)

## stopifnot(all(!is.na(uh_df$uh)))


## seqSetFilter(gds,variant.id=snp_id)



## split(uh_df,uh_df$fgeneid) %>% walk(function(x){
##   mutate(snpinfo,Z=x$uh) %>% write_delim(path =outf[x$fgeneid[1]],delim="\t")
## })
