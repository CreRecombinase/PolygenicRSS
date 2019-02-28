# save.image("ldi.RData")
# stop()
library(tidyverse)
library(SeqArray)
library(EigenH5)
library(ldshrink)
library(progress)
# data("break_df")
# evdf <- "/run/media/nwknoblauch/Data/EVD_H5/chr19_bd_bd_F_omni_T_F.h5"
#
# outf <- paste0("/home/nwknoblauch/Desktop/scratch/polyg_scratch/eur_w_ld_chr_sscombined/",1:22,".l2.ldscore.gz")
# soutf <- paste0("/home/nwknoblauch/Desktop/scratch/polyg_scratch/eur_w_ld_chr_sscombined/",1:22,".l2.M_5_50")
evdf <- snakemake@input[["evdf"]]
outf <- snakemake@output[["outf"]]
soutf <- snakemake@output[["soutf"]]


file.create(outf)
file.create(soutf)

tdf <-tibble(CHR=character(),
                 SNP=character(),
                 BP=integer(),
                 CM=numeric(),
                 MAF=numeric(),
                 L2=numeric())
walk(outf,write_delim,x=tdf,delim="\t")
walk(soutf,write,x=0L)
stopifnot(length(soutf)==length(outf),
          length(soutf)==22L)
walk(evdf,function(x){
    snp_df <- read_df_h5(x,"LDinfo") %>% rename(CHR=chr,BP=pos,CM=map)

                                        #select(break_df,chr,region_id) %>%
                                        # inner_join(read_df_h5(evdf,"LDinfo")) %>%
                                        #    %>% select(-snp_id)

    ## if(is.null(snp_df[["MAF"]])){
    ##     stopifnot(!is.null(snp_df[["AF"]]))
    ##     snp_df <- mutate(snp_df,MAF=AF)
    ## }

                                        #pb <- progress_bar$new(total=length(outf))
                                        # split(snp_df,snp_df$CHR) %>% walk(function(df){
                                        #   pb$tick()

    tchr <- as.integer(snp_df$CHR[1])
    ldscoref <- outf[tchr]
    countf <- soutf[tchr]
    snp_df %>% select(CHR,SNP,BP,CM,MAF,L2) %>% readr::write_delim(path=outf[tchr],delim="\t")
    nc <- dplyr::filter(snp_df,MAF>0.05) %>% nrow()
    write(x=nc,file=soutf[tchr])
    stopifnot(all(file.exists(c(ldscoref,countf))))

})
#})


#Rprof(NULL)
