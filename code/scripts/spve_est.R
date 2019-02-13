library(SeqSupport)
library(EigenH5)
library(tidyverse)
library(RSSp)

estimation_file <- snakemake@input[["estf"]]
evdf  <- snakemake@input[["evdf"]]
quhf  <- snakemake@input[["quhf"]]
uhf  <- snakemake@input[["est_betaf"]]
## true_betaf  <- snakemake@input[["true_betaf"]]

est_spvef <- snakemake@output[["est_spvef"]]
## true_spvef <- snakemake@output[["true_spvef"]]



res_df  <- read_delim(estimation_file,delim="\t") %>% mutate(fgeneid=as.character(fgeneid))



tparam_df <- EigenH5::read_df_h5(uhf,"SimulationInfo")


snp_df <- EigenH5::read_df_h5(quhf,"SNPinfo")
snp_df_u <- EigenH5::read_df_h5(uhf,"SNPinfo",subcols=c("chr","pos","allele"))

stopifnot(nrow(snp_df)==nrow(snp_df_u))
stopifnot(all(snp_df$chr==snp_df_u$chr),
              all(snp_df$pos==snp_df_u$pos),
          all(snp_df$allele==snp_df_u$allele),
          nrow(res_df)==nrow(tparam_df))

tparam_df <- EigenH5::read_df_h5(uhf,"SimulationInfo")
stopifnot(all(res_df$fgeneid==tparam_df$fgeneid))
N  <- unique(tparam_df$n)
stopifnot(length(N)==1)
p <- nrow(snp_df)
g <- nrow(tparam_df)

uh_d <- get_dims_h5(uhf,"uh")
stopifnot(uh_d[1]==p,
          uh_d[2]==g,
          tparam_df$p[1]==p)


## se_l <- split(1:p,snp_df$region_id) %>% imap(~list(subset_rows=.x,
##                                                   filename=uhf,
##                                                   datapath="se",region_id=.y))
uh_l <- split(1:p,snp_df$region_id) %>% imap(~list(subset_rows=.x,
                                                  filename=uhf,
                                                  datapath="uh",region_id=.y))
d_l <- map(uh_l,~list_modify(.x,subset_rows=NULL,filename=evdf,datapath=paste0("EVD/",.x$region_id,"/D")))
quh_l <- map(uh_l,~list_modify(.x,filename=quhf,datapath="quh"))

## true_beta_l  <- map(se_l,~list_modify(.x,filename=true_spvef,datapath="Beta_Mean"))
## est_beta_l  <- map(true_beta_l,~list_modify(.x,filename=est_spvef))
## oracle_beta_l  <- map(se_l,~list_modify(.x,subset_cols=.x$subset_rows,subset_rows=NULL,filename=true_betaf,datapath="Beta"))

## create_matrix_h5(est_spvef,"/","Beta_Mean",numeric(),dims=c(p,g),chunksizes = c(pmin(1024,p),1))
## create_matrix_h5(true_spvef,"/","Beta_Mean",numeric(),dims=c(p,g),chunksizes = c(pmin(1024,p),1))




pl <- snakemake@wildcards
if(is.null(pl[["simulation"]])){
  pl[["simulation"]] <-"gwas"
}
pl <- as_tibble(pl[names(pl)!=""])


est_res  <-est_spve_h5(file_l=list(uh=uh_l,
                                   D=d_l,
                                   quh=quh_l),
                       N=N,sigu=res_df$sigu)
colnames(est_res) <- res_df$fgeneid
stopifnot(nrow(est_res)==length(uh_l))
est_df  <- as_tibble(est_res) %>% mutate(region_id=names(uh_l)) %>% gather(fgeneid,pve_chunk,-region_id)

write_df_h5(est_df,"spve_chunks",est_spvef)
write_df_h5(tparam_df,groupname="SimulationInfo",filename=est_spvef)
write_df_h5(pl,groupname = "Wildcards",filename=est_spvef)



## true_res  <-est_spve_h5(file_l=list(se=se_l,
##                                    uh=uh_l,
##                                    Q=q_l,
##                                    R=r_l,
##                                    D=d_l,
##                                    quh=quh_l,
##                                    Beta=true_beta_l),
##                         N=N,sigu=res_df$tsigu)
## colnames(true_res) <- res_df$fgeneid
## stopifnot(nrow(true_res)==length(se_l))
## true_df  <- as_tibble(true_res) %>% mutate(region_id=names(se_l)) %>% gather(fgeneid,pve_chunk,-region_id)

## write_df_h5(true_df,"spve_chunks",true_spvef)
## write_df_h5(pl,groupname = "Wildcards",filename=true_spvef)
## write_df_h5(tparam_df,groupname="SimulationInfo",filename=true_spvef)
## oracle_res <- true_spve_h5(file_l=list(se=se_l,
##                                    uh=uh_l,
##                                    Q=q_l,
##                                    R=r_l,
##                                    quh=quh_l,
##                                    Beta=oracle_beta_l),
##                            N=N)

## colnames(oracle_res) <- res_df$fgeneid
## stopifnot(nrow(oracle_res)==length(uh_l))
## oracle_df  <- as_tibble(oracle_res) %>% mutate(region_id=names(se_l)) %>% gather(fgeneid,pve_chunk,-region_id)
## write_df_h5(oracle_df,"oracle_spve_chunks",true_spvef)
## write_df_h5(oracle_df,"oracle_spve_chunks",est_spvef)











