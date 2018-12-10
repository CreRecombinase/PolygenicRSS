library(tidyverse)
library(EigenH5)
library(RSSp)
library(glue)
library(progress)

input <- snakemake@input[["evdf"]]
traitf <- snakemake@input[["traitf"]]
output <- snakemake@output[["quhf"]]
trait <- as.character(snakemake@params[["trait"]])
tnum <- as.integer(snakemake@threads)

cat(glue("Using {tnum} threads\n"))
snp_df <- EigenH5::read_df_h5(input,"LDinfo",subcols=c("region_id","gwas_snp_id"))
if(!is.null(traitf)){
    if(length(traitf)==1){
        data_snp_df <-  fst::read_fst(traitf,columns=c("beta_hat","se","gwas_snp_id")) %>%
            dplyr::slice(snp_df$gwas_snp_id) %>%
            mutate(region_id=snp_df$region_id) %>% mutate(uh=beta_hat/se)
    }else{
        data_snp_df <- list(uh=do.call(cbind,purrr::map(traitf,~dplyr::slice(fst::read_fst(.x,columns=c("beta_hat","se")),snp_df$gwas_snp_id) %>% mutate(uh=beta_hat/se) %>% pull(uh))),
                            region_id=snp_df$region_id)
    }

}else{
    data_snp_df <- EigenH5::read_df_h5(input,"LDinfo",subcols=c("region_id","gwas_snp_id","beta_hat","se")) %>% mutate(uh=beta_hat/se)
}
#stopifnot(all.equal(data_snp_df$gwas_snp_id,snp_df$gwas_snp_id))

quh_fun <- function(x){
    read_matrix_h5(input,glue::glue("EVD/{x}/Q"))
}

snp_fun <- gensplit_uh_uf(data_snp_df[["uh"]],data_snp_df$region_id,na.impute=T)

regid <- unique(data_snp_df$region_id)
#head(snp_df)
#cat(regid)
pb <- dplyr::progress_estimated(length(regid))
Dvec <- unlist(map(glue("EVD/{regid}/D"),~{
    pb$tick()$print()
    read_vector_h5(input,.x)
}))
quhl <-map_convert_quh(regid,snp_fun,quh_fun)
if(NCOL(quhl[[1]])>1){
    rssp_l <- list(quh = do.call(rbind,quhl),
                   region_id=data_snp_df[["region_id"]],
                   D=Dvec)
    cat("Writing\n")
    write_matrix_h5(rssp_l$quh,filename=output,datapath="quh_df/quh")
}else{
    rssp_l <- list(quh = unlist(quhl),
                   region_id=data_snp_df[["region_id"]],
                   D=Dvec)
    cat("Writing\n")
    write_vector_h5(rssp_l$quh,filename=output,datapath="quh_df/quh")
}

write_vector_h5(rssp_l$region_id,filename=output,datapath="quh_df/region_id")
write_vector_h5(rssp_l$D,filename=output,datapath="quh_df/D")
write_vector_h5(trait,filename=output,datapath="TraitInfo/trait")
#rssp_l <- list(quh=,region_id=data_snp_df$region_id,D=Dvec)
stopifnot(all(!is.na(rssp_l$quh)),all(Dvec>0),
          NROW(rssp_l$quh)==length(Dvec),
          NROW(rssp_l$quh)==length(rssp_l$region_id))

pl <- snakemake@wildcards
pl <- as_data_frame(pl[names(pl)!=""])
write_df_h5(pl, filename = output, datapath = "Wildcards")
