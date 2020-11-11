library(RSSp)
library(dplyr)
library(EigenH5)
library(ldmap)
library(purrr)


input_f <- snakemake@input[["h5f"]]
sumstat_h5 <- snakemake@input[["gwasf"]]
                                        #snplist_f <- snakemake@input[["snp_list"]]
                                        #snplist_df <- tibble(rsid=rsid2int(scan(snplist_f,what=character())))

offset_df <- read_df_h5(sumstat_h5,"ldmap_region_offset") %>% mutate(ldmr=ldetect_EUR[value])
panel_size <- as.integer(snakemake@params[["samplesize"]] %||% 10000)

oldf <- snakemake@output[["oldf"]]

ncol_Z <- dim_h5(sumstat_h5,"Z")[2]

read_snp_h5 <- function(file, ldmr_id) {
  tibble::tibble(snp_struct = read_vector_h5(file, paste0(ldmr_id, "/snp_struct")))
}

all_reg_df <- map_dfr(set_names(input_f, input_f), ~tibble(ldmr = ls_h5(.x)), .id = "file")


read_r <- function(tldmr,tf) {

  sub_offset_df <- slice(offset_df,which(as.character(ldetect_EUR[offset_df$value])==tldmr))
  sumstat_df <- tibble(snp_struct=read_vector_h5(sumstat_h5,datapath="/snp_struct",offset=sub_offset_df$offset,datasize=sub_offset_df$datasize)) %>% 
    mutate(snp_pos=clear_alleles(snp_struct)) %>%
    mutate(Z=read_matrix_h5(sumstat_h5,"Z",offset=c(sub_offset_df$offset,0L),datasize=c(sub_offset_df$datasize,ncol_Z)))
  ld_gwas_df <- read_snp_h5(tf,tldmr) %>% mutate(snp_pos=clear_alleles(snp_struct)) %>% 
    left_join(sumstat_df,by="snp_pos",suffix=c("_ld","_gwas")) %>%
    mutate(ams=if_else(allele_match(snp_struct_ld,snp_struct_gwas)=="perfect_match",1,-1))  %>% 
    mutate(Z=ams*Z)
  if(sum(is.na(ld_gwas_df$Z))>0){
    ld_gwas_df$Z[is.na(ld_gwas_df$Z)] <- sample(ld_gwas_df$Z[!is.na(ld_gwas_df$Z)],sum(is.na(ld_gwas_df$Z)),replace=F)
    stopifnot(sum(is.na(ld_gwas_df$Z))==0)
  }

  ntr <- nrow(ld_gwas_df)
  cat(as.character(tldmr),"\n")

  if (exists_h5(tf, as.character(tldmr), "Q") && dim_h5(tf,paste0(tldmr,"/Q"))[1]==ntr) {
    sid <- read_vector_h5(tf, paste0(as.character(tldmr), "/D"))
    Q <- read_matrix_h5(tf, paste0(as.character(tldmr), "/Q"))
    D <- zapsmall(read_vector_h5(tf, paste0(as.character(tldmr), "/D")))
    Q <- Q[,seq_along(D)]
    ## R <- Q%*%(D*t(Q))
    ## diag(R) <- 1
  } else {
    save.image("quh.RData")
    stop(paste0("dim mismatch ",tf,": ",ntr,"!=",dim_h5(tf,paste0(tldmr,"/Q"))[1]," in ",as.character(tldmr)))
    ## R <- read_matrix_h5(tf, paste0(as.character(tldmr), "/R"))
    ## stopifnot(all(tdf$ld_id) %in% seq_len(nrow(R)))
    ## mR <- RSSp:::modify_R(R[tdf$ld_id, tdf$ld_id, drop = FALSE], tdf$Z, panel_size)
    ## ldvr <- eigen(mR)
    ## Q <- ldvr$vectors
    ## D <- ldvr$values
  }
  quh <- RSSp::convert_quh(ld_gwas_df$Z, Q)
  gc()
  if(any(is.na(quh))){
    save.image("quh.RData")
    stopifnot(all(!is.na(quh)))
  }
  tibble(quh = quh, D = zapsmall(D), ldmr = as.character(tldmr))
}

rssp_df <- map2_dfr(all_reg_df$ldmr,all_reg_df$file, read_r) %>% filter(D > 1e-5)
for(i in seq_len(ncol(rssp_df$quh))){
  x <- rssp_df$quh[,i]
  tibble(quh=x,D=rssp_df$D,ldmr=rssp_df$ldmr) %>%   qs::qsave(oldf[i])
}
