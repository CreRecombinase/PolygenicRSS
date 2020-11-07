## save.image("quh.RData")
  ## stop("oops!")

## load("quh.RData")
## sumstatf <- "/run/media/nwknoblauch/Backup412/ukb_gwas/100001_raw.gwas.imputed_v3.male.tsv.bgz"
## input_f <- fs::dir_ls("/run/media/nwknoblauch/BackupData/ld_shrink_polym/",glob="*ind*")

  library(RSSp)
  library(dplyr)
  library(EigenH5)
  library(ldmap)
  library(purrr)


  input_f <- snakemake@input[["h5f"]]
  sumstatf <- snakemake@input[["gwasf"]]
                                          #snplist_f <- snakemake@input[["snp_list"]]
                                          #snplist_df <- tibble(rsid=rsid2int(scan(snplist_f,what=character())))
sumstat_df <- vroom::vroom(sumstatf) %>% transmute(snp_struct=as_ldmap_snp(variant),samplesize=n_complete_samples,Z=tstat)

gc()
sumstat_df <- mutate(sumstat_df,ldmr=ldetect_EUR[snp_in_region(snp_struct,ldetect_EUR)])
sumstat_l <- split(sumstat_df,sumstat_df$ldmr)
panel_size <- mean(sumstat_df$samplesize)

oldf <- snakemake@output[["outputf"]]


  read_snp_h5 <- function(file, ldmr_id) {
    tibble::tibble(snp_struct = read_vector_h5(file, paste0(ldmr_id, "/snp_struct")))
  }

  all_reg_df <- map_dfr(set_names(input_f, input_f), ~tibble(ldmr = ls_h5(.x)), .id = "file")



progressr::with_progress({
  p <- progressr::progressor(along=sumstat_l)
  rssp_df <- map2_dfr(all_reg_df$ldmr,all_reg_df$file,function(tldmr,tf) {
    p()
    t_sumstat_df <- sumstat_l[[tldmr]] %>% 
      mutate(snp_pos=clear_alleles(snp_struct))
    ## mutate(Z=read_matrix_h5(sumstat_h5,"Z",offset=c(sub_offset_df$offset,0L),datasize=c(sub_offset_df$datasize,ncol_Z)))
    ld_gwas_df <- read_snp_h5(tf,tldmr) %>% mutate(snp_pos=clear_alleles(snp_struct)) %>% 
      left_join(t_sumstat_df,by="snp_pos",suffix=c("_ld","_gwas")) %>%
      mutate(ams=if_else(allele_match(snp_struct_ld,snp_struct_gwas)=="perfect_match",1,-1))  %>% 
      mutate(Z=ams*Z)
    while(sum(is.na(ld_gwas_df$Z))>0){
      ld_gwas_df$Z[is.na(ld_gwas_df$Z)] <- sample(ld_gwas_df$Z[!is.na(ld_gwas_df$Z)],sum(is.na(ld_gwas_df$Z)),replace=F)
      stopifnot(sum(is.na(ld_gwas_df$Z))==0)
    }

    ntr <- nrow(ld_gwas_df)

    if (exists_h5(tf, as.character(tldmr), "Q") && dim_h5(tf,paste0(tldmr,"/Q"))[1]==ntr) {
      sid <- read_vector_h5(tf, paste0(as.character(tldmr), "/D"))
      Q <- read_matrix_h5(tf, paste0(as.character(tldmr), "/Q"))
      D <- zapsmall(read_vector_h5(tf, paste0(as.character(tldmr), "/D")))
      Q <- Q[,seq_along(D)]
    } else {
      save.image("quh.RData")
      stop(paste0("dim mismatch ",tf,": ",ntr,"!=",dim_h5(tf,paste0(tldmr,"/Q"))[1]," in ",as.character(tldmr)))
    }
    quh <- RSSp::convert_quh(ld_gwas_df$Z, Q)
    gc()
    if(any(is.na(quh))){
      save.image("quh.RData")
      stopifnot(all(!is.na(quh)))
    }
    tibble(quh = quh, D = zapsmall(D), ldmr = as.character(tldmr))
  }) %>% filter(D > 1e-5)
},enable=TRUE)
res_df <- RSSp::RSSp_estimate(
                  rssp_df$quh,
                  rssp_df$D,
                  sample_size = N,
                  p = nrow(rssp_df)
                )
saveRDS(res_df, oldf)
