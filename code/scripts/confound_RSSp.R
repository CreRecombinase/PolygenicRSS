## outf <-  "~/Downloads/PolygenicRSS/output/RSSp_snakemake/est_ntr_chr1-22AF0SNP0N0_ntr_EUR_T_ntr_F_RSSp_res_data_0.01_T_0_F_1.txt.gz"
## res_df <- read_delim(outf,delim="\t") %>% slice(1:10)
## varu <- res_df$sigu^2
## save.image("ws.RData")
## stop()
library(RSSp)


                                        #library(SeqSupport)
library(EigenH5)
library(tidyverse)


quh_dff <- snakemake@input[["quhf"]]
samp_dff <- snakemake@input[["samp_f"]]
outf <- snakemake@output[["dff"]]
doConfound <- as.integer(snakemake@params[["doConfound"]] %||% "1")
pvv <- snakemake@params[["pvv_min"]] %||% 0
trait <- snakemake@params[["trait"]] %||% NA_character  %>% as.character()

#y_grp <- snakemake@params[["y_grp"]] %||% "SimulationInfo"



pvv <- as.numeric(pvv)

stopifnot(!is.na(pvv),
          !is.na(as.numeric(pvv)),
          length(pvv)==1)
## doConfound <- doConfound =="T"


stopifnot(all(file.exists(quh_dff)))
pl <- snakemake@wildcards
pl <- as_data_frame(pl[names(pl)!=""]) %>% mutate(ttca=NA)

if("TraitInfo" %in% ls_h5(quh_dff[1])){
    trait_id_vec <- read_vector_h5(quh_dff[1],"TraitInfo/trait")
    stopifnot(sum(trait %in% trait_id_vec)==1)
    trait_id <- which(trait==trait_id_vec)[1]
    read_quh_fun <- function(x,y){
        data_frame(D=read_vector_h5(x,"quh_df/D"),
                   region_id=read_vector_h5(x,"quh_df/region_id"),
                   quh=c(read_matrix_h5(x,"quh_df/quh",subset_cols=trait_id)),
                   tc=y)
    }
}else{
    read_quh_fun <- function(x,y){
        read_df_h5(x,"quh_df") %>% mutate(tc=y)
    }
}

samp_df <- read_delim(samp_dff,delim="\t") %>% mutate(trait=as.character(trait))
quh_df <- imap_dfr(quh_dff,read_quh_fun)
wc_df <- imap_dfr(quh_dff,~read_df_h5(.x,"Wildcards") %>% select(-chrom))  %>% distinct() %>% mutate(ttca=NA) %>% inner_join(samp_df) %>% inner_join(pl) %>% select(-ttca)

p <- sum(quh_df$D)
n <- wc_df$median_N[1]
p_n <- p/n

filter_normD <- function(D_df,normD=0,add_region_id=F){
  if(pvv==0){
    return(D_df)
  }
  stopifnot(!is.null(D_df$D))
  if(is.null(D_df[["norm_D"]])){
    if(is.null(D_df$region_id)){
      if(add_region_id){
        D_df <- dplyr::mutate(D_df,region_id=1)
      }else{
        stop("column `region_id` missing from D_df, (use `add_region_id=TRUE` to use only one block)")
      }
    }
    stopifnot(!is.null(D_df$region_id))
    D_df <- dplyr::group_by(D_df,region_id) %>% mutate(norm_D=D/mean(D)) %>% ungroup()
  }
  stopifnot(!is.null(D_df[["norm_D"]]))
  ret_df <- D_df %>% dplyr::filter(norm_D>=normD) %>%
    dplyr::ungroup()
  stopifnot(nrow(ret_df)>0)
  return(ret_df)
}

quh_df <- group_by(quh_df,tc) %>% do(filter_normD(.,normD=pvv)) %>% ungroup()
wc_df <- mutate(wc_df,trait_id=paste0(consortia,"_",trait))
trait_id <- unique(wc_df$trait_id)
stopifnot(length(trait_id)==1)


rssp_res <- RSSp_estimate(quh=quh_df$quh,
                          D=quh_df$D,
                          sample_size=n,
                          trait_id=trait_id,
                          nterms=doConfound,
                          pve_bounds=c(.Machine$double.eps, 6 - .Machine$double.eps),
                          eigenvalue_cutoff=0) %>% select(-one_of(c("log_params","useGradient","optim","method"))) %>% inner_join(wc_df) #%>% mutate(new_pve=estimate_pve(dvec=quh_df$D,quh=quh_df$quh,cvec=

saveRDS(rssp_res,outf)





## new_pve <- array_branch(quh,2) %>% map2_dbl(varu,function(tquh,cvec){
##     estimate_pve(D,tquh,cvec,n)
##     })

cat("Done!")
