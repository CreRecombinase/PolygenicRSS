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
doConfound <- snakemake@params[["doConfound"]] %||% "F"
pvv <- snakemake@params[["pvv_min"]] %||% 0

#y_grp <- snakemake@params[["y_grp"]] %||% "SimulationInfo"


pvv <- as.numeric(pvv)

stopifnot(!is.na(pvv),
          !is.na(as.numeric(pvv)),
          length(pvv)==1)
doConfound <- doConfound =="T"


stopifnot(all(file.exists(quh_dff)))
pl <- snakemake@wildcards
pl <- as_data_frame(pl[names(pl)!=""]) %>% mutate(ttca=NA)

samp_df <- read_delim(samp_dff,delim="\t")
quh_df <- imap_dfr(quh_dff,~read_df_h5(.x,"quh_df") %>% mutate(tc=.y))
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
                          p_n=p_n,
                          trait_id=trait_id,
                          pve_bounds=c(.Machine$double.eps, 6 - .Machine$double.eps),
                          eigenvalue_cutoff=0) %>% select(-log_params,-useGradient,-bias,-optim,-method,-contains("confound")) %>% inner_join(wc_df)



## exec_fn <-  function(quhl,tpdfl,D,region_id){
##     max_sigu <- calc_sigu(1,tpdfl$p/tpdfl$n)[1]
##     tbias  <- tpdfl[["tbias"]] %||% NA_real_
##     tpve  <- tpdfl[["tpve"]] %||% NA_real_
##     tsigu <- tpdfl[["tsigu"]] %||% NA_real_
##     data_frame(fgeneid=tpdfl$fgeneid,
##                region_id=region_id,
##                quh=quhl,
##                D=D,
##                p_n=tpdfl$p/tpdfl$n,
##                tbias=tbias,
##                tpve=tpve,
##                tsigu=tsigu) %>%filter_normD(normD=pvv) %>%
##     RSSp_estimate(pve_bounds=c(.Machine$double.eps, 4 - .Machine$double.eps), doConfound=doConfound, eigenvalue_cutoff=0) %>%
##         inner_join(as_data_frame(tpdfl),by="fgeneid")
## }


## il_a <- array_branch(read_matrix_h5(rdsf, "quh", subset_cols=trait), margin=2)
## il_b <- transpose(tparam_df)
## D <- read_vector_h5(rdsf,"D")
## reg <- read_vector_h5(rdsf,"SNPinfo/region_id")

## ## save.image("all_d.RData")
## ## stop()
## rss_res <- map2_df(il_a, il_b,exec_fn,D=D,
##   region_id=reg)


write_delim(rssp_res,outf,delim="\t")




## estimate_pve <- function(dvec,quh,cvec,N,n_samples=0){

##   num_c <- length(cvec)

##   if(num_c>1){
##     stop("multiple cvec terms not yet implemented")
##   }
##   if(n_samples!=0){
##     stop("sampling based pve estimate not yet implemented")
##   }
##   # lambda_init <- 0


##   tmp <- 1+1/(cvec*dvec)
##   rto <- (quh^2)/((dvec)^2)
##   rto <- dvec*rto
##   rto <- rto/(tmp^2)

##   pve_mean <- sum(1/tmp)+sum(rto)
##   pve_mean <- pve_mean/N


##   return(pve_mean)

## }

## new_pve <- array_branch(quh,2) %>% map2_dbl(varu,function(tquh,cvec){
##     estimate_pve(D,tquh,cvec,n)
##     })

cat("Done!")
