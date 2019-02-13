## outf <-  "~/Downloads/PolygenicRSS/output/RSSp_snakemake/est_ntr_chr1-22AF0SNP0N0_ntr_EUR_T_ntr_F_RSSp_res_data_0.01_T_0_F_1.txt.gz"
## res_df <- read_delim(outf,delim="\t") %>% slice(1:10)
## varu <- res_df$sigu^2
library(RSSp)
library(dplyr)

                                        #library(SeqSupport)
library(EigenH5)
library(tidyverse)


rdsf <- snakemake@input[["rdsf"]]
outf <- snakemake@output[["dff"]]
doConfound <- snakemake@params[["doConfound"]] %||% "F"
pvv <- snakemake@params[["pvv_min"]] %||% 0
pvv_max <-

y_grp <- snakemake@params[["y_grp"]] %||% "SimulationInfo"


pvv <- as.numeric(pvv)

stopifnot(!is.na(pvv),
          !is.na(as.numeric(pvv)),
          length(pvv)==1)
doConfound <- doConfound =="T"


stopifnot(file.exists(rdsf),
          !is.null(rdsf))
pl <- snakemake@wildcards
pl <- as_tibble(pl[names(pl)!=""]) %>% mutate(ttca=NA)
g <- dim_h5(rdsf,paste0(y_grp,"/",ls_h5(rdsf,y_grp)[1]))
trait_start <- snakemake@params[["trait"]] %||% 1 %>% as.integer()
num_traits <- snakemake@params[["numtraits"]] %||% g %>% as.integer()

trait <- trait_start:min(c(trait_start+num_traits-1),g)

tparam_df <- read_df_h5(rdsf,y_grp,subset=trait) %>% mutate(ttca=NA) %>% inner_join(pl) %>% select(-ttca) %>% mutate(fgeneid=as.character(fgeneid))
tparam_df[["p"]]  <- tparam_df[["p"]]  %||% dim_h5(rdsf,"SNPinfo/pos")
tparam_df[["n"]]  <- tparam_df[["n"]]  %||% snakemake@params[["samp_size"]]




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






exec_fn <-  function(quhl,tpdfl,D,region_id){
    max_sigu <- calc_sigu(1,tpdfl$p/tpdfl$n)[1]
    tbias  <- tpdfl[["tbias"]] %||% NA_real_
    tpve  <- tpdfl[["tpve"]] %||% NA_real_
    tsigu <- tpdfl[["tsigu"]] %||% NA_real_
    tibble(fgeneid=tpdfl$fgeneid,
               region_id=region_id,
               quh=quhl,
               D=D,
               p_n=tpdfl$p/tpdfl$n,
               tbias=tbias,
               tpve=tpve,
               tsigu=tsigu) %>%filter_normD(normD=pvv) %>%
    RSSp_estimate(pve_bounds=c(.Machine$double.eps, 4 - .Machine$double.eps), doConfound=doConfound, eigenvalue_cutoff=0) %>%
        inner_join(as_tibble(tpdfl),by="fgeneid")
}


il_a <- array_branch(read_matrix_h5(rdsf, "quh", subset_cols=trait), margin=2)
il_b <- transpose(tparam_df)
D <- read_vector_h5(rdsf,"D")
reg <- read_vector_h5(rdsf,"SNPinfo/region_id")

## save.image("all_d.RData")
## stop()
rss_res <- map2_df(il_a, il_b,exec_fn,D=D,
  region_id=reg)


write_delim(rss_res,outf,delim="\t")




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
