# save.image("sim.RData")
# stop()
library(RSSp)


                                        #library(SeqSupport)
library(EigenH5)
library(tidyverse)


quh_dff <- snakemake@input[["quhf"]]
samp_dff <- snakemake@input[["samp_f"]]
outf <- snakemake@output[["dff"]]
nterms <- snakemake@params[["nterms"]] %||% 1
pvv <- snakemake@params[["pvv_min"]] %||% 0
trait <- snakemake@params[["trait"]] %||% NA_character_  %>% as.character()

#y_grp <- snakemake@params[["y_grp"]] %||% "SimulationInfo"



pvv <- as.numeric(pvv)

stopifnot(!is.na(pvv),
          !is.na(as.numeric(pvv)),
          length(pvv)==1)
#doConfound <- doConfound =="T"


stopifnot(all(file.exists(quh_dff)))
pl <- snakemake@wildcards
pl <- as_tibble(pl[names(pl)!=""]) %>% mutate(ttca=NA)

if("TraitInfo" %in% ls_h5(quh_dff[1])){
    trait_id_vec <- read_vector_h5(quh_dff[1],"TraitInfo/trait")
    # stopifnot(sum(trait %in% trait_id_vec)==1)
    # trait_id <- which(trait==trait_id_vec)[1]
    read_quh_fun <- function(x){
        tibble(D=read_vector_h5(x,"quh_df/D"),
                   region_id=read_vector_h5(x,"quh_df/region_id"),
                   quh=read_matrix_h5(x,"quh_df/quh"))
    }
}else{
    if("quh" %in% ls_h5(quh_dff[1])){
        read_quh_fun <- function(x){



            tibble(
              quh=read_matrix_h5(x,"quh"),
              region_id=read_vector_h5(x,"SNPinfo/region_id"),
              D=read_vector_h5(x,"D"))
        }
    }else{
        read_quh_fun <- function(x){
            read_df_h5(x,"quh_df")
        }
    }
}

if(fs::path_ext(samp_dff)=="h5"){
    samp_df <- read_df_h5(samp_dff,"SimulationInfo") %>% mutate(trait=as.character(fgeneid))
}else{
    samp_df <- read_delim(samp_dff,delim="\t") %>% mutate(trait=as.character(trait))
}
if(TRUE){
tf <- tempfile()
concat_mats(newfile = tf,newpath = "t_quh",transpose(list(filename=quh_dff,datapath=rep("quh",length(quh_dff)))))
quh_df <- tibble(quh=read_matrix_h5(tf,"t_quh"),
                 region_id=unlist(map(quh_dff,~read_vector_h5(.x,"SNPinfo/region_id"))),
                 D=unlist(map(quh_dff,~read_vector_h5(.x,"D"))))
}else{
quh_df <- map_dfr(quh_dff,read_quh_fun)
}
wc_df <- imap_dfr(quh_dff,~read_df_h5(.x,"Wildcards") %>%
                    select(-chrom))  %>%
  distinct() %>%
  mutate(ttca=NA) %>%
  inner_join(mutate(samp_df,ttca=NA)) %>% select(-ttca)

p <- sum(quh_df$D)
n <- wc_df[["median_N"]] %||% wc_df[["n"]]
p_n <- p/n



if("consortia" %in% colnames(wc_df)){
  wc_df <- mutate(wc_df,trait_id=paste0(consortia,"_",trait))
}else{

}
trait_id <- unique(wc_df$trait_id)
# stopifnot(length(trait_id)==1)
if("matrix" %in% class(quh_df$quh)){
  rssp_res <-  map2_dfr(array_branch(quh_df$quh,margin = 2),wc_df$trait_id, ~RSSp:::RSSp_estimate(quh=.x,
                          D=quh_df$D,
                          sample_size=n,
                          trait_id=.y,
                          pve_bounds=c(.Machine$double.eps, 6 - .Machine$double.eps),nterms=nterms,
                          eigenvalue_cutoff=0,useGradient=F))%>% inner_join(mutate(wc_df,trait_id=as.character(trait_id)))
}else{
  stopifnot(length(trait_id)==1)
  rssp_res <- RSSp:::RSSp_estimate(quh=quh_df$quh,
                                   D=quh_df$D,
                                   sample_size=n,
                                   trait_id=trait_id,
                                   pve_bounds=c(.Machine$double.eps, 6 - .Machine$double.eps),nterms=nterms,
                                   eigenvalue_cutoff=0,useGradient=F)%>% inner_join(wc_df)

}
#%>% select(-one_of(c("log_params","useGradient","bias","optim","method")),-contains("confound"))  #%>% mutate(new_pve=estimate_pve(dvec=quh_df$D,quh=quh_df$quh,cvec=

rssp_res %>% saveRDS(outf)





## new_pve <- array_branch(quh,2) %>% map2_dbl(varu,function(tquh,cvec){
##     estimate_pve(D,tquh,cvec,n)
##     })

cat("Done!")
