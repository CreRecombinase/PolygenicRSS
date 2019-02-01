library(tidyverse)
library(EigenH5)
library(SeqSupport)
library(glue)
library(RSSp)



evd_f <- glue("/scratch/t.cri.nknoblauch/polyg_scratch/EVD_H5/ukb/chr{1:22}AF0SNP0N0_ukb_EUR_T_EUR_T.h5")

evd_grps <- map(evd_f,~ls_h5(.x,"EVD"))
names(evd_grps) <- evd_f

D_l <- imap(evd_grps,function(grp,file){
    flatten_dbl(map(grp,~read_vector_h5(file,paste0("EVD/",.x,"/D"))))
})

N <- 250000
p <- sum(lengths(D_l))
D <- flatten_dbl(D_l)

tparam_df <- gen_tparamdf_norm(pve=c(0.06,0.08,.1,.3,.5,.7),bias=0,nreps=5,n=N,p=p)
quh <- do.call("cbind",purrr::map(tparam_df$tsigu,~rnorm(n = p,mean = 0,sd = sqrt(.x^2*D^2+D))))


tparam_df <- rename(tparam_df,trait_id=fgeneid)
res_df <- array_branch(quh,2) %>% imap_dfr(~RSSp_estimate(quh=.x,
                                                          D=D,sample_size=N,
                                                          trait_id=as.character(.y),
                                                          pve_bounds=c(.Machine$double.eps, 6 - .Machine$double.eps),
                                                          eigenvalue_cutoff=0))  %>%
    inner_join(tparam_df)

saveRDS(res_df,"~/Downloads/PolygenicRSS/data/direct_sim_1kg.RDS")

## rssp_res <- RSSp:::RSSp_estimate(quh=quh,
##                           D=D,
##                           p_n=p_n,
