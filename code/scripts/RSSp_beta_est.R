library(RSSp)
library(dplyr)
library(SeqSupport)
library(SeqArray)
library(purrr)
library(readr)

rdsf <- snakemake@input[["rdsf"]]
gdsf <- snakemake@input[["test_gdsf"]]
evdf <- snakemake@input[["evdf"]]
estf <- snakemake@input[["dff"]]

outf <- snakemake@output[["beta_hf"]]


LDchunk <- as.integer(snakemake@params[["LDchunk"]])

X <- read_region_id(gdsf,LDchunk,center=T)

quh_mat <- read_2d_mat_h5(rdsf,LDchunk,"quh")
se_mat <- read_2d_mat_h5(rdsf,LDchunk,"se")
Q <- read_2d_mat_h5(evdf,"EVD","Q")
D <- read_vec(evdf,"EVD/D")
rss_res <- read_delim(estf,delim="\t") %>% filter(method=="Confound") %>%arrange(as.integer(fgeneid))

beta_est <- posterior_mean_Beta(
    sigu=rss_res$sigu,
    confound=rss_res$bias,
    dvec=D,
    quh=quh_mat,
    Q=Q,
    se=se_mat)

beta_oracle <- posterior_mean_Beta(
    sigu=rss_res$tsigu,
    confound=rss_res$tbias,
    dvec=D,
    quh=quh_mat,
    Q=Q,
    se=se_mat)
y_est <- X%*%beta_est
y_oracle <- X%*%beta_oracle

write_mat_h5(outf,"","beta_est",beta_est)
write_mat_h5(outf,"","beta_oracle",beta_oracle)
write_mat_h5(outf,"","y_oracle",y_oracle)
write_mat_h5(outf,"","y_est",y_est)
