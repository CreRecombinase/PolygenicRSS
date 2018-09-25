
library(tidyverse)
library(RSSp)
library(LDshrink)
library(EigenH5)



evdf <- snakemake@input[["evdf"]]
quhf <- snakemake@output[["quhf"]]
uhf <- snakemake@input[["uhf"]]
y_grp <- snakemake@params[["y_grp"]]

if(is.null(y_grp)){
    y_grp  <- "SimulationInfo"
}
stopifnot(!is.null(quhf),
          !is.null(evdf),
          file.exists(evdf),
          file.exists(uhf),
          !file.exists(quhf),
          file.exists(uhf))


snp_df <- EigenH5::read_df_h5(evdf,"LDinfo")
snp_df_u <- EigenH5::read_df_h5(uhf,"SNPinfo",subcols=c("chr","pos","allele"))

stopifnot(nrow(snp_df)==nrow(snp_df_u))
stopifnot(all(snp_df$chr==snp_df_u$chr),
              all(snp_df$pos==snp_df_u$pos),
              all(snp_df$allele==snp_df_u$allele))

tparam_df <- EigenH5::read_df_h5(uhf,y_grp)


p <- nrow(snp_df)
g <- nrow(tparam_df)

uh_d <- get_dims_h5(uhf,"uh")
stopifnot(uh_d[1]==p,
          uh_d[2]==g)


uh_l <- split(1:p,snp_df$region_id) %>% imap(~list(subset_rows=.x,
                                                  filename=uhf,
                                                  datapath="uh",region_id=.y))
q_l <- map(uh_l,~list_modify(.x,subset_rows=NULL,filename=evdf,datapath=paste0("EVD/",.x$region_id,"/Q")))
quh_l <- map(uh_l,~list_modify(.x,filename=quhf,datapath="quh"))


D <- map(q_l,~read_vector_h5(filename = .x$filename,datapath=paste0("EVD/",.x$region_id,"/D"))) %>% as_vector(.type="double")
stopifnot(length(D)==p)


EigenH5::write_df_h5(tparam_df, quhf, y_grp)
pl <- snakemake@wildcards
if(is.null(pl[["simulation"]])){
  pl[["simulation"]] <-"gwas"
}
pl <- as_data_frame(pl[names(pl)!=""])
write_df_h5(pl, filename = quhf, datapath = "Wildcards")


EigenH5::write_vector_h5(D,quhf,"D")
create_matrix_h5(quhf,"quh",numeric(),dims=c(p,g),chunksizes = c(pmin(1024,p),1))

EigenH5::write_df_h5(snp_df,quhf,"SNPinfo")
SeqSupport::crossprod_quh_h5(list(Q=q_l,uh=uh_l,quh=quh_l),TRUE)

cat("Done!\n")
