library(EigenH5)
library(SeqSupport)
library(tidyverse)
library(progress)
library(SeqArray)


evdf <- snakemake@input[["evdf"]]
traitf <- snakemake@input[["traitf"]]
stopifnot(!is.null(evdf),
          !is.null(traitf))


out_h5f <- snakemake@output[["hdff"]]
stopifnot(!is.null(out_h5f))


#n <- calc_N_h5(c(gdsf,"/","dosage"))

tparam_df <- read_df_h5(traitf,"SimulationInfo")
snp_df <- read_df_h5(evdf,"LDinfo")
EigenH5::write_df_h5(snp_df,"SNPinfo",out_h5f)
p <- nrow(snp_df)
n <- unique(tparam_df$n)
a_grp <- get_objs_h5(evdf,"EVD")
create_matrix_h5(out_h5f,"/","dosage",dims=c(p,n),data=numeric(),chunksizes = c(100,n))
mr <- group_by(snp_df,region_id) %>% do({
  x <-.$region_id[1]
  Q <- read_matrix_h5(evdf,paste0("EVD/",x),"Q")
  D <- read_vector_h5(evdf,paste0("EVD/",x),"D")
  tp <- length(D)
  vm <-matrix(rnorm(tp*n),tp,n)
  tdosage <- evd_rnorm_i(Q = Q,s = sqrt(D),vm = vm)
  offset <- min(.$ld_snp_id)-1
  write_matrix_h5(out_h5f,"/","dosage",tdosage,offsets = c(offset,0))
  data_frame(off=.$ld_snp_id[1])
})

th5 <- read_matrix_h5()

# ymat <- gen_sim_phenotype_h5(snp_df,gdsf,beta_h5file,tparam_df)
# EigenH5::write_matrix_h5(out_h5f,"trait","ymat",ymat)
# EigenH5::write_df_h5(tparam_df,groupname = "SimulationInfo",outfile = out_h5f)

