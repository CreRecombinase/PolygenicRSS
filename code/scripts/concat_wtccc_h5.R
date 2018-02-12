library(EigenH5)
library(tidyverse)
library(progress)

inf <- snakemake@input[["input_h5"]]
outf <- snakemake@output[["outf"]]

nchunks <- length(inf)
pb <- progress_bar$new(total = nchunks)
common_snpf <- intersect_snpinfo_h5(inf)
all_N <- map_int(inf,~get_dims_h5(.x,groupname="/",dataname="dosage")[2])
tot_N <- sum(all_N)
p <- length(common_snpf[[1]])
create_matrix_h5(filename = outf,
                 groupname = "/",
                 dataname = "dosage",
                 data = integer(),
                 doTranspose = F,
                 dims = c(p,tot_N),
                 chunksizes = c(10,tot_N))
tsnp_df <- read_df_h5(h5filepath = inf[1],groupname = "SNPinfo",filtervec =common_snpf[[inf[1]]])
tsnp_df <- mutate(tsnp_df,snp_id=1:n())
stopifnot(nrow(tsnp_df)==p)
write_df_h5(tsnp_df,"SNPinfo",outf)
cur_N <- 0
for(i in 2:length(inf)){
  tinf <- inf[i]
  tdosage <-read_matrix_h5(filename = tinf,groupname = "/",dataname = "dosage",subset_rows = common_snpf[[tinf]])
  write_matrix_h5(filename = outf,groupname = "/",dataname = "dosage",data = tdosage,doTranspose = F,offsets = c(0L,cur_N))
  cur_N <- cur_N+ncol(tdosage)
  gc()
  pb$tick()
}
