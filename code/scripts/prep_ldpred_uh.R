library(SeqSupport)
library(dplyr)

uhf <- snakemake@input[["uhf"]]

hdf5f <- snakemake@input[["hdf5f"]]

fgeneids <- as.character(snakemake@params[["fgeneid"]])

snp_info <- read_df_h5(hdf5f,"SNPinfo")


snp_id <- data_frame(snp_id=read_vec(uhf,"snp_id"))

stopifnot(all.equal(snp_info$snp_id,snp_id$snp_id))






