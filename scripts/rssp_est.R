library(RSSp)
library(dplyr)
library(ldmap)
library(ldshrink)
library(EigenH5)
library(purrr)

input_f <- snakemake@input[["gwasf"]]
stopifnot(!is.null(input_f))

rssp_df <- map_df(input_f, qs::qread) %>% filter(D > 1e-5)
Nf <- snakemake@input[["N"]]
if(!is.null(Nf)){
  N <- qs::qread(Nf)
}else{
  N <- snakemake@params[["samplesize"]]
  stopifnot(!is.null(N))
  N <- as.integer(N)
}


res_df <- RSSp::RSSp_estimate(
                  rssp_df$quh,
                  rssp_df$D,
                  sample_size = N,
                  p = nrow(rssp_df)
                )
saveRDS(res_df, snakemake@output[["est_rdsf"]])
