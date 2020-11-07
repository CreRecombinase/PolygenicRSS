library(RSSp)
library(dplyr)
library(ldmap)
library(ldshrink)
library(purrr)


input_f <- snakemake@input[["qsf"]]
stopifnot(!is.null(input_f))

rssp_df <- map_df(input_f, qs::qread) %>% filter(D > 1e-5)
Nf <- snakemake@input[["samplesize_rdsf"]]
if(!is.null(Nf)){
  N <- qs::qread(Nf)
}else{
  N <- snakemake@params[["samplesize"]]
  stopifnot(!is.null(N))
  N <- as.numeric(N)
}


res_df <- RSSp::RSSp_estimate(
                  rssp_df$quh,
                  rssp_df$D,
                  sample_size = N,
                  p = nrow(rssp_df)
                )
saveRDS(res_df, snakemake@output[["est_rdsf"]])
