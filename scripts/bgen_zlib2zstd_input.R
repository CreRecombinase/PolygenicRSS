library(ldmap)
library(stringr)
library(purrr)

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
output_f <- snakemake@output[["ldmap_rs"]]
num_chunks <- length(output_f)
chromosome <- as.integer(snakemake@params[["chrom"]])
num_c  <- ldetect_EUR
chrom_name <- stringr::str_pad(as.character(chromosome),2,pad="0")
lds <- split(ldetect_EUR,chromosomes(ldetect_EUR))

stopifnot(min(lengths(lds))>=num_chunks)

sub_ld <- lds[[chromosome]]
ret_vec <- map_chr(chunk2(sub_ld,num_chunks),function(x){
  cvh <- convex_hull(x)
  paste0(chrom_name,":",starts(cvh),"-",ends(cvh))
})

walk2(output_f,ret_vec,function(outf,rv){
  readr::write_lines(rv,outf)
  })
