library(EigenH5)
library(dplyr)


h5f <- snakemake@input[["h5f"]]
ohf <- snakemake@output[["ohf"]]
chrom  <- as.integer(snakemake@params[["chrom"]])

write_df_h5(read_df_h5(h5f,"SampleInfo"),ohf,"SampleInfo")
write_df_h5(read_df_h5(h5f,"Wildcards"),ohf,"Wildcards")


snp_df <- read_df_h5(h5f,"SNPinfo") %>% filter(chr==chrom) %>% mutate(snp_id=1:n())

write_df_h5(snp_df,ohf,"SNPinfo")
write_matrix_h5(read_matrix_h5(h5f,"dosage",subset_rows=snp_df$snp_id),ohf,"dosage")
