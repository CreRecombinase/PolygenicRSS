
library(EigenH5)
library(tidyverse)

sample_idf <- snakemake@input[["sampleidf"]]
dosagef <- snakemake@input[["dosagef"]]
chrom  <- as.integer(snakemake@params[["chrom"]])
mat_path <- snakemake@params[["mat_path"]] %||% "dosage"
stopifnot(all(!is.na(chrom)))

dosage_df <- tibble(filename=dosagef,chrom=chrom) %>%
    arrange(chrom)
outf <- snakemake@output[["outf"]]

num_files <- length(dosagef)
stopifnot(all(file.exists(dosagef)))
file_l <- list(filename = dosage_df$filename,
               datapath = rep(mat_path, num_files))
snp_df <- map_df(file_l$filename,
                 read_df_h5,
                 "SNPinfo",
                 subcols=c("MAF","allele","chr","pos","snp_id"))  %>%
    rename(orig_snp_id=snp_id)
                                        #sample_df <- map_df(file_l$filename,read_df_h5,"SampleInfo")


file_p <- transpose(file_l)
concat_mats(outf,mat_path,file_p,margin = 1)
snp_df <- mutate(snp_df,snp_id=1:n(),chr=as.integer(chr),pos=as.integer(pos))
write_df_h5(snp_df,outf,"SNPinfo")

all_1 <- do.call("cbind",map(file_l$filename,~read_matrix_h5(.x,mat_path,subset_rows=13L)))



read_delim(sample_idf, delim = "\t") %>% write_df_h5(outf,"SampleInfo")

fi <- read_matrix_h5(outf,mat_path,subset_rows=13L)

#write_df_h5(read_df_h5(file_l$filename,"SampleInfo"),outf,"SampleInfo")
stopifnot(all.equal(all_1,fi))
