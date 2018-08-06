library(EigenH5)
library(tidyverse)


dosagef <- snakemake@input[["dosagef"]]


dosage_df <- data_frame(filename=dosagef,chrom=chrom) %>% arrange(as.integer(chrom))
outf <- snakemake@output[["outf"]]

num_files <- length(dosagef)
stopifnot(all(file.exists(dosagef)))
file_l <- list(filename = dosage_df$filename, datapath = rep("dosage", num_files))
snp_df <- map_df(file_l$filename,read_df_h5,"SNPinfo",subcols=c("MAF","allele","chr","pos"))
                                        #sample_df <- map_df(file_l$filename,read_df_h5,"SampleInfo")
snp_df <- mutate(snp_df,snp_id=1:n())
write_df_h5(snp_df,"SNPinfo",outf)

file_p <- transpose(file_l)
concat_mats(outf,"dosage",file_p,margin = 1)

all_1 <- do.call("cbind",map(file_l$filename,~read_matrix_h5(.x,"/","dosage",subset_rows=13L)))

fi <- read_matrix_h5(outf,"/","dosage",subset_rows=13L)

stopifnot(all.equal(all_1,fi))
