## save.image("vc.RData")
## stop()
## load("vc.RData")
library(EigenH5)
library(tidyverse)

sample_idf <- snakemake@input[["sampleidf"]]
dosagef <- snakemake@input[["dosagef"]]
chrom  <- as.integer(snakemake@params[["chrom"]])
stopifnot(all(!is.na(chrom)))

dosage_df <- data_frame(filename=dosagef,chrom=chrom) %>%
    arrange(chrom)
outf <- snakemake@output[["outf"]]

num_files <- length(dosagef)
stopifnot(all(file.exists(dosagef)))
file_l <- list(filename = dosage_df$filename, datapath = rep("dosage", num_files))
snp_df <- map_df(file_l$filename,read_df_h5,"SNPinfo",subcols=c("MAF","allele","chr","pos"))
                                        #sample_df <- map_df(file_l$filename,read_df_h5,"SampleInfo")


file_p <- transpose(file_l)
concat_mats(outf,"dosage",file_p,margin = 1)
snp_df <- mutate(snp_df,snp_id=1:n(),chr=as.integer(chr),pos=as.integer(pos))
write_df_h5(snp_df,outf,"SNPinfo")

all_1 <- do.call("cbind",map(file_l$filename,~read_matrix_h5(.x,"dosage",subset_rows=13L)))

fi <- read_matrix_h5(outf,"dosage",subset_rows=13L)

sample_df <- read_delim(sample_idf, delim = "\t", col_names = c("SampleName")) %>% mutate(sample_id=1:n())
write_df_h5(sample_df,outf,"SampleInfo")

#write_df_h5(read_df_h5(file_l$filename,"SampleInfo"),outf,"SampleInfo")
stopifnot(all.equal(all_1,fi))
