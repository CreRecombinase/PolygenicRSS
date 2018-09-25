library(tidyverse)
library(EigenH5)

sample_idf <- snakemake@input[["expf"]]
h5f <- snakemake@output[["h5f"]]
sample_f <- snakemake@input[["sample_list"]]

gene_df <- read_delim(sample_idf, delim=" ")


gene_df %>% select(covariate=ID)  %>%  mutate(covariate_id=1:n())  %>% write_df_h5(h5f,"CovarInfo")

sample_df <- read_delim(sample_f, delim = "\t", col_names = c("SampleID"))
sample_cdf <- data_frame(SampleID=colnames(gene_df)[-1])
stopifnot(all(sample_df$SampleID==sample_cdf$SampleID))

select(gene_df,-ID) %>% data.matrix() %>% write_matrix_h5(h5f,"covariates")
