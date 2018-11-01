library(tidyverse)
library(EigenH5)

sample_idf <- snakemake@input[["expf"]]
h5f <- snakemake@output[["h5f"]]
sample_f <- snakemake@input[["sample_list"]]

gene_df <- read_delim(sample_idf, delim=" ")


gene_df %>% select(covariate=ID)  %>%  mutate(covariate_id=1:n())  %>% write_df_h5(h5f,"CovarInfo")

sample_df <- read_delim(sample_f, delim = "\t")

sample_cdf <- data_frame(SampleID=colnames(gene_df)[-1])
stopifnot(all(sample_df$SubjectID==sample_cdf$SampleID))

cvrt_matrix  <- select(gene_df, -ID) %>% data.matrix()

stopifnot(all(colnames(cvrt_matrix) %in% sample_df$SubjectID))

cvrt_matrix[, as.character(sample_df$SubjectID)] %>% write_matrix_h5(h5f, "covariates")
