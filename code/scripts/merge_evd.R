library(EigenH5)
library(LDshrink)
library(tidyverse)

library(progress)

input_file <- snakemake@input[["evdfl"]]
output_file <- snakemake@output[["evdf"]]

stopifnot(!is.null(input_file),
          all(file.exists(input_file)))






wildcard_df <-map_df(input_file, ~read_df_h5(.x, "Wildcards")) %>% mutate(input_file=input_file) 

head(wildcard_df)
wc_df <- wildcard_df%>% slice(1)

wildcard_df <- wildcard_df%>% slice(-1)

write_vector_h5(read_vector_h5(wc_df$input_file,"D"),output_file,
                       paste0("EVD/",wc_df$rid,"/D"))
write_matrix_h5(read_matrix_h5(wc_df$input_file[1],"Q"),output_file,
                       paste0("EVD/",wc_df$rid,"/Q"))
write_vector_h5(read_vector_h5(wc_df$input_file[1],"L2"),output_file,
                       paste0("L2/",wc_df$rid,"/L2"))

num_reg <- length(input_file)-1
stopifnot(nrow(wildcard_df)==num_reg)
stopifnot(all(file.exists(wildcard_df$input_file)))
link_objects_h5(filename_from = wildcard_df$input_file,
                filename_to = output_file,
                datapath_from = rep("/D", num_reg),
                datapath_to = paste0("/EVD/", wildcard_df$rid, "/D"))
link_objects_h5(filename_from = wildcard_df$input_file,
                filename_to = output_file,
                datapath_from = rep("/Q", num_reg),
                datapath_to = paste0("/EVD/", wildcard_df$rid, "/Q"))
link_objects_h5(filename_from = wildcard_df$input_file,
                filename_to = output_file,
                datapath_from = rep("/Q", num_reg),
                datapath_to = paste0("/L2/", wildcard_df$rid, "/L2"))

all_snp_df <- map_df(input_file, ~read_df_h5(.x, "LDinfo"))
write_df_h5(all_snp_df, output_file, "LDinfo")
