## if(interactive()){
##     load("ws.RData")
## }else{
##     save.image("ws.RData")
##     stop()
## }

library(EigenH5)
library(ldshrink)
library(tidyverse)
library(glue)
library(progress)

input_file <- snakemake@input[["evdfl"]]
output_file <- snakemake@output[["evdf"]]

stopifnot(!is.null(input_file),
          all(file.exists(input_file)))






wildcard_df <-map_df(input_file, ~read_df_h5(.x, "Wildcards")) %>% mutate(input_file=input_file)
wildcard_df <- filter(wildcard_df,map_lgl(input_file,~isObject(.x,"LDinfo")))

head(wildcard_df)
all_snp_df <- map_df(unique(wildcard_df$input_file), ~read_df_h5(.x, "LDinfo"))
write_df_h5(all_snp_df, output_file, "LDinfo")
rm(all_snp_df)
gc()
wc_df <- wildcard_df%>% slice(1)

wildcard_df <- wildcard_df%>% slice(-1)


walk(ls_h5(wc_df$input_file,"EVD"),~write_vector_h5(read_vector_h5(wc_df$input_file,glue("EVD/{.x}/D")),output_file,glue("EVD/{.x}/D")))
walk(ls_h5(wc_df$input_file,"EVD"),~write_matrix_h5(read_matrix_h5(wc_df$input_file,glue("EVD/{.x}/Q")),output_file,glue("EVD/{.x}/Q")))
walk(ls_h5(wc_df$input_file,"L2"),~write_vector_h5(read_vector_h5(wc_df$input_file,glue("L2/{.x}/L2")),output_file,glue("L2/{.x}/L2")))


## write_vector_h5(
##     read_vector_h5(wc_df$input_file,
##                    glue("EVD/{ls_h5(inf,'EVD')}/D",inf=wc_df$input_file)),
##                 output_file,
##                 paste0("EVD/",wc_df$rid,"/D"))

## write_matrix_h5(read_matrix_h5(wc_df$input_file[1],
##                                paste0("EVD/",wc_df$rid,"/Q")),
##                 output_file,
##                 paste0("EVD/",wc_df$rid,"/Q"))

## write_vector_h5(read_vector_h5(wc_df$input_file[1],
##                                paste0("L2/",wc_df$rid,"/L2")),
##                 output_file,
##                 paste0("L2/",wc_df$rid,"/L2"))


num_reg <- length(input_file)-1

stopifnot(all(file.exists(wildcard_df$input_file)))

datapath_fun <- function(filename,prefix,postfix,full_names=F){
    paste0("/",ls_h5(filename,prefix,full_names=full_names),"/",postfix)
}


dp_df <- unnest(mutate(wildcard_df,
                       dp=map(input_file,datapath_fun,prefix="EVD",postfix="D",full_names=T)))
dp_qf <- unnest(mutate(wildcard_df,
                       dp=map(input_file,datapath_fun,prefix="EVD",postfix="Q",full_names=T)))
dp_lf <- unnest(mutate(wildcard_df,
                       dp=map(input_file,datapath_fun,prefix="L2",postfix="L2",full_names=T)))

link_objects_h5(filename_from = dp_df$input_file,
                filename_to = output_file,
                datapath_from = dp_df$dp,
                datapath_to = dp_df$dp)

link_objects_h5(filename_from = dp_qf$input_file,
                filename_to = output_file,
                datapath_from = dp_qf$dp,
                datapath_to = dp_qf$dp)

link_objects_h5(filename_from = dp_lf$input_file,
                filename_to = output_file,
                datapath_from = dp_lf$dp,
                datapath_to = dp_lf$dp)
