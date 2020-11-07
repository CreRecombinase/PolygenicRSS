library(ldmap)
library(dplyr)

famf <- snakemake@input[["famf"]]
#grm_id <- snakemake@input[["grm_id"]]
#samplesize <- 12000
samplesize <- as.integer(snakemake@params[["samplesize"]])
read_fam <- function(x,col_names = TRUE){
vroom::vroom(x, col_names = col_names)
}

fam_df <- read_fam(famf,TRUE) %>% rename(fid=ID_1,iid=ID_2) %>% filter(fid!=0)
#ind_df <- read_fam(grm_id,FALSE) %>% 
#rename(fid=X1,iid=X2) %>% 
#filter(fid>0) %>% 
                                        #mutate(fid=fid,iid=iid)
in_df <- fam_df %>% sample_n(samplesize, replace = FALSE)

readr::write_tsv(in_df, snakemake@output[["grm_id"]], col_names = FALSE)


anti_join(fam_df, in_df) %>% sample_n(samplesize,replace=FALSE) %>% 
  readr::write_tsv(snakemake@output[["sub_f"]],col_names=FALSE)
