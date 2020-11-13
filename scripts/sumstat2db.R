library(duckdb)
library(vroom)
library(dplyr)
library(ldmap)
gwasf <- snakemake@input[["gwasf"]]
output_d <- snakemake@output[["output_d"]]

con = dbConnect(duckdb::duckdb(), dbdir=output_d,read_only=FALSE)
sumstat_df <-  vroom::vroom(gwasf) %>%  
    dplyr::transmute(SNP=SNP,A1=A1,A2=A2,N=N,Z=BETA/SE) %>%
    mutate(rsid = rsid2int(SNP))

dbWriteTable(con,"gwas",sumstat_df)
gc()
dbDisconnect(con,shutdown=TRUE)
gc()
