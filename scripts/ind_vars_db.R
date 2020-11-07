library(dplyr)
library(RSQLite)
library(dbplyr)
input_dbf <- snakemake@input[["a_dbf"]]


all_rsids_uniq <- unlist(purrr::map(input_dbf,function(tdbf){

con <- DBI::dbConnect(RSQLite::SQLite(), dbname = tdbf,flags=RSQLite::SQLITE_RO)

  var_db <- tbl(con, "Variant")
  uniq_rsid <- count(var_db,rsid) %>% filter(n==1) %>% select(rsid) %>% collect() %>% pull(1)
  DBI::dbDisconnect(con)
  return(uniq_rsid)
}),use.names=FALSE)

if(sum(duplicated(all_rsids_uniq))>0){
  all_rsids_uniq <- tibble(rsid=all_rsids_uniq) %>% filter(n==1) %>% select(rsid)  %>% pull(1)
}

tibble(rsid=all_rsids_uniq) %>% readr::write_tsv(snakemake@output[["singlef"]], col_names=FALSE)
