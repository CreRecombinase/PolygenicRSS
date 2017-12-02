library(tidyverse)
# save.image()
# stop()
summf <- snakemake@input[["summf"]]
params<- snakemake@params

params <- keep(params,nchar(names(params))>0)
stopifnot(all(lengths(params)==length(summf)))

outf <- snakemake@output[["of"]]
indf <-map(params,as.character) %>% as_data_frame() %>%mutate(temp_t_ind=1:n())
indf <- map_df(summf,function(x){
  read.table(x,header=T) %>% nest()
  }) %>%mutate(temp_t_ind=1:n())%>% inner_join(indf) %>% unnest() %>% select(-temp_t_ind)
write_delim(indf,outf,delim="\t")
