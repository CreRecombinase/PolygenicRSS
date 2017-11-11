library(tidyverse)

summf <- snakemake@input[["summf"]]
pv_shr<- snakemake@params[["pv_shr"]]

outf <- snakemake@output[["dff"]]


summ_df <- map2_df(summf,pv_shr,function(x,y){
    return(
        mutate(
            read_delim(x,delim="\t"),
            pv_shr=y
        ) %>%
        separate(pv_shr,c("pv","shrinkage"),sep=",")
    )
})
write_delim(summ_df,outf,delim="\t")
