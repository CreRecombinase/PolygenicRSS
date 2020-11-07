library(dplyr)
vroom::vroom(snakemake@input[["vecf"]]) %>% 
  dplyr::transmute(SNP=SNP,A1=A1,A2=A2,N=N,Z=BETA/SE) %>% 
vroom::vroom_write(snakemake@output[["tempf"]],delim="\t")
