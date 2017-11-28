                                        #Rprof(filename=snakemake@output[["proff"]],append=F)
library(dplyr)
outf <- snakemake@output[["gdsf"]]
out_genof <- snakemake@output[["geno_gdsf"]]
mapdf <- readRDS(snakemake@input[["mapf"]])
toutf <- snakemake@input[["temp_gds"]]
tout_genof <- snakemake@input[["geno_gdsf"]]
breakf <- snakemake@input[["breakf"]]




SeqSupport::import_panel_data(temp_gds=toutf,
                  map_df=mapdf,
                  output_file = outf,
                  ld_break_file=breakf,
                  overwrite=T)

SeqSupport::import_panel_data(temp_gds=tout_genof,
                  map_df=mapdf,
                  output_file = out_genof,
                  ld_break_file=breakf,
                  overwrite=T)
                                                

#Rprof(NULL)
