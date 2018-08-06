library(EigenH5)
library(LDshrink)
library(tidyverse)
library(SeqArray)
library(progress)
#library(fst)
## input_file <- "/scratch/midway2/nwknoblauch/polyg_scratch/gds/ALL_CEU_geno.gds"
## mapf <- "/project/xinhe/eQTL/1kg/1000-genomes-genetic-maps/CEU_map.RDS"
## subsnpf <- "../../data/Snakemake_inputs/chr11AF0.05SNP0N0_CEU_CEU_LD.txt.gz"
## subldf <- "../../data/Snakemake_inputs/chr11AF0.05SNP0N0_CEU_CEU_LD.RDS"
## pop <- "CEU"
## bdf <- "../../data/Snakemake_inputs/ldetect_F.txt.gz"
## output_file <- "/project2/xinhe/LD/CEU/LD_DF/chr11AF0.05SNP0N0_CEU_CEU_F_CEU_T_0.01.RDS"
## snpdff <- "/project2/xinhe/LD/CEU/LD_DF_SNPList/chr11AF0.05SNP0N0_CEU_CEU_F_CEU_T_0.01.RDS"
## r2_cutoff <- 0.01
## cutoff <- 1e-3

writefn <- saveRDS


input_file <- snakemake@input[["input_file"]]
mapf <- snakemake@input[["mapf"]]
useLDshrink <- T

subsnpf <- snakemake@input[["subsnpf"]]
subldf <- snakemake@input[["subldf"]]
bdf <- snakemake@input[["bdf"]]
output_file <- snakemake@output[["snpdff"]]

useLDshrink <- snakemake@params[["useLDshrink"]]=="T"
useLDT <- snakemake@params[["useLDetect"]]=="T"
pop <- snakemake@params[["pop"]]
r2_cutoff <- snakemake@params[["r2c"]]
stopifnot(!is.null(r2_cutoff))
r2_cutoff <- as.numeric(r2_cutoff)
stopifnot(!is.na(r2_cutoff))
chunk_a <- as.integer(snakemake@params[["chunk_a"]])
chunk_tot <- as.integer(snakemake@params[["chunk_tot"]])

data(map_parameters)
my_map <- filter(map_parameters, Pop==pop)

stopifnot(nrow(my_map)==1)

cutoff <- formals(LDshrink::LDshrink)[["cutoff"]]

m <- my_map$m
Ne <- my_map$Ne

stopifnot(!is.null(input_file), !is.null(output_file), !is.null(mapf), !is.null(bdf), !useLDT)

stopifnot(file.exists(input_file), !file.exists(output_file), file.exists(mapf), file.exists(bdf))
break_df <- read_delim(bdf, delim="\t")


ind_v <- readRDS(subldf)
snp_df <- read_delim(subsnpf, delim="\t")
op <- nrow(snp_df)
tsnp_df <- snp_df
snp_df <- chunk_genome(snp_df, n_chunk=chunk_tot)
p <- nrow(snp_df)

stopifnot(sorted_snp_df(snp_df))

map_df <- readRDS(mapf)
cat("Assigning Map\n")
snp_df <- assign_map(snp_df, map_df)
cat("Removing Map\n")
rm(map_df)
stopifnot(!is.unsorted(snp_df$chr))
stopifnot(!is.unsorted(snp_df$snp_id, strictly = T))

stopifnot(group_by(snp_df, chr) %>%
            summarise(sorted=!is.unsorted(map)) %>%
            summarise(sorted=all(sorted)) %>%
            pull(1))

stopifnot(file.exists(input_file),
          !file.exists(output_file),
          !is.null(snp_df[["region_id"]]),
          !is.null(snp_df[["map"]]))

snp_df <- dplyr::mutate(snp_df, ld_snp_id=snp_id)
#cat("Writing LDinfo\n")

snp_dfl <- split(snp_df, snp_df$region_id)

                                        #cat("Estimating LD ")
snp_dfa <- snp_dfl[[chunk_a]]
saveRDS(snp_dfa, output_file)
