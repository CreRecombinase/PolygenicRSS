library(EigenH5)
library(LDshrink)
library(tidyverse)
library(SeqArray)
library(progress)


cutoff <- 1e-3
input_file <- snakemake@input[["input_file"]]
mapf <- snakemake@input[["mapf"]]

subsnpf <- snakemake@input[["subsnpf"]]
subldf <- snakemake@input[["subldf"]]
bdf <- snakemake@input[["bdf"]]
output_file <- snakemake@output[["evdf"]]
snpdff  <- snakemake@output[["snpif"]]
useLDshrink <- snakemake@params[["useLDshrink"]]=="T"
pop <- snakemake@params[["pop"]]
r2_cutoff <- snakemake@params[["r2c"]]
stopifnot(!is.null(r2_cutoff))
r2_cutoff <- as.numeric(r2_cutoff)
stopifnot(!is.na(r2_cutoff))

data(map_parameters)
my_map <- filter(map_parameters, Pop==pop)

stopifnot(nrow(my_map)==1)

cutoff <- formals(LDshrink::LDshrink)[["cutoff"]]

m <- my_map$m
Ne <- my_map$Ne



stopifnot(!is.null(input_file), !is.null(output_file), !is.null(mapf), !is.null(bdf),!is.null(snpdff))

stopifnot(file.exists(input_file), !file.exists(output_file), file.exists(mapf), file.exists(bdf))
break_df <- read_delim(bdf, delim="\t")


ind_v <- readRDS(subldf)
snp_df <- read_delim(subsnpf, delim="\t")
op <- nrow(snp_df)
tsnp_df <- snp_df
snp_df <-assign_snp_block(snp_df, break_df, assign_all = T)
semi_join(break_df, snp_df) %>% distinct(chr)
stopifnot(nrow(snp_df)==op)
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




## p <- nrow(snp_df)
## dosage_dims <-EigenH5::dim_h5(input_file, "dosage")
## tch <- EigenH5::dim_h5(input_file, "SNPinfo/chr")


## SNPfirst <-  dosage_dims[1]==tch
## N <- length(ind_v)
## if(!SNPfirst){
##   stopifnot(dosage_dims[2]==tch)
## }

snp_df <- dplyr::mutate(snp_df, ld_snp_id=snp_id)
cat("Writing LDinfo\n")
saveRDS(snp_df, snpdff)
pl <- snakemake@wildcards
pl <- as_data_frame(pl[names(pl)!=""])
#write_df_h5(pl, groupname = "Wildcards", filename=output_file)
snp_dfl <- split(snp_df, snp_df$region_id)
cat("Estimating LD")
num_b <- length(snp_dfl)
pb <- progress::progress_bar$new(total = num_b)
tgds_a <- seqOpen(input_file)
for(i in 1:num_b){
    subsnp_df <- snp_dfl[[i]]
    seqSetFilter(tgds_a,variant.id=subsnp_df$ld_snp_id,sample.sel=ind_v)

    dosage <- seqGetData(tgds_a,"$dosage")

    mrid <- unique(subsnp_df$region_id)
    stopifnot(length(mrid)==1)

    saveRDS(LDshrink_df(dosage,
                                     subsnp_df$map,
                                     subsnp_df$SNP,
                                     m,
                                     Ne,
                                     cutoff,
                                     r2_cutoff,
                                     useLDshrink = useLDshrink,progress=T),output_file)
    pb$tick()
}
