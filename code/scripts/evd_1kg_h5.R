# save.image("evd.RData")
# stop()
# if(file.exists("evd.RData")){
#     file.remove("evd.RData")
# }

#library(profvis)
library(EigenH5)
library(ldshrink)
library(tidyverse)
# library(future)
library(progress)
# plan(sequential)

                                        #stop()



input_file <- snakemake@input[["input_file"]]
mapf <- snakemake@input[["mapf"]]

subsnpf <- snakemake@input[["subsnpf"]]
subldf <- snakemake@input[["subldf"]]
bdf <- snakemake@input[["bdf"]]
output_file <- snakemake@output[["evdf"]]
useLDshrink <- snakemake@params[["useLDshrink"]] == "T"
chunking <- snakemake@params[["useLDetect"]]
stopifnot(!is.null(chunking))

useLDetect <- chunking == "T"
useChunking <- !is.na(as.numeric(chunking))
if (useChunking) {
        chunking <- as.integer(chunking)
}

## cutoff <- formals(ldshrink::ldshrink)[["cutoff"]]
## m <- formals(ldshrink::ldshrink)[["m"]]
## Ne <- formals(ldshrink::ldshrink)[["Ne"]]


stopifnot(!is.null(input_file), !is.null(output_file), !is.null(mapf), !is.null(bdf))

stopifnot(file.exists(input_file), !file.exists(output_file), file.exists(mapf), file.exists(bdf))
break_df <- read_delim(bdf, delim = "\t") %>% mutate(
        chr = as.integer(chr),
        start = as.integer(start),
        stop = as.integer(stop),
        region_id = as.integer(region_id)
)
head(break_df)


if(TRUE){
    break_df <- group_by(break_df, chr) %>%
        mutate(start = if_else(start == min(start), 0L, start),
               stop = if_else(stop == max(stop), .Machine$integer.max, stop)) %>%
        ungroup()
}


stopifnot(is.integer(break_df$chr))

if (!is.null(subsnpf)) {
        if (fs::path_ext(subsnpf) == "h5") {
                snp_df <- read_df_h5(subsnpf, "SNPinfo")
        } else {
                snp_df <- read_tsv(subsnpf) %>% arrange(chr, pos)
        }
} else {
        snp_df <- read_df_h5(input_file, "SNPinfo")
        all_alleles <- outer(c("A", "C", "T", "G"),
                             c("A", "C", "T", "G"),
                             function(x, y) paste(x, y, sep = ","))
        all_alleles <- tibble(allele = c(all_alleles[upper.tri(all_alleles)],
                                             all_alleles[lower.tri(all_alleles)])
                                  )
        snp_df <- semi_join(snp_df, all_alleles)
}
op <- nrow(snp_df)
if (!useChunking) {
        snp_df <- assign_snp_block(snp_df, break_df, assign_all = T)
        semi_join(break_df, snp_df) %>% distinct(chr)
} else {
        snp_df <- chunk_genome(snp_df, chunk_size = chunking)
}

stopifnot(nrow(snp_df) == op)
p <- nrow(snp_df)
## stopifnot(sorted_snp_df(snp_df))
map_df <- read_df_h5(mapf, "SNPinfo") #%>%
#    distinct(chr, map, .keep_all = T) %>% arrange(chr, pos)
cat("Assigning Map\n")
snp_df <- ldshrink::assign_genetic_map(snp_df, map_df, FALSE)
cat("Removing Map\n")
rm(map_df)
stopifnot(!is.unsorted(snp_df$chr))
stopifnot(!is.unsorted(snp_df$snp_id, strictly = T))
stopifnot(group_by(snp_df, chr) %>%
        summarise(sorted = !is.unsorted(map)) %>%
        summarise(sorted = all(sorted)) %>%
        pull(1))

stopifnot(
        file.exists(input_file),
        !file.exists(output_file),
        !is.null(snp_df[["region_id"]]),
        !is.null(snp_df[["map"]])
)
p <- nrow(snp_df)
dosage_dims <- EigenH5::dim_h5(input_file, "dosage")
tch <- EigenH5::dim_h5(input_file, "SNPinfo/chr")
stopifnot(max(snp_df$snp_id)<=tch)
SNPfirst <- dosage_dims[1] == tch
if (!SNPfirst) {
    stopifnot(dosage_dims[2] == tch,)

}
if (!is.null(subldf)) {
        if (fs::path_ext(subldf) == "h5") {
                ind_v <- read_vector_h5(subldf, "SampleInfo/sample_id")
        } else {
                ind_v <- readRDS(subldf)
        }
        N <- length(ind_v)
} else {
        if (SNPfirst) {
                N <- dosage_dims[2]
        } else {
                N <- dosage_dims[1]
        }
        ind_v <- 1:N
}



snp_df <- dplyr::mutate(snp_df, ld_snp_id = 1:n(),L2=rep(NA_real_,n()))

write_df_h5(snp_df, output_file, "LDinfo")
pl <- snakemake@wildcards
pl <- as_tibble(pl[names(pl) != ""])
write_df_h5(pl, filename = output_file, datapath = "Wildcards")

snp_dfl <- split(snp_df, snp_df$region_id)

cat("Estimating LD")
num_b <- length(snp_dfl)
pb <- progress::progress_bar$new(total = num_b)


for (i in 1:num_b) {
        tdf <- snp_dfl[[i]]

        if (SNPfirst) {
            dosage <- EigenH5::read_matrix_h5(
                                   input_file,
                                   "dosage",
                                   subset_rows = tdf$snp_id,
                                   subset_cols = ind_v,
                                   doTranspose = T)
        } else {
            dosage <- EigenH5::read_matrix_h5(
                                   input_file,
                                   "dosage",
                                   subset_cols = tdf$snp_id,
                                   subset_rows = ind_v)
        }
        mrid <- unique(tdf$region_id)
        stopifnot(length(mrid) == 1)
        retl <- ldshrink_evd(
                reference_panel = dosage,
                map = tdf$map,
                useldshrink = useLDshrink,
                na.rm = F
        )
        EigenH5::write_vector_h5(retl$L2, output_file, paste0("LDinfo/L2"),subset=tdf$ld_snp_id)
        stopifnot(length(retl$D) == length(tdf$pos))
        EigenH5::write_vector_h5(retl$D, output_file, paste0("EVD/", mrid, "/D"))

        EigenH5::write_matrix_h5(retl$Q, output_file, paste0("EVD/", mrid, "/Q"))
        #  EigenH5::write_vector_h5(retl$L2,output_file,paste0("L2/", mrid, "/L2"))
        pb$tick()
}
#   },packages=c("EigenH5","ldshrink")))
# },m=m,Ne=Ne,cutoff=cutoff,useldshrink=useldshrink,SNPfirst=SNPfirst,ind_l=ind_v)


# mmret <- SeqSupport::waitr(retl)
# mret <- purrr::map(retl,future::value)
# })
