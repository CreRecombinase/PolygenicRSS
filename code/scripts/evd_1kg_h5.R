#Rprof(filename=snakemake@output[["proff"]],append=F)
library(EigenH5)
library(LDshrink)
library(tidyverse)


data("break_df")

hdff <- snakemake@input[["hdff"]]
mapf <- snakemake@input[["mapf"]]
outf <- snakemake@output[["evdf"]]
cutoff <- as.numeric(snakemake@params[["cutoff"]])
if(length(cutoff)==0){
    cutoff <- formals(SeqSupport::chunkwise_LDshrink)[["cutoff"]]
}

stopifnot( !is.null(hdff), !is.null(outf),!is.null(mapf))

stopifnot(file.exists(hdff), !file.exists(outf),file.exists(mapf))

snp_df <- read_df_h5(normalizePath(hdff),"SNPinfo")

stopifnot(sorted_snp_df(snp_df))

break_dfl <- split(break_df,break_df$chr)
break_df <- map_df(break_dfl,function(x){
mutate(x,start=ifelse(start==min(start),0,start),stop=ifelse(stop==max(stop),.Machine$integer.max,stop))
})




ld_r <- LDshrink::set_ld_region(ld_regions = break_df,snp_info = snp_df)
stopifnot(all(!is.na(ld_r)))


map_df <- read_df_h5(mapf,"SNPinfo")
snp_dfl <- split(snp_df,snp_df$chr)
map_dfl <- split(map_df,map_df$chr)

snp_df <- map2_df(map_dfl,snp_dfl,~mutate(.y,map=LDshrink::interpolate_map(.x$map,.x$pos,.y$pos)))


SeqSupport::chunkwise_LDshrink_h5(hdf5_file = hdff, outfile = outf,region_id=ld_r,map=snp_df$map,cutoff=cutoff)
#
# # tgdsf <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/gds/scombined_19.gds"
# gds_t1d <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/gds/t1d_19.gds"
# tt_gds <- seqOpen(gds_t1d)
# snp_id <- 156068:156073
# #seqSetFilter(tt_gds,variant.sel = snp_id)

# tx <- seqBlockApply(gdsfile = tt_gds,var.name = "$dosage",FUN = function(x){
#   sum(is.na(x))
#   },margin = c("by.variant"),as.is="unlist",.progress=T,parallel = 10)
# t_snp_df <- read_SNPinfo_gds(tt_gds) %>% mutate(chr=as.integer(chr))
# nt_snp_df <- inner_join(tsnp_df,t_snp_df,by=c("SNP","chr","pos"))
# seqSetFilter(tt_gds,variant.sel=nt_snp_df$snp_id.x)
# ttX <- seqGetData(tt_gds,"$dosage")
# tgds <- seqOpen(tgdsf)
# si <- seqGetData(tgds,"sample.id")
# np <- calc_p(tgds)
#
# tsnp_df <- dplyr::filter(snp_df,ld_r==831)
# seqSetFilter(tgds,variant.sel = which(ld_r==831))
# tD <- seqGetData(tgds,"$dosage")
# tG <- seqGetData(tgds,"genotype")
# tH <-read_matrix_h5(hdff,"/","dosage",subset_rows=which(ld_r==831))
# tbad_snp <- which(tH<0,arr.ind=T)
# bSNPs <- duplicated(tH,MARGIN=2)
