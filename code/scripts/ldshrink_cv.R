library(EigenH5)
library(LDshrink)
library(tidyverse)


hdff <- snakemake@input[["hdff"]]
mapf <- snakemake@input[["mapf"]]

subsnpf <- snakemake@input[["subsnpf"]]
bdf <- snakemake@input[["bdf"]]
outf <- snakemake@output[["evdf"]]

ldetect <- snakemake@params[["ldetect"]]=="T"
dataset <- snakemake@params[["dataset"]]
snpset <- snakemake@params[["snpset"]]
mapsource <- snakemake@params[["mapsource"]]

m <- as.numeric(snakemake@params[["m"]])
Ne <- as.numeric(snakemake@params[["Ne"]])
cutoff <- as.numeric(snakemake@params[["cutoff"]])
## if(length(cutoff)==0){
##     cutoff <- formals(LDshrink::chunkwise_LDshrink_h5)[["cutoff"]]
## }

stopifnot( !is.null(hdff), !is.null(outf),!is.null(mapf),!is.null(bdf))

stopifnot(file.exists(hdff), !file.exists(outf),file.exists(mapf),file.exists(bdf))
break_df <- read_delim(bdf,delim="\t")



snp_df <- read_delim(subsnpf,delim="\t")
snp_df <-assign_snp_block(snp_df,break_df,assign_all = F) %>% filter(!is.na(region_id))
p <- nrow(snp_df)

stopifnot(sorted_snp_df(snp_df))

map_df <- read_df_h5(mapf,"SNPinfo")
snp_df <- assign_map(snp_df,map_df)
rm(map_df)
stopifnot(group_by(snp_df,chr) %>% summarise(sorted=!is.unsorted(map)) %>% summarise(sorted=all(sorted)) %>% pull(1))


#snp_dfl <- split(snp_df,snp_df$region_id)

cv_ld <- function(sdf,split_l,m,Ne,cutoff){
    train_i <- split_l[["train"]]
    test_i <- split_l[["test"]]

    tH <- EigenH5::read_matrix_h5(hdff,"/","dosage",subset_rows=sdf$snp_id)
    trainH <- t(tH[,train_i])
    testS <- cor(t(tH[,test_i]))
    baseM <- cor(trainH)-testS
    trainS <- calcLD_prel(hmata=trainH,
                          mapa=sdf$map,
                          m=m,
                          Ne=Ne,
                          cutoff=cutoff)
    
    dM <-  trainS-testS
    baseF <- norm(baseM,"F")
    baseE <- max(abs(baseM))
    baseS <- norm(baseM,"2")
    
    dF <- norm(dM,"F")
    dE <- max(abs(dM))
    dS <- norm(dM,"2")
    return(tibble(Frob_LDshrink=dF,
                      maxDist_LDshrink=dE,
                      Spectral_LDshrink=dS,
                      Frob_cor=baseF,
                      maxDist_cor=baseE,
                      Spectral_cor=baseS,
                      p=nrow(sdf)))
}
    
N <- EigenH5::get_dims_h5(hdff,"/","dosage")[2]
ind <- 1:N
train_size <- N/2
train_ind <- sample(ind,train_size,replace=F)
test_ind <- ind[-train_ind]
split_l <- list(train=train_ind,test=test_ind)
cat("Starting")
ret_df <- group_by(snp_df,region_id) %>% do(cv_ld(.,split_l=split_l,m=m,Ne=Ne,cutoff=cutoff)) %>% mutate(snpset=snpset,dataset=dataset,m=m,Ne=Ne,cutoff=cutoff,ldetect=ldetect,mapsource=mapsource)
write_delim(ret_df,outf,delim="\t")

## LDshrink::chunkwise_LDshrink_h5(input_file = hdff,output_file = outf,snp_df = snp_df,evd = T,svd = T,cutoff=cutoff,df=F)
