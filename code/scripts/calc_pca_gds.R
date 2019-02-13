library(SeqArray)
library(SNPRelate)
library(EigenH5)
library(dplyr)


gdsf <- snakemake@input[["snpgdsf"]]
covarf <- snakemake@output[["covarf"]]
subsnpf <- snakemake@input[["subsnpf"]]
subldf <- snakemake@input[["subgwasf"]]

cores <- snakemake@threads

ind_v <- readRDS(subldf)
snp_df <- readr::read_delim(subsnpf,delim="\t")
#gds <- seqOpen(gdsf)
#ngdsf <- seqGDS2SNP(gds,"~/Dropbox/scratch/polyg_scratch/snpgds/bd_19.gds",compress.geno="LZ4_RA.fast",compress.annotation = "LZ4_RA.fast")
                                        # tf <- tempfile(seqGDS2SNP())

gds <- snpgdsOpen(gdsf)
sample_info <- tibble::tibble(sample_id=gdsfmt::read.gdsn(index.gdsn(gds,"sample.id")))  %>% dplyr::slice(ind_v)
pcao <- snpgdsPCA(gds,sample.id=sample_info$sample_id,snp.id=snp_df$snp_id,num.thread = cores,eigen.cnt = 20,algorithm="randomized")
vecs <-pcao$eigenvect
sample_info <- tibble::tibble(sample_id=gdsfmt::read.gdsn(index.gdsn(gds,"sample.id")))
EigenH5::write_df_h5(sample_info,covarf,"SampleInfo")
covarinfo <- tibble(covariate_id=as.character(1:ncol(vecs)))
write_matrix_h5(vecs,filename = covarf,"covariates")
write_df_h5(covarinfo,covarf,"CovarInfo")
