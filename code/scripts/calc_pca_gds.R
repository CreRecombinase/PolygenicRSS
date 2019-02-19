library(SeqArray)
library(SNPRelate)
library(EigenH5)
library(dplyr)


gdsf <- snakemake@input[["snpgdsf"]]
covarf <- snakemake@output[["covarf"]]
subsnpf <- snakemake@input[["subsnpf"]]
subldf <- snakemake@input[["subgwasf"]]

cores <- snakemake@threads


#gds <- seqOpen(gdsf)
#ngdsf <- seqGDS2SNP(gds,"~/Dropbox/scratch/polyg_scratch/snpgds/bd_19.gds",compress.geno="LZ4_RA.fast",compress.annotation = "LZ4_RA.fast")
                                        # tf <- tempfile(seqGDS2SNP())

gds <- seqOpen(gdsf)

if(is.null(subldf)){
  ind_v <- 1:SeqSupport::calc_N(gds)
}else{
  ind_v <- readRDS(subldf)
}
N <- length(ind_v)

if(!is.null(subsnpf)){
  snp_df <- readr::read_delim(subsnpf,delim="\t")
}else{
  snp_df <- SeqSupport::read_SNPinfo_gds(gds)
}


sample_info <- tibble::tibble(sample_id=gdsfmt::read.gdsn(index.gdsn(gds,"sample.id")))  %>% dplyr::slice(ind_v)

pcao <- snpgdsPCA(gds,sample.id=sample_info$sample_id,snp.id=snp_df$snp_id,num.thread = cores,eigen.cnt = 20)
# singular_values <- sqrt(pcao$eigenval*(N-1))
vecs <-pcao$eigenvect
vals <- pcao$eigenval[!is.na(pcao$eigenval)]
sample_info <- tibble::tibble(sample_id=gdsfmt::read.gdsn(index.gdsn(gds,"sample.id")))
EigenH5::write_df_h5(sample_info,covarf,"SampleInfo")
covarinfo <- tibble(covariate_id=as.character(1:ncol(vecs)),covariate_weight=vals)
write_matrix_h5(vecs,filename = covarf,"covariates")
write_df_h5(covarinfo,covarf,"CovarInfo")
