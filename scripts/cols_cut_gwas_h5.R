library(EigenH5)
library(dplyr)
library(ldmap)
library(vroom)

gwas_files <- snakemake@input[["vecf"]]
h2 <- snakemake@params[["h2"]]
traits <- snakemake@params[["traits"]]
output_h5 <- snakemake@output[["outf"]]
rsidpath <- "SNP"
snppath <- "snp_struct"
z_mat <- "Z"
param_df <- tibble(h2=h2,traits=traits)

read_fun <- vroom::vroom
co <- vroom::cols_only(CHR=vroom::col_character(),
                       SNP=vroom::col_character(),
                       POS=vroom::col_integer(),
                       A1=vroom::col_character(),
                       A2=vroom::col_character(),
                       BETA=vroom::col_double(),
                       SE=vroom::col_double())
oco <- vroom::cols_only(BETA=vroom::col_double(),
                        SE=vroom::col_double())

EigenH5::write_df_h5(param_df,output_h5,"params")
tf <- read_fun(gwas_files[1],delim="\t",col_types = co)
snp_s <- new_ldmap_snp(tf$CHR,tf$POS,tf$A1,tf$A2)
ldr <- ldmap::snp_in_region(snp_s,ldetect_EUR)
rle2o <- EigenH5::rle2offset(ldr)
EigenH5::write_df_h5(rle2o,output_h5,"ldmap_region_offset")
EigenH5::write_vector_h5(snp_s,output_h5,snppath)
EigenH5::write_vector_h5(tf$SNP,output_h5,rsidpath)

EigenH5::create_matrix_h5(filename = output_h5,datapath = z_mat,data=numeric(),dim=c(nrow(tf),length(gwas_files)))
names(gwas_files) <- seq_along(gwas_files)
progressr::with_progress({
  p <- progressr::progressor(along=gwas_files)
  purrr::iwalk(gwas_files,function(tf,offset){
    p()
    in_df <- read_fun(tf,delim="\t",col_types=oco)
    in_mat <- as.matrix(in_df$BETA/in_df$SE)
    EigenH5::write_matrix_h5(in_mat,filename = output_h5,datapath=z_mat,offset=c(0L,as.integer(offset)-1L),datasize=dim(in_mat))
  })
})
