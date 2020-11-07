library(dplyr)
library(ldmap)
library(ldshrink)
library(purrr)
library(EigenH5)
library(BBmisc)
library(duckdb)

shrink <- snakemake@params[["shrink"]] 
ldchunk  <- as.integer(snakemake@params[["ldchunk"]]  %||% 1L)
n_ldchunk <- as.integer(snakemake@params[["n_ldchunks"]]  %||% 1L)
sumstat_dbd <- snakemake@input[["sumststat_dbd"]]
stopifnot(!is.null(sumstat_dbd))
if(is.null(shrink)){
  doshrink <- TRUE
}else{
  doshrink <- shrink=="shrink"
}
snplist_f <- snakemake@input[["snp_list"]]
snplist_df <- tibble(SNP=scan(snplist_f,what=character()))
con  <-  dbConnect(duckdb::duckdb(), dbdir=sumstat_dbd,read_only=TRUE)

sumstat_df <- tbl(con,"gwas")
bim_df <- read_plink_bim(snakemake@input[["bimf"]]) %>% 
  mutate(snp_id = 1:n(),
         ldmr = snp_overlap_region(snp_struct, ldetect_EUR),
         SNP=rsid) %>%
  semi_join(snplist_df) %>%
  semi_join(sumstat_df, by = "SNP", copy = TRUE) %>%
  collect() %>% 
  distinct(SNP, .keep_all = TRUE)
fam_df <- read_plink_fam(snakemake@input[["famf"]])
N <- nrow(fam_df)
bim_l <- chunk(split(bim_df, bim_df$ldmr),n.chunks=n_ldchunk)[[ldchunk]]

purrr::walk(bim_l, function(df){
  cat(as.character(ldetect_EUR[unique(df$ldmr)]),"\n")
  gl <- read_plink_bed(snakemake@input[["bedf"]], subset = df$snp_id, N = N)
  Xm <- gt2matrix(gl)
  cS <- colSums(Xm, na.rm = TRUE)
  cAF <- cS / (N * 2)
  cM <- cS/(N-1)
  bad_snps <- cAF==0 | cAF==1
  indx <- which(is.na(Xm), arr.ind = TRUE)
  Xm[indx] <- cM[indx[,2]]
  sXm <- Xm[,!bad_snps]
  sdf <- filter(df,!bad_snps)
  cor2  <- function(x){1/(NROW(x) -1) * crossprod(scale(x, TRUE , TRUE))}
  if(!doshrink){
    R  <- cor2(sXm)
    ## svdX <- svd(sXm, nu = 0, nv = ncol(sXm))
    ## d <- (svdX$d^2) / (nrow(fam_df)-1)
    ## q <- svdX$v
  }else{
    R <- ldshrink::ldshrink(sXm, sdf$map, isGeno = TRUE, na.rm=FALSE)
  }
  lvdR <- eigen(R)
  d <- lvdR$values
  q <- lvdR$vectors
  ldmr_id <- as.character(ldetect_EUR[unique(sdf$ldmr)])
  write_matrix_h5(q, snakemake@output[["h5f"]], paste0(ldmr_id, "/Q"))
  write_vector_h5(d, snakemake@output[["h5f"]], paste0(ldmr_id, "/D"))
  write_vector_h5(sdf$snp_id, snakemake@output[["h5f"]], paste0(ldmr_id, "/snp_id"))
  write_vector_h5v(sdf$SNP, snakemake@output[["h5f"]], paste0(ldmr_id, "/rsid"))
  write_vector_h5(sdf$snp_struct, snakemake@output[["h5f"]], paste0(ldmr_id, "/snp_struct"))
})
