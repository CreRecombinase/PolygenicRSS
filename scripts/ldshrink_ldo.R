library(dplyr)
library(ldmap)
library(ldshrink)
library(EigenH5)
shrink <- snakemake@params[["shrink"]]
if(is.null(shrink)){
  doshrink <- TRUE
}else{
  doshrink <- shrink=="shrink"
}
snplist_f <- snakemake@input[["snp_list"]]
snplist_df <- tibble(rsid=rsid2int(scan(snplist_f,what=character())))

bim_df <- read_plink_bim(snakemake@input[["bimf"]]) %>%
  mutate(snp_id = 1:n(),
         ldmr = chromosomes(snp_struct),
         rsid = rsid2int(rsid))  %>% 
         semi_join(snplist_df)
fam_df <- read_plink_fam(snakemake@input[["famf"]])
N <- nrow(fam_df)
bim_l <- split(bim_df, bim_df$ldmr)

purrr::walk(bim_l, function(df){
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
  if(!doshrink){
    svdX <- svd(sXm, nu = 0, nv = ncol(sXm))
    d <- (svdX$d^2) / (nrow(fam_df)-1)
    q <- svdX$v
  }else{
    R <- ldshrink::ldshrink(sXm, sdf$map, isGeno = TRUE, na.rm=FALSE)
    lvdR <- eigen(R)
    d <- lvdR$values
    q <- lvdR$vectors
  }
  ldmr_id <- as.character(unique(sdf$ldmr))

  write_matrix_h5(q,snakemake@output[["h5f"]], paste0(ldmr_id, "/Q"))
  write_vector_h5(d, snakemake@output[["h5f"]], paste0(ldmr_id, "/D"))
  write_vector_h5(sdf$snp_id, snakemake@output[["h5f"]], paste0(ldmr_id, "/snp_id"))
  write_vector_h5(sdf$rsid, snakemake@output[["h5f"]], paste0(ldmr_id, "/rsid"))
})
