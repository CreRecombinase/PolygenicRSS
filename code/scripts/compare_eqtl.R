library(tidyverse)
library(EigenH5)

eqtl_file <- "/gpfs/data/xhe-lab/Fram/eqtl-gene-annot.txt.gz"
expf <- "/scratch/t.cri.nknoblauch/polyg_scratch/trait/fram_exp.h5"
snpif <- "/home/t.cri.nknoblauch/Downloads/PolygenicRSS/data/Snakemake_inputs/chr10AF0.05SNP0N0_fram_fram_trait_snps.h5"
uhf <- "/scratch/t.cri.nknoblauch/polyg_scratch/gwas_uh/chr10AF0.05SNP0N0_fram_fram_data_10_sim.h5"
eqtl_file <-snakemake@input[["eqtl_genef"]]
expf <-snakemake@input[["exph5f"]]
#expif <-snakemake@input[["exp_infof"]]
snpif <- snakemake@input[["snp_dff"]]
uhf <- snakemake@input[["uhf"]]

eqtl_df <-read_delim(eqtl_file,delim=",",n_max = 10000)
exp_df <- read_df_h5(expf,"TraitInfo") %>%mutate(uh_exp_id=1:n())
uh_dims <- dim_h5(uhf,"uh")
snp_df <- read_df_h5(snpif,"SNPinfo") %>% mutate(uh_snp_id=1:n())
stopifnot(nrow(snp_df)==uh_dims[1],nrow(exp_df)==uh_dims[2])
