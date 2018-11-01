library(tidyverse)
library(EigenH5)

#eqtl_file <- "/gpfs/data/xhe-lab/Fram/eqtl-gene-annot.txt.gz"
eqtl_h5 <- "/scratch/t.cri.nknoblauch/polyg_scratch/eqtl-gene-annot_sub.h5"
expf <- "/scratch/t.cri.nknoblauch/polyg_scratch/trait/fram_exp.h5"
snpif <- "/home/t.cri.nknoblauch/Downloads/PolygenicRSS/data/Snakemake_inputs/chr22AF0.05SNP0N0_fram_fram_trait_snps.h5"
uhf   <- "/scratch/t.cri.nknoblauch/polyg_scratch/gwas_uh/chr22AF0.05SNP0N0_fram_fram_data_0_sim.h5"
                                        #uhf <- "/scratch/t.cri.nknoblauch/polyg_scratch/gwas_uh/chr10AF0.05SNP0N0_fram_fram_data_10_sim.h5"
subset_number <- 10000
subset_number <- as.integer(snakemake@param[["sample_number"]])

eqtl_file <-snakemake@input[["eqtl_genef"]]
expf <-snakemake@input[["exph5f"]]
#expif <-snakemake@input[["exp_infof"]]
snpif <- snakemake@input[["snp_dff"]]
uhf <- snakemake@input[["uhf"]]

#eqtl_df <-read_delim(eqtl_file,delim=",",n_max = 10000)
uh_dims <- dim_h5(uhf,"uh")
exp_df <- read_df_h5(expf,"Traitinfo") %>%mutate(uh_exp_id=1:n(),transcript_cluster_id=as.integer(fgeneid)) %>% rename(chr_exp=chr)

snp_df <- read_df_h5(snpif,"SNPinfo") %>% mutate(uh_snp_id=1:n())  %>% select(-AF,-region_id)

stopifnot(nrow(snp_df)==uh_dims[1],nrow(exp_df)==uh_dims[2])

chr_vec <- read_vector_h5(eqtl_h5,"eQTLinfo/chr")

un_chr <- unique(snp_df$chr)
stopifnot(length(un_chr)==1)

sub_chr <- which(chr_vec%in%un_chr)
sub_sub_chr <- sort(sample(sub_chr,min(c(subset_number,length(sub_chr))),replace=F))
eqtl_df <- read_df_h5(eqtl_h5,"eQTLinfo",subset=sub_sub_chr)  %>% inner_join(snp_df)  %>% inner_join(exp_df)


unique_exp <- sort(unique(eqtl_df$uh_exp_id))
unique_snp <- sort(unique(eqtl_df$uh_snp_id))
uh_df <- read_matrix_h5(uhf,"uh",subset_rows=unique_snp,subset_cols=unique_exp) %>% set_colnames(as.character(unique_exp))  %>% as_data_frame()  %>% mutate(uh_snp_id=unique_snp)  %>% gather("uh_exp_id","uh",-uh_snp_id) %>% mutate(uh_exp_id=as.integer(uh_exp_id))
uhm <- read_matrix_h5(uhf,"uh",subset_rows=unique_snp,subset_cols=unique_exp)
sem  <- read_matrix_h5(uhf,"se",subset_rows=unique_snp,subset_cols=unique_exp)
bhm  <- uhm*sem
uh_df <- uhm*sem %>% set_colnames(as.character(unique_exp))  %>% as_data_frame()  %>% mutate(uh_snp_id=unique_snp)  %>% gather("uh_exp_id","bh",-uh_snp_id) %>% mutate(uh_exp_id=as.integer(uh_exp_id))
