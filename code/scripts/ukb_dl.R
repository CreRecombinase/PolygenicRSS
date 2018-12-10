all_gwasf <- "UKBB_mf.csv"
trait <- "21021"
## outf <- "/scratch/t.cri.nknoblauch/polyg_scratch/summary_statistics_ukb/21021_summary_statistics.fst"
library(readr)
library(dplyr)
library(magrittr)
library(tidyr)
library(fst)


all_gwasf <- snakemake@input[["input_f"]]
trait <- snakemake@params[["trait"]]
outf <- snakemake@output[["fstf"]]

gwas_df <- read_csv(all_gwasf)

suffix <- "_irnt"
gwas_sub_df <- filter(gwas_df,`Phenotype Code`==paste0(trait,suffix),Sex=="both_sexes") %>% mutate(new_url=gsub("dl=0","dl=1",`Dropbox File`))
if(nrow(gwas_sub_df)==0){
    suffix <- ""
    gwas_sub_df <- filter(gwas_df,`Phenotype Code`==paste0(trait,suffix),Sex=="both_sexes") %>% mutate(new_url=gsub("dl=0","dl=1",`Dropbox File`))
}
stopifnot(nrow(gwas_sub_df)==1)
out_file <- glue::glue("/dev/shm/{trait}{suffix}.gwas.imputed_v3.both_sexes.tsv.bgz")
wget_cmd <- gsub("-O ","-O /dev/shm/",gwas_sub_df[["wget command"]])
cat("wget_cmd: ",wget_cmd,"\n")
system(wget_cmd)

c_spec <- cols(
  variant = col_character(),
  minor_allele = col_character(),
  minor_AF = col_double(),
  low_confidence_variant = col_character(),
  n_complete_samples = col_integer(),
  AC = col_double(),
  ytx = col_double(),
  beta = col_double(),
  se = col_double(),
  tstat = col_double(),
  pval = col_double()
)

readr::read_tsv(gzfile(out_file),col_types=c_spec) %>%
    tidyr::separate(variant,c("chrom","pos","ref_allele","alt_allele"),sep=":",convert=T) %>%
        mutate(ukb_snp_id=1:n(),chrom=as.integer(chrom),pos=as.integer(pos)) %>%
        select(chrom,pos,ref_allele,alt_allele,beta_hat=beta,se,sample_size=n_complete_samples,p_value=pval) %>% mutate(gwas_snp_id=as.integer(1:n())) %>%
        fst::write_fst(path=outf)

file.remove(out_file)
