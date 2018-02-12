library(dplyr)
library(SeqSupport)
library(SeqArray)
library(readr)
library(purrr)
library(progress)

rss_rdsf <- snakemake@input[["rdsf"]]
mfgeneid <- as.character(snakemake@params[["fgeneid"]])
outf <- snakemake@output[["ldscf"]]
names(outf) <- mfgeneid
gds <- seqOpen(snakemake@input[["gdsf"]])


tparam_df <- read_df_h5(rss_rdsf,"tparam_df")
nfgeneid <- tparam_df$fgeneid
tparam_df <- mutate(tparam_df,ldscf=outf[as.character(fgeneid)])



N <- unique(tparam_df$n)
snp_id <- read_vec(rss_rdsf,"snp_id")
uh_l <- array_branch(read_2d_mat_h5(rss_rdsf,"/","uh"),margin=2)

          


seqSetFilter(gds,variant.id=snp_id)
snpinfo <- read_SNPinfo_ldsc(gds)

walk2(uh_l,transpose(tparam_df),
      function(Zd,tp_dfl,snpinfo,Nd){
          cat(tp_dfl$ldscf,"\n")
          mutate(snpinfo,Z=Zd,N=Nd) %>%
              write_delim(path=tp_dfl$ldscf,delim="\t")
          },snpinfo=snpinfo,Nd=N)
#Rprof(NULL)
