library(SeqSupport)
library(dplyr)
library(tidyr)
library(SeqArray)
library(purrr)
library(readr)
uhf <- snakemake@input[["uhf"]]
Sf <- snakemake@input[["Sf"]]
gdsf <- snakemake@input[["gdsf"]]

fgeneids <- as.character(snakemake@params[["fgeneid"]])
stopifnot(!is.null(fgeneids))
ldpf <- snakemake@output[["ldpredf"]]
stopifnot(!is.null(ldpf))
names(ldpf) <- fgeneids

gds <- seqOpen(gdsf)
n <- calc_N(gds)
snp_info <- read_SNPinfo_gds(gds,alleles=T,MAF=T,region_id=F,map=F,info=T) %>% separate(allele,c("ref","alt"),sep=",") %>% rename(reffrq=MAF,rs=SNP) %>% mutate(chr=paste0("chr",chr))
seqClose(gds)

snp_id <- data_frame(snp_id=read_vec(uhf,"snp_id"))
stopifnot(all.equal(snp_info$snp_id,snp_id$snp_id))
out_files<-  ldpf[transpose(read_df_h5(uhf,"tparam_df")) %>% map_chr("fgeneid")]
stopifnot(all(!is.na(out_files)))
S<-read_vec(Sf,"S")
uhl <- list()


write_ldpred <- function(uh,fgeneidf,S,n,df){
    mutate(df,effalt=uh*S,pval=2*pt(abs(uh),df=n-1,lower.tail=F)) %>% select(chr,pos,ref,alt,reffrq,info,rs,pval,effalt) %>% write_delim(path=fgeneidf,delim="\t")
}

array_branch(read_2d_mat_h5(uhf,"/","uh"),2) %>% set_names(out_files) %>% iwalk(write_ldpred,S=S,n=n,df=snp_info)


                       








