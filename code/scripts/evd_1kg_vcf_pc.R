
## my_chrom <- 17
input_file <- "/scratch/t.cri.nknoblauch/polyg_scratch/vcf/EUR/EUR.chr2.vcf.gz"
input_tbi  <- paste0(input_file,".tbi")
c_pref <- "/scratch/t.cri.nknoblauch/polyg_scratch/vcf/chunk_snplist/bc/baso-p_EUR.chr2_00"
subsnpf <- "/scratch/t.cri.nknoblauch/polyg_scratch/vcf/chunk_snplist/bc/baso-p_EUR.chr2_01"
## subldf <- "~/Downloads/PolygenicRSS/data/Snakemake_inputs/EUR.samples"
## asnp_df <- read_delim("~/Downloads/PolygenicRSS/data/Snakemake_inputs/ntr_snps.txt",col_names=c("chrom","snp","pos","map","snp_id"),delim=" ",,col_types=cols(chrom="i",
##                                                                                                       snp="c",
##                                                                                                       map="d",pos="i",snp_id="i"))




library(EigenH5)
library(LDshrink)
library(tidyverse)
library(progress)







cutoff <- 1e-3
input_file <- snakemake@input[["input_file"]]
input_tbi  <- paste0(input_file,".tbi")
my_chrom <- snakemake@params[["chrom"]]

subsnpf <- snakemake@input[["subsnpf"]]
subldf <- snakemake@input[["subldf"]]
useLDetect <- snakemake@params[["useLDetect"]]
output_file <- snakemake@output[["evdf"]]
useLDshrink <- snakemake@params[["useLDshrink"]]=="T"
my_vcf <- Rsamtools::TabixFile(input_file,input_tbi)
#ld_ind <- scan(subldf,what=character())

cn_df <- scan(c_pref,what=character(),sep="\t",nlines=1)
subsnp_df <- read_delim(subsnpf,delim="\t",col_names=cn_df)
mchr <- as.character(unique(as.integer(subsnp_df$chrom)))
reg <- range(subsnp_df$pos)







cutoff <- formals(LDshrink::LDshrink)[["cutoff"]]
m <- formals(LDshrink::LDshrink)[["m"]]
Ne <- formals(LDshrink::LDshrink)[["Ne"]]

snp_df <- dplyr::mutate(snp_df, snp_id=as.integer(snp_id), ld_snp_id = as.integer(ld_snp_id))

write_df_h5(snp_df, output_file, "LDinfo")
pl <- snakemake@wildcards
pl <- as_data_frame(pl[names(pl)!=""])

write_df_h5(pl, filename = output_file, datapath = "Wildcards")


snp_dfl <- split(snp_df, snp_df$region_id)

cat("Estimating LD")
num_b <- length(snp_dfl)
pb <- progress::progress_bar$new(total = num_b)

i <- names(snp_dfl)[[1]]
for(i in names(snp_dfl)){

    tdf <- snp_dfl[[i]]
    tldetect_1 <- ldetect_regions[[i]]

    params <- VariantAnnotation::ScanVcfParam(samples=ld_ind,which=unlist(tldetect_1))
    vcf <- VariantAnnotation::readVcf(my_vcf,"hg19",param=params)
    ## vcf_df <- SummarizedExperiment::rowRanges(vcf) %>% as.data.frame(row.names=NULL) %>%
    ##     as_data_frame() %>%
    ##     mutate(chr=as.integer(as.character(seqnames)),
    ##            pos=start,
    ##            snp=rownames(vcf),vcf_id=1:n()) %>%
    ##     dplyr::select(chr,snp,pos,vcf_id)

    ## vcf_idx <- semi_join(vcf_df,tdf,by="pos") %>% pull(vcf_id)
    vcf_idx <- tdf$snp
    vcf <- vcf[vcf_idx,]
    CEU_mat <- VariantAnnotation::genotypeToSnpMatrix(vcf)# First convert to `snpStats` genotype matrix
    dosage <- as(CEU_mat$genotypes,"numeric") #convert `snpStats` matrix to numeric type

    mrid <- unique(tdf$region_id)
    stopifnot(length(mrid)==1)
    if(nrow(dosage)==ncol(dosage)){
        dosage <- t(dosage)
    }
  retl <- LDshrink::LDshrink_evd(dosage,
                       tdf$map,
                       m,
                       Ne,
                       cutoff,
                       useLDshrink=useLDshrink,
                       na.rm=F)
  stopifnot(length(retl$D)==length(tdf$pos))
  EigenH5::write_vector_h5(retl$D,output_file,paste0("EVD/", mrid, "/D"))
  EigenH5::write_matrix_h5(retl$Q,output_file,paste0("EVD/", mrid, "/Q"))
  EigenH5::write_vector_h5(retl$L2,output_file,paste0("L2/", mrid, "/L2"))
  pb$tick()
}
#   },packages=c("EigenH5","LDshrink")))
# },m=m,Ne=Ne,cutoff=cutoff,useLDshrink=useLDshrink,SNPfirst=SNPfirst,ind_l=ind_v)


# mmret <- SeqSupport::waitr(retl)
# mret <- purrr::map(retl,future::value)
# })
