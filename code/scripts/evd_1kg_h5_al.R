## input_file <- "/scratch/t.cri.nknoblauch/polyg_scratch/vcf/EUR/EUR.chr20.vcf.gz"
## input_tbi  <- paste0(input_file, ".tbi")
## #c_pref <- "/scratch/t.cri.nknoblauch/polyg_scratch/vcf/snplist/bc/baso-p_EUR.chr_2_00.txt"
## subsnpf <- "/scratch/t.cri.nknoblauch/polyg_scratch/vcf/snplist/bc/rdw_EUR.chr20.fst"
## chunksize <- 30000

library(EigenH5)
library(hdf5r)
library(ldshrink)
library(dplyr)
library(purrr)
library(tidyr)
library(progress)
library(feather)

input_file <- snakemake@input[["input_file"]]
my_chrom <- snakemake@params[["chrom"]]

subsnpf <- snakemake@input[["subsnpf"]]
subldf <- snakemake@input[["subldf"]]

output_file <- snakemake@output[["evdf"]]
useLDshrink <- snakemake@params[["useLDshrink"]]=="T"

pl <- snakemake@wildcards
pl <- as_data_frame(pl[names(pl)!=""])
write_df_h5(pl, filename = output_file, datapath = "Wildcards")


subsnp_df <- read_feather(subsnpf)


if(nrow(subsnp_df)>0){
    hf <- H5File$new(input_file,"r")
    snprs <- hf[["variants/ID"]][]
    snppos <- hf[["variants/POS"]][]
    snpref <- hf[["variants/REF"]][]

    hid_df <- data_frame(snp=snprs,pos=snppos) %>% mutate(data_snp_id=1:n())
    hd <- hf[["calldata/GT"]]

    if(is.null(subsnp_df[["region_id"]])){
        stop("please assign LD blocks")
    }

## cutoff <- formals(LDshrink::LDshrink)[["cutoff"]]
## m <- formals(LDshrink::LDshrink)[["m"]]
## Ne <- formals(LDshrink::LDshrink)[["Ne"]]

    snp_df <- dplyr::mutate(subsnp_df, snp_id=as.integer(ld_snp_id), ld_snp_id = as.integer(ld_snp_id))

    write_df_h5(snp_df, output_file, "LDinfo")


    snp_dfl <- split(snp_df, snp_df$region_id)


    num_b <- length(snp_dfl)
#    pb <- progress::progress_bar$new(total = num_b)
    i <- 1

for(i in 1:length(snp_dfl)){
    tdf <- snp_dfl[[i]]
    vcf_idx <- inner_join(hid_df,tdf) %>% distinct(snp,pos,.keep_all=T) %>% pull(data_snp_id)
    stopifnot(length(vcf_idx)==nrow(tdf))
    dosage <- hd[,,vcf_idx,drop=F]
    N <- dim(dosage)[2]
    tp <- length(vcf_idx)
    stopifnot(dim(dosage)[3]==tp)
    dim(dosage) <- c(N*2,tp)
    storage.mode(dosage) <- "numeric"
    ## vcf <- vcf[vcf_idx,]
    ## CEU_mat <- VariantAnnotation::genotypeToSnpMatrix(vcf)# First convert to `snpStats` genotype matrix
    ## dosage <- as(CEU_mat$genotypes,"numeric") #convert `snpStats` matrix to numeric type
    dvar <- apply(dosage,2,var)
    if(min(dvar)==0){
        options(tibble.width = Inf)
        print(dplyr::filter(tdf,dvar==0))
        xt <- table(dosage[,dvar==0])
        cat(xt,names(xt),"!!\n")
        stop("non-variant sites")
    }
    stopifnot(all(!is.na(c(dosage))),
              all(dvar>0))

    mrid <- unique(tdf$region_id)
    stopifnot(length(mrid)==1)
    if(nrow(dosage)==ncol(dosage)){
        dosage <- t(dosage)
    }
    cat("Estimating LD")
    retl <- ldshrink::ldshrink_evd(dosage,
                                   tdf$map,
                                   useldshrink=useLDshrink,
                                   na.rm=F,progress=T)
    stopifnot(length(retl$D)==nrow(tdf))
    EigenH5::write_vector_h5(retl$D,output_file,paste0("EVD/", mrid, "/D"))
    EigenH5::write_matrix_h5(retl$Q,output_file,paste0("EVD/", mrid, "/Q"))
    EigenH5::write_vector_h5(retl$L2,output_file,paste0("L2/", mrid, "/L2"))
#    EigenH5::write_df_h5(tdf,output_file,"LDinfo")

}
}
