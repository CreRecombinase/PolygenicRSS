
## my_chrom <- 17
## input_file <- "/scratch/t.cri.nknoblauch/polyg_scratch/vcf/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
## input_tbi  <- paste0(input_file,".tbi")
## subsnpf <- "~/Downloads/PolygenicRSS/code/snakemake_files/17_matched_kg.txt.gz"
## subldf <- "~/Downloads/PolygenicRSS/data/Snakemake_inputs/EUR.samples"
## asnp_df <- read_delim("~/Downloads/PolygenicRSS/data/Snakemake_inputs/ntr_snps.txt",col_names=c("chrom","snp","pos","map","snp_id"),delim=" ",,col_types=cols(chrom="i",
##                                                                                                       snp="c",
##                                                                                                       map="d",pos="i",snp_id="i"))



#library(profvis)
library(EigenH5)
library(LDshrink)
library(tidyverse)
# library(future)
library(progress)
# plan(sequential)
                                        # mb <-profvis({



read_LDetect <- function(subset=c("all",paste0("chr",1:22)),pop="EUR"){
  if(length(subset)>1){
    warning("In function read_LDetect: length(subset)>1, taking first element of `subset`\n",call.=F)
    subset <- subset[1]
  }
  base_url <- "https://bitbucket.org/nygcresearch/ldetect-data/raw/ac125e47bf7ff3e90be31f278a7b6a61daaba0dc/"
  dl_url <- paste0(base_url,pop,"/fourier_ls-",subset,".bed")
  quietly_read <- purrr::quietly(read_delim) #read_delim is really chatty
  return(quietly_read(dl_url,delim="\t",trim_ws = T) %>% magrittr::extract2("result") %>%  dplyr::mutate(chrom=as.integer(gsub("chr","",chr))) %>% dplyr::select(chrom,start,stop))
}


                                        # load("ssSNP.RData")
## save.image("evd.RData")
## stop()

#file.remove("evd.RData")
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
ld_ind <- scan(subldf,what=character())








subsnp_df <- readr::read_delim(subsnpf,col_names=T,delim="\t") %>%
    mutate(chrom=as.integer(chrom),pos=as.integer(pos)) %>% mutate(map=if_else(map<0,if_else(snp_id==max(snp_id),
                                                                                             max(map)+.Machine$double.eps,
                                                                                             0),map))
p <- nrow(subsnp_df)

if(useLDetect=="T"){
    break_df <- read_LDetect() %>%   dplyr::mutate(region_id=1:n(),start=as.integer(start),stop=as.integer(stop)) %>% dplyr::filter(chrom %in% my_chrom) %>% group_by(chrom) %>% mutate(start=if_else(start==min(start),0L,start),stop=if_else(stop==max(stop),536870912L,stop)) %>% ungroup()
    snp_df <- LDshrink::assign_snp_block(dplyr::rename(subsnp_df,chr=chrom),dplyr::rename(break_df,chr=chrom),assign_all = T)
break_df <- semi_join(break_df, distinct(snp_df, region_id))
}else{

    chunksize <- as.integer(useLDetect)
    num_chunks <- ceiling(p / chunksize)
    stopifnot(length(unique(subsnp_df$chrom))==1, !is.na(chunksize))
    snp_df <- mutate(subsnp_df, region_id= as.integer(gl(num_chunks, chunksize, length=p))) %>%
        rename(chr=chrom)

    break_df <- group_by(snp_df, region_id) %>% summarise(start=min(pos)-1, stop = max(pos)+1)
}

ldetect_regions <- break_df  %>%
  purrr::pmap(function(chrom,start,stop,region_id){ # for each list apply the `pmap` function, which if passed a dataframe, will be call the function once for each row of a dataframe
    GenomicRanges::GRanges(chrom,IRanges::IRanges(start,stop,name=region_id))
  })
names(ldetect_regions) <-  as.character(break_df$region_id)


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
