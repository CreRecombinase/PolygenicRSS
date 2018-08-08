

library(tidyverse)
library(EigenH5)
# library(GWASTools)


#imputedDosageFile(input.files = c(snpdf,snpif,tf),filename=ngds,chromosome=21,input.type="MaCH",input.dosage=T,file.type="gds")

# sample_idf <- "/home/nwknoblauch/Desktop/Fram/phg000679.v1.FHS_SHARe_Imputed_1000G.genotype-imputed-data.c2/all_ind.txt"
# AFc <-0.01
# chrom <-22
# snpif <- "/home/nwknoblauch/Desktop/Fram/phg000679.v1.FHS_SHARe_Imputed_1000G.genotype-imputed-data.c2/imputed-metrics/machout.chr22.info.gz"
# snpdf <- "/home/nwknoblauch/Desktop/Fram/phg000679.v1.FHS_SHARe_Imputed_1000G.genotype-imputed-data.c2/machout.chr22.dose_NPU.gz"
 # output_f <- "/home/nwknoblauch/Desktop/Fram/chr22_NPU.h5"
 # file.remove(output_f)


sample_idf <- snakemake@input[["sampleidf"]]
                                        #infof <- snakemake@input[["sampleinfo"]]
AFc  <- as.numeric(snakemake@params[["AF"]])
chrom <-as.integer(snakemake@params[["chrom"]])
snpif <- snakemake@input[["snpinfo"]]
snpdf <- snakemake@input[["snpdosage"]]
output_f <- snakemake@output[["dosagef"]]
snp_anno_f <-snakemake@output[["snp_anno_f"]]
scan_anno_f <-snakemake@output[["scan_anno_f"]]


sample_i <- scan(sample_idf,what=character())
sample_df <- read_delim(sample_idf, delim = "\t", col_names = c("SampleID"))
N <- nrow(sample_df)
snp_df <- read_delim(snpif, delim = "\t") %>%
    separate(SNP, c("chr", "pos", "more"))  %>%
    mutate(Genotyped = Genotyped == "Genotyped") %>%
    unite(allele, Al1, Al2, sep = ",") %>%
    select(-LooRsq, -EmpR, -EmpRsq, -Dose1, -Dose2) %>%
    mutate(orig_snp_id = 1:n())

# tf <- tempfile()

# snp_df %>% separate(SNP,c("chr","position","extra"),remove = F,convert=T) %>% write_delim(path = tf,delim="\t")
p <- nrow(snp_df)
#x <- scan(snpdf, what = character(), sep = "\n", n = 10)



good_snp_df <- snp_df %>% filter(MAF >= AFc|!is.na(more)) %>% select(-more)

index_cols <- good_snp_df$orig_snp_id

stopifnot(length(index_cols)>1,
          length(sample_i)>1)
EigenH5::mach2h5(dosagefile = snpdf,
                 h5file = output_f,
                 datapath = "dosage",
                 snp_idx = index_cols-1,
                 names = sample_i,
                 p = p,options=list(
                           buffer_size = 50000*6,
                           SNPfirst = F,progress=T))

good_snp_df %>% mutate(Genotyped=as.integer(Genotyped)) %>% mutate(snp_id=1:n()) %>% write_df_h5(groupname = "SNPinfo",filename = output_f)
#
# nmdos <-matrix(scan(snpdf,what=character(),sep="\t",nlines = 1634),nrow = 1634,byrow = T)[,-c(1,2)][,index_cols]
# class(nmdos) <- "numeric"
# tdos <-read_matrix_h5(output_f,"/","dosage")
# testthat::expect_equal(nmdos,tdos)
# # imputedDosageFile(input.files = c(snpdf,snpif,tf),\
# #                   filename=output_f,
#                   chromosome=chrom,input.type="MaCH",
#                   input.dosage=T,
#                   file.type="gds",
#                   snp.annot.filename=snp_anno_f,
#                   snp.exclude=index_cols,
#                   scan.annot.filename=scan_anno_f)

