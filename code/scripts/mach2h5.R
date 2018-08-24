library(tidyverse)
library(EigenH5)


sample_idf <- snakemake@input[["sampleidf"]]
#AFc  <- as.numeric(snakemake@params[["AF"]])
chrom <-as.integer(snakemake@params[["chrom"]])
snpif <- snakemake@input[["snpinfo"]]
snpdf_a <- snakemake@input[["snpdosage_a"]]
snpdf_b <- snakemake@input[["snpdosage_b"]]
output_f <- snakemake@output[["dosagef"]]
snp_anno_f <-snakemake@output[["snp_anno_f"]]
scan_anno_f <-snakemake@output[["scan_anno_f"]]


sample_i <- scan(sample_idf,what=character())
sample_df <- read_delim(sample_idf, delim = "\t", col_names = c("SampleName")) %>% mutate(sample_id=1:n())
write_df_h5(sample_df,output_f,"SampleInfo")


N <- nrow(sample_df)
snp_df <- read_delim(snpif, delim = "\t") %>%
    separate(SNP, c("chr", "pos", "more"))  %>%
    mutate(Genotyped = Genotyped == "Genotyped") %>%
    unite(allele, Al1, Al2, sep = ",") %>%
    select(-LooRsq, -EmpR, -EmpRsq, -Dose1, -Dose2,-MAF) %>%
    mutate(orig_snp_id = 1:n())

p <- nrow(snp_df)

all_alleles  <- outer(c("A", "C", "T", "G"),c("A", "C", "T", "G"), function(x,y)paste(x, y, sep=","))
all_alleles <- data_frame(allele=c(all_alleles[upper.tri(all_alleles)], all_alleles[lower.tri(all_alleles)]))

good_snp_df <- snp_df %>% filter(!is.na(more)) %>% select(-more)  %>% semi_join(all_alleles)

index_cols <- good_snp_df$orig_snp_id

stopifnot(length(index_cols)>1,
          length(sample_i)>1)
good_snp_df %>% mutate(Genotyped=as.integer(Genotyped)) %>% mutate(snp_id=1:n()) %>% write_df_h5(filename = output_f,"SNPinfo")

n_p <- length(index_cols)
create_dataset_h5(output_f,"dosage", numeric(), list(dims=c(N, n_p)))

EigenH5::mach2h5(dosagefile = snpdf_a,
                 h5file = output_f,
                 datapath = "dosage",
                 snp_idx = index_cols-1,
                 names = sample_i,
                 p = p, options=list(
                           buffer_size = 50000*6,
                           progress=T))

cat("Finished 1/2\n")
EigenH5::mach2h5(dosagefile = snpdf_b,
                 h5file = output_f,
                 datapath = "dosage",
                 snp_idx = index_cols-1,
                 names = sample_i,
                 p = p, options=list(
                           buffer_size = 50000*6,
                           progress=T))
cat("Finished 2/2!\n")
