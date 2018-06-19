library(tidyverse)
# library(EigenH5)
library(GWASTools)


# imputedDosageFile(input.files = c(snpdf,snpif,tf),filename=ngds,chromosome=21,input.type="MaCH",input.dosage=T,file.type="gds")


# sample_idf <- snakemake@input[["sampleidf"]]
                                        #infof <- snakemake@input[["sampleinfo"]]
AFc  <- as.numeric(snakemake@params[["AF"]])
chrom <-as.integer(snakemake@params[["chrom"]])
snpif <- snakemake@input[["snpinfo"]]
snpdf <- snakemake@input[["snpdosage"]]
output_f <- snakemake@output[["dosagef"]]
snp_anno_f <-snakemake@output[["snp_anno_f"]]
scan_anno_f <-snakemake@output[["scan_anno_f"]]


## sample_df <- read_delim(infof, delim = "\t", comment = "#", col_names = T)

# sample_df <- read_delim(sample_idf, delim = "\t", col_names = c("SampleID"))
# N <- nrow(sample_df)
snp_df <- read_delim(snpif, delim = "\t") %>%
  separate(SNP,c("chr","position","extra"),remove = F,convert=T) %>%
    mutate(Genotyped = Genotyped == "Genotyped") %>%
    unite(allele, Al1, Al2, sep = ",") %>%
    select(-LooRsq, -EmpR, -EmpRsq, -Dose1, -Dose2) %>%
    mutate(orig_snp_id = 1:n())

tf <- paste0(tempfile(),".txt.gz")

snp_df %>% write_delim(path = tf,delim="\t")

bad_snp_df <- snp_df %>% filter(MAF <= AFc | !is.na(extra)) %>% select(-extra)

index_cols <- bad_snp_df$orig_snp_id

imputedDosageFile(input.files = c(snpdf,snpif,tf),
                  filename=output_f,
                  chromosome=chrom,input.type="MaCH",
                  input.dosage=T,
                  file.type="gds",
                  snp.annot.filename=snp_anno_f,
                  snp.exclude=index_cols,
                  scan.annot.filename=scan_anno_f)

