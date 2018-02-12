# Create a new GDS file
library(SNPRelate)
library(SeqArray)
library(EigenH5)
library(tidyverse)

inf <- dir("/media/nwknoblauch/Data/wtcc_input",full.names=T,pattern="^[a-z0-9.]+h5$")
mapf <- "/media/nwknoblauch/Data/1kg/1000-genomes-genetic-maps/interpolated_hapmap.h5"
## marchf <- "/media/nwknoblauch/Data/RsMergeArch.h5"
gn <- gsub(".+input/(.+).h5","\\1",inf)
outf_snpgds <- paste0("/media/nwknoblauch/Data/wtcc_input/",gn,".gds")
outf_seqgds <- paste0("/media/nwknoblauch/Data/wtcc_input/",gn,"_seq_.gds")

out_pref_plink <- paste0("/media/nwknoblauch/Data/wtcc_input/plink/",gn)
out_pref_vcf <- paste0("/media/nwknoblauch/Data/wtcc_input/",gn,".vcf")

# inf <- "/media/nwknoblauch/Data/wtcc_input/bd.h5"
# outf_snpgds <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/gds/bd_snp_wtcc_geno.gds"
# outf_seqgds <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/gds/bd_seq_wtcc_geno.gds"
# out_pref_plink <- "/home/nwknoblauch/Desktop/scratch/polyg_scratch/plink/bd"
# gn <- "bd"
## #out_pref_plink <- snakemake@params[["out_pref_plink"]]



inf <- snakemake@input[["input_h5"]]
gn <- snakemake@params[["gn"]]
outf_snpgds <- snakemake@output[["outf_snpgds"]]
outf_seqgds <- snakemake@output[["outf_seqgds"]]
out_pref_plink <- snakemake@params[["outf_pref_plink"]]



in_dims <- map(inf,~get_dims_h5(.x,"/","X"))

newfile_l <- map(outf_snpgds,createfn.gds)
#merge_df <- read_df_h5(marchf,"mergeInfo")
## map_df <- read_df_h5(mapf,"SNPinfo")

make_snpdf <- function(xf,gn){
  col_names <- c("labels","major","minor","chr","pos")
  map_dfc(col_names,~set_names(data_frame(c(read_matrix_h5(xf,groupname = "/",dataname = .x))),.x)) %>%
    mutate(labels=abs(as.integer(labels)),
           chr=as.integer(chr),
           pos=as.integer(pos),
           major=R.oo::intToChar(major),
           minor=R.oo::intToChar(minor),
           allele=paste0(major,"/",minor)) %>% rename(SNP=labels) %>% select(-major,-minor) %>% mutate(snp_id=1:n(),trait_id=gn)
}

snp_dfl <- map2(inf,gn,make_snpdf)
ops <- map_int(snp_dfl,nrow)
inter_snp_df <- map(snp_dfl,~select(.x,-snp_id,-trait_id)) %>%
  reduce(semi_join) %>% distinct() %>%
  distinct(chr,pos,.keep_all = T) %>%
  distinct(chr,pos,allele,.keep_all = T) %>%
  arrange(chr,pos)
good_snp_df <-map(snp_dfl,~inner_join(.x,inter_snp_df))


good_ps<- map_int(good_snp_df,nrow)
stopifnot(length(unique(good_ps))==1)
Ns <- map2_int(in_dims,ops,~.x[.x!=.y])
new_dims <- map2(good_ps,Ns,c)
# add a flag
walk(newfile_l,~put.attr.gdsn(.x$root, "FileFormat", "SNP_ARRAY"))
sample.id <- map2(Ns,gn,~paste0(1:.x,"_",.y))
walk2(newfile_l,sample.id,~add.gdsn(.x, "sample.id", .y))
# sample.id <-as.character(paste0(1:N,"_",gn))
# Add variables

walk2(newfile_l,good_snp_df,function(newfile,snp_df){
  snp_df <- arrange(snp_df,chr,pos)
  add.gdsn(newfile, "snp.id", snp_df$SNP)
  add.gdsn(newfile,"snp.rs.id",paste0("rs",snp_df$SNP))
  add.gdsn(newfile, "snp.position", snp_df$pos)
  add.gdsn(newfile, "snp.chromosome", snp_df$chr)
  add.gdsn(newfile, "snp.allele", snp_df$allele)
})

#####################################################################
# Create a snp-by-sample genotype matrix

# Add genotypes

var.geno_l <- map2(newfile_l,new_dims,~add.gdsn(.x, "genotype",
                                               valdim=.y, storage="bit2"))

# Indicate the SNP matrix is snp-by-sample
walk(var.geno_l,~put.attr.gdsn(.x, "snp.order"))

# Write SNPs into the file sample by sample
pwalk(list(var.geno=var.geno_l,snp_df=good_snp_df,infile=inf),function(var.geno,snp_df,infile){
  p <-nrow(snp_df)
  p_chunks <-BBmisc::chunk(1:p,chunk.size = 10000)
  pb <- progress::progress_bar$new(total=length(p_chunks))
  for (i in 1:length(p_chunks)){
    idl <- p_chunks[[i]]
    g <- read_matrix_h5(infile,"/","X",subset_rows = idl)
    write.gdsn(var.geno, g, start=c(idl[1],1), count=c(length(idl),-1))
    pb$tick()
  }
})

# Get a description of chromosome codes
#   allowing to define a new chromosome code, e.g., snpgdsOption(Z=27)
option <- snpgdsOption()
walk(newfile_l,function(newfile){
  var.chr <- index.gdsn(newfile, "snp.chromosome")
  put.attr.gdsn(var.chr, "autosome.start", option$autosome.start)
  put.attr.gdsn(var.chr, "autosome.end", option$autosome.end)
  for (i in 1:length(option$chromosome.code))
  {
    put.attr.gdsn(var.chr, names(option$chromosome.code)[i],
                  option$chromosome.code[[i]])
  }
})

  # # Add your sample annotation
samp.annotl <- map2(gn,Ns,~data.frame(trait=rep(.x,.y)))
walk2(newfile_l,samp.annotl,~add.gdsn(.x, "sample.annot", .y))

# Add your SNP annotation
# snp.annot <- data.frame(pass=rep(T,p))
# add.gdsn(newfile, "snp.annot", snp.annot)
#seqSNP2GDS(gds,out.fn = )
walk2(newfile_l,out_pref_plink,snpgdsGDS2BED,snpfirstdim=T)
#walk2(newfile_l,snpgdsGDS2BED(newfile,out_pref_plink)
walk(newfile_l,closefn.gds)
#
# walk2(outf_snpgds,outf_seqgds,seqSNP2GDS,storage.option = "LZ4_RA.fast")
# merge_f <-"/media/nwknoblauch/Data/wtcc_input/combined_seq.gds"
# seqMerge(outf_seqgds,merge_f,storage.option = "LZ4_RA.fast")
# seqGDS2VCF()
# Close the GDS file

