library(tidyverse)
library(EigenH5)
off_c1_expf <- snakemake@input[["off_expf_c1"]] #"~/Desktop/dbGaP/files/phe000002.v6.FHS_SABRe_project3.expression-data-matrixfmt.c1/FinalFile_Gene_OFF_2446_Adjusted_c1.txt.gz"
off_c2_expf <- snakemake@input[["off_expf_c2"]] #"~/Desktop/dbGaP/files/phe000002.v6.FHS_SABRe_project3.expression-data-matrixfmt.c2/FinalFile_Gene_OFF_2446_Adjusted_c2.txt.gz"
gen3_c1_expf <- snakemake@input[["gen3_expf_c1"]]# "~/Desktop/dbGaP/files/phe000002.v6.FHS_SABRe_project3.expression-data-matrixfmt.c1/Final_Gene_GENIII_3180_c1.txt.gz"
gen3_c2_expf <- snakemake@input[["gen3_expf_c2"]]#~/Desktop/dbGaP/files/phe000002.v6.FHS_SABRe_project3.expression-data-matrixfmt.c2/Final_Gene_GENIII_3180_c2.txt.gz"
exp_if <- snakemake@input[["expif"]] #"~/Desktop/dbGaP/files/phe000002.v6.FHS_SABRe_project3.sample-info.MULTI/phe000002.v6_release_manifest.txt"
snp_if <- snakemake@input[["snpif"]] # ~/Desktop/dbGaP/files/phg000679.v1.FHS_SHARe_Imputed_1000G.sample-info.MULTI/
exp_infof <- snakemake@input[["expinfof"]]

expf <- snakemake@output[["h5f"]]
explf <- snakemake@output[["sample_list"]]

off_c1 <- read_delim(off_c1_expf,delim="\t",comment="#")
off_c2 <- read_delim(off_c2_expf,delim="\t",comment="#")

sample_off <-  rbind(tibble(FileName=c(colnames(off_c1)[-1]),grp="c1",data="OFF"),
                     tibble(FileName=colnames(off_c2)[-1],grp="c2",data="OFF"))

gen3_c1 <- read_delim(gen3_c1_expf,delim="\t",comment="#")
gen3_c2 <- read_delim(gen3_c2_expf,delim="\t",comment="#")

sample_gen3 <- rbind(tibble(FileName=c(colnames(gen3_c1)[-1]),grp="c1",data="GENIII"),
                     tibble(FileName=colnames(gen3_c2)[-1],grp="c2",data="GENIII"))


sample_exp <- rbind(sample_gen3,sample_off)

sample_df_snp <- read_delim(snp_if,delim="\t",comment = '#')
sample_df_exp <- read_delim(exp_if,delim="\t",comment="#")

share_df_exp <- select(sample_df_exp,SubjectID,FileName) %>% inner_join(select(sample_df_snp,SubjectID))

subset_exp <- inner_join(sample_exp,share_df_exp) %>% mutate(sample_id=1:n())


all_exp <- cbind(cbind(gen3_c1,
                 select(gen3_c2,-transcript_cluster_id)),
                 cbind(select(off_c1,-transcript_cluster_id),
                       select(off_c2,-transcript_cluster_id)))[,c("transcript_cluster_id",subset_exp$FileName)]

exp_info_df <- read_delim(exp_infof,delim=",")
exp_info_df <- exp_info_df %>% select(fgeneid=transcript_cluster_id,chr=seqname,start,stop) %>% mutate(chr=gsub("chr","",chr))

n_all_exp <- inner_join(exp_info_df,rename(all_exp,fgeneid=transcript_cluster_id))
#all.equal(n_exp_info_df$fgeneid,all_exp$transcript_cluster_id)

select(n_all_exp,fgeneid,chr,start,stop) %>% mutate(trait_id=1:n()) %>% write_df_h5(expf,"Traitinfo")
select(n_all_exp,-fgeneid,-chr,-start,-stop) %>% data.matrix() %>% write_matrix_h5(filename = expf,"ymat/trait")
write_df_h5(subset_exp,expf,"SampleInfo")
write_delim(subset_exp,explf,delim="\t")













## sample_idf <- snakemake@input[["expf"]]
## exp_if  <- snakemake@input[["expif"]]
## h5f <- snakemake@output[["h5f"]]
## sample_f <- snakemake@output[["sample_list"]]

## gene_df <- read_delim(sample_idf, delim=" ")
## exp_idf <- read_csv(exp_if)

## new_gene_df <- gene_df %>% select(ID)  %>%
##     separate(ID, into=c("fgeneid","symbol")) %>%
##     mutate(trait_id =1:n(), transcript_cluster_id = as.integer(fgeneid))
## new_gene_df <-
##     inner_join(new_gene_df,exp_idf)  %>%
##     arrange(trait_id)
## new_gene_df <- select(new_gene_df,transcript_cluster_id,seqname,start,stop,fgeneid,trait_id)

## #%>%
## #    select(fgeneid, symbol,trait_id, seqname, strand, start, stop, mrna_assignment, pathway, category)

## cn <- colnames(new_gene_df)
## for( i in sample(cn)){
##     cat(paste0("about to write: ",i,"\n"))
##     write_vector_h5(new_gene_df[[i]],h5f,paste0("TraitInfo/",i))
##     gc()
## }

## sample_df <- tibble(sample_id=colnames(gene_df)[-1])
## sample_df %>% write_df_h5(h5f,"SampleInfo")
## sample_df$sample_id %>% write(sample_f, ncolumns=1)

## select(gene_df,-ID) %>% data.matrix() %>% write_matrix_h5(h5f,"ymat/trait")
