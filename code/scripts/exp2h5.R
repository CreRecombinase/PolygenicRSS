library(tidyverse)
library(EigenH5)

sample_idf <- snakemake@input[["expf"]]
exp_if  <- snakemake@input[["expif"]]
h5f <- snakemake@output[["h5f"]]
sample_f <- snakemake@output[["sample_list"]]

gene_df <- read_delim(sample_idf, delim=" ")
exp_idf <- read_csv(exp_if)

new_gene_df <- gene_df %>% select(ID)  %>%
    separate(ID, into=c("fgeneid","symbol")) %>%
    mutate(trait_id =1:n(), transcript_cluster_id = as.integer(fgeneid))
new_gene_df <-
    inner_join(new_gene_df,exp_idf)  %>%
    arrange(trait_id)
new_gene_df <- select(new_gene_df,transcript_cluster_id,seqname,start,stop,fgeneid,trait_id)

#%>%
#    select(fgeneid, symbol,trait_id, seqname, strand, start, stop, mrna_assignment, pathway, category)

cn <- colnames(new_gene_df)
for( i in sample(cn)){
    cat(paste0("about to write: ",i,"\n"))
    write_vector_h5(new_gene_df[[i]],h5f,paste0("TraitInfo/",i))
    gc()
}

sample_df <- data_frame(sample_id=colnames(gene_df)[-1])
sample_df %>% write_df_h5(h5f,"SampleInfo")
sample_df$sample_id %>% write(sample_f, ncolumns=1)

select(gene_df,-ID) %>% data.matrix() %>% write_matrix_h5(h5f,"ymat/trait")
