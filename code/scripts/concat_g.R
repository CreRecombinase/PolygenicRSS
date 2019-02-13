library(EigenH5)
library(tidyverse)

## save.image()
## stop()

input_f <- snakemake@input[["input_f"]]
idx <- snakemake@params[["idx"]]
idx_pth <- snakemake@params[["idx_pth"]]

concat_mats <- snakemake@params[["concat_mat"]]
keep_mats <- snakemake@params[["keep_mat"]]

concat_vecs <- snakemake@params[["concat_v"]]
keep_vecs <- snakemake@params[["keep_v"]]

keep_df <- snakemake@params[["keep_df"]]
concat_df <- snakemake@params[["concat_df"]]


outf <- snakemake@output[["outf"]]

stopifnot(all(file.exists(input_f)),
          length(input_f)==length(idx),
          !file.exists(outf))





concat_fn <- function(input_files, output_file, datapath, keep=F, chunksizes=10000L){
    if(is.null(datapath)){
        return(F)
    }
    tinput_file <- input_files[1]
    if(keep){
        input_files <- tinput_file
        datapath <- datapath[1]
    }
    if(isGroup(input_files[1], datapath)){
        temp_df <- read_df_h5(tinput_file, datapath, subset=1L)  %>% slice(integer())
        write_df_h5(temp_df, output_file, datapath, chunksizes=chunksizes, max_dims=NA_integer_)
        walk(input_files, ~write_df_h5(read_df_h5(.x, datapath), output_file, datapath, append=T))
        return(T)
    }
    data_dim <- dim_h5(tinput_file, datapath)
    if(length(data_dim)==1){
        tempv <- read_vector_h5(tinput_file, datapath, subset = 1L)[-1]
        write_vector_h5(tempv, output_file, datapath, chunksizes=chunksizes, max_dims=c(NA_integer_))
        walk(input_files, ~write_vector_h5(read_vector_h5(.x, datapath), output_file, datapath, append=T))
        return(T)
    }
    ## if(length(data_dim)==1){
    ##     tempv <- read_vector_h5(tinput_file, datapath, subset = 1L)[-1]
    ##     write_vector_h5(tv, output_file, datapath, chunksizes=chunksizes, max_dims=c(NA_integer_))
    ##     walk(input_files, ~write_vector_h5(read_vector_h5(.x, datapath), output_file, datapath, append=T))
    ##     return(T)
    ## }
    if(length(data_dim)!=2){
        stop("can only concatenate dataframes,vectors or matrices")
    }
    file_l <- transpose(list(filename = input_files,
                             datapath = rep(datapath, length(input_files))))

    dim_vec <- map(input_files, ~dim_h5(.x, datapath))
    if(length(unique(map_int(dim_vec, 1))==1)){
        margin <- 0
    }else{
        if(length(unique(map_int(dim_vec, 2))!=1)){
            stop("there must be a common dimension")
        }
        margin <- 1
    }
    concat_mats(output_file, datapath, file_l, margin=margin)
    return(T)
}



map(concat_mats, ~concat_fn(input_f, outf, .x))
map(keep_mats, ~concat_fn(input_f, outf, .x, keep=T))

map(concat_vecs, ~concat_fn(input_f, outf, .x))
map(keep_vecs, ~concat_fn(input_f, outf, .x, keep=T))

map(concat_df, ~concat_fn(input_f, outf, .x))
map(keep_df, ~concat_fn(input_f, outf, .x, keep=T))


## dosage_df <- tibble(filename=dosagef,chrom=chrom) %>%
##     arrange(chrom)

## num_files <- length(dosagef)
## stopifnot(all(file.exists(dosagef)))
## file_l <- list(filename = dosage_df$filename,
##                datapath = rep(mat_path, num_files))
## snp_df <- map_df(file_l$filename,
##                  read_df_h5,
##                  "SNPinfo",
##                  subcols=c("MAF","allele","chr","pos","snp_id"))  %>%
##     rename(orig_snp_id=snp_id)
##                                         #sample_df <- map_df(file_l$filename,read_df_h5,"SampleInfo")


## file_p <- transpose(file_l)
## concat_mats(outf,mat_path,file_p,margin = 1)
## snp_df <- mutate(snp_df,snp_id=1:n(),chr=as.integer(chr),pos=as.integer(pos))
## write_df_h5(snp_df,outf,"SNPinfo")

## all_1 <- do.call("cbind",map(file_l$filename,~read_matrix_h5(.x,mat_path,subset_rows=13L)))



## read_delim(sample_idf, delim = "\t") %>% write_df_h5(outf,"SampleInfo")

## fi <- read_matrix_h5(outf,mat_path,subset_rows=13L)

## #write_df_h5(read_df_h5(file_l$filename,"SampleInfo"),outf,"SampleInfo")
## stopifnot(all.equal(all_1,fi))
