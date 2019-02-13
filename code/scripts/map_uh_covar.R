# save.image("sim.RData")
# stop()

library(ldshrink)
library(tidyverse)
library(SeqSupport)
library(EigenH5)


h5f <- snakemake@input[["h5f"]]
covarf <- snakemake@input[["covarf"]]
uhf <- snakemake@output[["uhf"]]
subsnpf <- snakemake@input[["subsnpf"]]
subgwasf <- snakemake@input[["subgwasf"]]
ncovar <- as.integer(snakemake@params[["ncovar"]] %||% 0)
y_grp  <- snakemake@params[["y_grp"]] %||% "TraitInfo"

snp_chunksize <- snakemake@params[["chunksize_snp"]]
if(is.null(snp_chunksize)){
    snp_chunksize <- 20000
}


stopifnot(!is.null(h5f),!is.null(uhf))
ymatf <- normalizePath(snakemake@input[["ymatf"]])
stopifnot(!is.null(ymatf))
cores <- snakemake@threads

cat("Using ",cores," cores\n")

map_eqtl_h5(snp_h5=h5f,
            exp_h5=ymatf,
            covar_h5=covarf,exp_info=y_grp,
            uh_h5=uhf,ncovar=ncovar,exp_path="ymat/trait",subsnpf=subsnpf,subgwasf=subgwasf,snp_chunksize=as.integer(snp_chunksize),threads=cores)
cat("Done!")


## if(is.null(y_grp)){
##     y_grp="SimulationInfo"
## }






## tparam_df <- EigenH5::read_df_h5(ymatf, y_grp)

## p <- nrow(snp_df)



## ## N <- as.integer(unique(tparam_df$n))

## g <- as.integer(nrow(tparam_df))


## dim_dos <- dim_h5(h5f,"dosage")
## SNPfirst <- dim_dos[1]==p



## ## stopifnot(length(ind_v)==N)





## chunksize <- as.integer(5000)
## num_chunks <- ceiling(p/chunksize)
## cat("Chunking SNP data\n")
## snp_df <- mutate(snp_df,snp_chunk=gl(n = num_chunks,k=chunksize,length = p),snp_chunk_id=1:n())
## snp_l <-  split(select(snp_df,snp_id,snp_chunk_id),snp_df$snp_chunk)

## snp_lff  <- snp_l %>% map(~list(subset_rows=.x$snp_id,
##                                 filename=h5f,
##                                 subset_cols=ind_v,
##                                 datapath="dosage"))


## exp_lff <- list(list(filename=ymatf,datapath="trait/ymat"))

## uh_lff <-  snp_l %>% map(~list(subset_rows=.x$snp_chunk_id,filename=uhf,datapath="uh"))
## se_lff <-  map(uh_lff,~update_list(.,datapath="se"))









## cat("Creating output matrices\n")
## EigenH5::create_matrix_h5(uhf,"uh",numeric(),dims=c(p,g),chunksizes=c(1000L,g))
## EigenH5::create_matrix_h5(uhf,"se",numeric(),dims=c(p,g),chunksizes=c(1000L,g))

## stopifnot(nrow(tparam_df)>0)
## #stopifnot(!is.null(snp_df$AF),
## #          all(!is.na(snp_df$AF)))

## cat("Writing simulation/data info\n")
## write_df_h5(snp_df,uhf,"SNPinfo")
## write_df_h5(tparam_df,uhf,y_grp)
## pl <- snakemake@wildcards
## pl <- as_tibble(pl[names(pl)!=""])
## write_df_h5(pl,uhf, "Wildcards")

## cat("Mapping traits\n")
## cat("SNPfirst=",SNPfirst,"\n")
## cat("EXPfirst=",EXPfirst,"\n")

## SeqSupport::map_eQTL_chunk_h5(snp_lff,exp_lff,uh_lff,se_lff,covarmat,
##                               list(SNPfirst=SNPfirst,
##                                    EXPfirst=EXPfirst))


## cat("Checking uh\n")
## tuh <- EigenH5::read_matrix_h5(uhf,"uh",
##                                subset_rows=sort(sample(1:p,min(p,100),replace=F)),
##                                subset_cols=sort(sample(1:g,min(g,100),replace=F)))
## stopifnot(all(!is.na(c(tuh))))
