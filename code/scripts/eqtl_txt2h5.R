library(tidyverse)
library(EigenH5)
library(progress)
input_file <- "/scratch/t.cri.nknoblauch/polyg_scratch/eqtl-gene-annot.txt"
output_file <- "/scratch/t.cri.nknoblauch/polyg_scratch/eqtl-gene-annot.h5"
eqtl_df <- read_csv(input_file,progress=T,n_max=189269666)

cn <- colnames(eqtl_df)
bad_cols <- character()
ext_ls <- character()


pb <- progress_bar$new(total=length(cn))
for(i in cn){
    if(!i %in% ext_ls){
        if(typeof(eqtl_df[[i]])=="character"){
            max_nc <- max(nchar(eqtl_df[[i]]),na.rm=T)
        }else{
            max_nc  <- 0
        }
        if(max_nc<254){
            write_vector_h5(eqtl_df[[i]],output_file,paste0("eQTLinfo/",i))
            ext_ls  <- ls_h5(output_file,"eQTLinfo")
        }else{
            bad_cols <- c(bad_cols,i)

        }
    }
    pb$tick()
}
cat(paste0(c("Bad Cols:",bad_cols,collapse="\n")))
saveRDS(bad_cols,"~/Downloads/bad_cols.RDS")
