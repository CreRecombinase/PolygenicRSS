
library(SeqSupport)
library(EigenH5)
library(tidyverse)

input_f <- snakemake@input[["input_f"]]
outputf  <- snakemake@output[["h5f"]]


pl <- snakemake@wildcards
pl <- as_tibble(pl[names(pl)!=""])
EigenH5::write_df_h5(pl,filename=outputf, "Wildcards")

input_f_a <- input_f[1]

tp_df <- read_df_h5(input_f_a,"SimulationInfo") %>% mutate(trait_id=1:n())
write_df_h5(tp_df,outputf,"SimulationInfo")
ymat <- read_matrix_h5(input_f_a,"genetic_trait/genetic_ymat")
if(length(input_f)>1){

    rest_f <- input_f[-1]

    for( tif in rest_f){
        ymat <- ymat+read_matrix_h5(tif,"genetic_trait/genetic_ymat")
    }

}


real_ymat <- gen_sim_resid(t(ymat),tp_df)

write_matrix_h5(real_ymat,outputf,"ymat/trait")
