                                        #Rprof(filename=snakemake@output[["proff"]],append=F)
# save.image()
# stop()
library(SeqArray)
library(SeqSupport)
library(RSSp)

gdsf <- snakemake@input[["gdsf"]]
ymatf <- snakemake@input[["ymatf"]]
bhf <- snakemake@output[["bhf"]]



gds <- seqOpen(gdsf, readonly = T)

## ## LD_region <- seqGetData(gds, "annotation/info/LD_chunk")
## ## good_LD <- LD_region == LDchunk
## ## stopifnot(sum(good_LD) > 0)
## ## seqSetFilter(gds, variant.sel = good_LD)
## snp_id <- seqGetData(gds,var.name="variant.id")
tparam_df <- read_df_h5(ymatf, "SimulationInfo")
ymat <- read_2d_mat_h5(ymatf, "trait", "ymat")
#p <- length(LD_region)

mat_res<- gen_bhat_se_block_gds(gds,
                                ymat,
                                cores,
                                tparam_df,
                                na.rm = T)

write_mat_h5(LDchunkf, as.character(LDchunk),
             "uh",
             data = mat_res[["bias_uh_mat"]],
             deflate_level = 0,
             doTranspose = F)

write_mat_h5(LDchunkf, as.character(LDchunk),
             "se",
             data = mat_res[["se_mat"]],
             deflate_level = 0,
             doTranspose = F)
#snp_id <- as.integer((1:p))[good_LD]
write_vec(LDchunkf,as.character(LDchunk),"snp_id",data = snp_id,deflate_level = 4L)
write_df_h5(tparam_df,"SimulationInfo",LDchunkf)
#Rprof(NULL)
