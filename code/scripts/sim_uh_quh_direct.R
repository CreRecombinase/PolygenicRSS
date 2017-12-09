# save.image()
# stop()
library(SeqSupport)
library(SeqArray)
library(dplyr)

evdf <- snakemake@input[["evdf"]]
gdsf <- snakemake@input[["gdsf"]]

quh_outf <- snakemake@output[["quh_hf"]]
uh_outf <- snakemake@output[["uh_hf"]]

LDchunk <- snakemake@params[["LDchunk"]]


mfgeneid <- as.character(snakemake@params[["fgeneid"]])
pve <- as.numeric(snakemake@params[["pve"]])
bias <- as.numeric(snakemake@params[["bias"]])
nreps <- as.integer(as.numeric(snakemake@params[["nreps"]]))
n <- as.numeric(snakemake@params[["N"]])

gds <- seqOpen(gdsf, readonly = T)

p <- length(seqGetData(gds, "variant.id"))
seqClose(gds)

tparam_df <- gen_tparamdf_norm(pve,bias,nreps,n = n,p = p,mfgeneid) %>% mutate(n=n,p=p)

stopifnot(!is.null(LDchunk))

ld_snpi <- RcppEigenH5::read_ivec(evdf,"LDinfo","snp_id")
Dvec <- RcppEigenH5::read_dvec(evdf,"EVD","D")
Q <- RcppEigenH5::read_2d_mat_h5(h5file = evdf,groupname = "EVD",dataname = "Q")
Qp <- ncol(Q)
csD <-cumsum(Dvec)/Qp
stopifnot(all.equal(csD[Qp],1))

uh_quh <- sim_uh_quh_dir(tparam_df$tsigu,
                         tparam_df$tbias,
                         tparam_df$fgeneid,
                         Q,
                         Dvec,NULL,ld_snpi)
write_mat_h5(quh_outf,
             LDchunk,
             "quh",
             uh_quh$quh,
             deflate_level=0,
             doTranspose=F)
se_chunk <- matrix(1.0,length(D),length(tparam_df$fgeneid))

write_mat_h5(uh_outf,
             LDchunk,
             "se",
             se_chunk,
             deflate_level=0,
             doTranspose=F)

write_mat_h5(quh_outf,
             LDchunk,
             "se",
             se_chunk,
             deflate_level=0,
             doTranspose=F)



write_vec(h5filename = quh_outf,groupname = LDchunk,dataname = "D",data = Dvec)

write_vec(h5filename = quh_outf,groupname = LDchunk,dataname = "snp_id",data = ld_snpi)
write_df_h5(tparam_df,"SimulationInfo",quh_outf)

write_mat_h5(uh_outf, LDchunk,
             "uh",
             data = uh_quh$uh,
             deflate_level = 0,
             doTranspose = F)
#snp_id <- as.integer((1:p))[good_LD]
write_vec(uh_outf,LDchunk,"snp_id",data = ld_snpi,deflate_level = 4L)
write_df_h5(tparam_df,"SimulationInfo",uh_outf)

