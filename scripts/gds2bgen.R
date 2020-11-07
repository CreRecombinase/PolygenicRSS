library(SeqArray)
library(gds2bgen)
library(BiocParallel)
input_f <- snakemake@input[["bgen"]]
output_f <- snakemake@output[["gds"]]
worker <- snakemake@threads

#worker <- BatchtoolsParam(workers=worker_n,cluster="torque", template="/home/t.cri.nknoblauch/torque_batchtools.tmpl",resources=list(walltime = 36000L,memory=8000))
seqBGEN2GDS(input_f,output_f,storage.option="LZ4_RA",float.type="packed8",geno=TRUE,dosage=FALSE,prob=FALSE,parallel=worker)
