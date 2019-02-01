library(SeqArray)
input_f_a  <- "/scratch/t.cri.nknoblauch/polyg_scratch/gds/FRAM/FRAM.chr19.gds"
input_f_b  <- "/scratch/t.cri.nknoblauch/polyg_scratch/gds/EUR/EUR.chr19.gds"


input_f_a <- snakemake@input[["gds_a"]]
input_f_b <- snakemake@input[["gds_b"]]

gds_a  <- seqOpen(input_f_a)
gds_b  <- seqOpen(input_f_b)
