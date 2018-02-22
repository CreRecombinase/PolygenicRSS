import allel
from allel import chunked
import h5py
vcf_f = snakemake.input[0]
out_hf = snakemake.output[0]
allel.vcf_to_hdf5(vcf_f,out_hf,compression='lzf',fields="*",compression_opts=None,overwrite=True,numbers={'variants/ALT': 1},log=sys.stdout,vlen = False)
hf=h5py.File(out_hf,'r+')
gg=hf['/']
td =hf['calldata']['GT']
hg=allel.GenotypeChunkedArray(td).to_haplotypes()
chunked.copy(hg, storage='hdf5', group=gg, name='dosage',compression="lzf")
hf.close()
