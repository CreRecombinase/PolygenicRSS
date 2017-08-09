#!/usr/bin/env python

# Command line arguments can be read with sys.argv. The first element
# is skipped because it is the name of the file.

# Usage:
#
# ./script.py ex1 ex2 ex3
#

import sys
import h5py
import numpy as np
import pandas as pd
import feather

args = sys.argv[1:]
matfile = args[0]
outfeather = args[1]
if len(args)==3:
    groupname = args[2]
else:
    groupname='snp_info'
    
f= h5py.File(matfile,'r')
si=f[groupname]

tdf=[f[si[i,0]][:] for i in range(3)]
tdfd={'rsid':[int(z) for z in tdf[0][0]],'ref':[chr(x) for x in tdf[1][0]],'alt':[chr(y) for y in tdf[2][0]]}

my_df=pd.DataFrame(tdfd)
my_df.to_feather(outfeather)



