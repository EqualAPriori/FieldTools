#!/usr/bin/env python3
import sys
sys.path.append("/home/lequieu/Work/tools/lib/")

import iotools as io
import pdb
import viztools as viz
import numpy as np
from domaintools import DomainAnalyzer

if __name__ == "__main__":
    
    infile="density.bin"
    #ndim, Nx, orthorhombic, M, nfields, AllCoords, AllFields = io.ReadBinFile(infile)
    AllCoords, AllFields = io.ReadBinFile(infile)

    np.save("coords.npy", AllCoords)
    np.save("densityfields.npy", AllFields)

    infile="fields.bin"
    #ndim, Nx, orthorhombic, M, nfields, AllCoords, AllFields = io.ReadBinFile(infile)
    AllCoords, AllFields = io.ReadBinFile(infile)

    #np.save("coords.npy", AllCoords)
    np.save("wfields.npy", AllFields)
