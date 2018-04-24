#!/usr/bin/env python3

import iotools as io
import pdb
import viztools as viz
import numpy as np
from domaintools import DomainAnalyzer

if __name__ == "__main__":
    
    infile="density.bin"
    #ndim, Nx, orthorhombic, M, nfields, AllCoords, AllFields = io.ReadBinFile(infile)
    AllCoords, AllFields = io.ReadBinFile(infile)

    domainanalyzer = DomainAnalyzer(AllCoords,AllFields)
    stats = domainanalyzer.getDomainStats(plotMesh=True)
    #stats = domainanalyzer.getDomainStats(plotMesh=False)
    
    np.savetxt("stats.dat",stats)

