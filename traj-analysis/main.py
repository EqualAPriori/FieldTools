#!/usr/bin/env python3

import FieldIO
import Configuration
import pdb
import PolyFTS_to_VTK
import numpy as np
from domaintools import DomainAnalyzer

if __name__ == "__main__":
    
    infile="density.bin"
    # TODO flip AllCoords and AllFields so that each grid point is in first dimension
    ndim, Nx, orthorhombic, M, nfields, AllCoords, AllFields = FieldIO.ReadBinFile(infile)

    # reshape AllCoords and Allfields 
    AllCoords = np.reshape(AllCoords.T ,list(Nx) + [ndim])
    AllFields = np.reshape(AllFields.T ,list(Nx) + [nfields])
    domainanalyzer = DomainAnalyzer(AllCoords,AllFields)
    stats = domainanalyzer.getDomainStats()
    
    np.savetxt("stats.dat",stats)

