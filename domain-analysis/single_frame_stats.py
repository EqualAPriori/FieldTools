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

    domainanalyzer = DomainAnalyzer(AllCoords,AllFields)

    ndomains, com,area,vol,IQ = domainanalyzer.getDomainStats(useMesh=False,plotMesh=True)
    stats = np.vstack((com.T,area,vol,IQ)).T
    np.savetxt("stats_nomesh.dat",stats, header="comx comy comz area vol IQ")

    ndomains, com,area,vol,IQ = domainanalyzer.getDomainStats(plotMesh=True)
    stats = np.vstack((com.T,area,vol,IQ)).T
    np.savetxt("stats_mesh.dat",stats, header="comx comy comz area vol IQ")


