#!/usr/bin/env python3
import sys
sys.path.append("/home/lequieu/Work/tools/lib/")

import iotools as io
import pdb
import viztools as viz
import numpy as np
from skimage import measure

from domaintools import DomainAnalyzer

if __name__ == "__main__":
    
    infile="density.bin"
    #ndim, Nx, orthorhombic, M, nfields, AllCoords, AllFields = io.ReadBinFile(infile)
    coords, fields = io.ReadBinFile(infile)

    domainanalyzer = DomainAnalyzer(coords,fields)
    domainanalyzer.meshAllDomains(datafile='contours.dat',plotfile="mesh_all.png")

