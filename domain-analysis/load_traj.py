#!/usr/bin/env python3

import sys
sys.path.append("/home/lequieu/Work/tools/lib/")

import iotools as io
import pdb
import viztools as viz
import numpy as np
from domaintools import DomainAnalyzer
import glob

def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    import re
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def mapDomainIndicies(comA, comB):
    print("HELLO!")
    nA = comA.shape[0]
    nB = comB.shape[0]
    if nA != nB:
        raise RuntimeError("Number of domains is not the same between two frames ({0} != {1})".format(nA,nB))
    
    domainMap = np.zeros(nA,dtype=np.int32)
    dist = np.zeros((nA,nB))
    for i in range(nA):
        dmin = 1e30
        for j in range(nB):
           d = np.abs(comA[i] - comB[i])
           if d < dmin:
               dmin = d
               jstar=j
           
        domainMap[i] = jstar 
    return domainMap

if __name__ == "__main__":
    
    #import argparse as ap
    #parser = ap.ArgumentParser(description='Analyze a Trajectory')
    #parser.add_argument('infiles',metavar='inputfiles',nargs='*',default=['./density.bin'],help='Input filename containing unformatted Field data')
    #parser.add_argument('-v','--verbose',action='store_true',default=False,help='Print lots of extra info')
    ## Parse the command-line arguments
    #args=parser.parse_args()

    #infiles = args.infiles
    infiles = glob.glob('density*bin')
    # parser
    infiles.sort(key=natural_keys)
    print(infiles)

    # figure out maximum frame number in inputfiles (used for zero padding)
    # OPTION: could just set this to be some big value (i.e 20)
    maxdigit=0
    for infile in infiles:
        ndigit=sum( c.isdigit() for c in infile)  
        if ndigit > maxdigit: maxdigit = ndigit
    
    nblocksize = 20 # number of frames to average together before
    nframes = len(infiles)
    nblocks = int(nframes / nblocksize)
    iblock = 0
    ndomains = np.zeros(nblocks,dtype=np.int32)

    count = 0
    for iframe, infile in enumerate(infiles):
        print("Analyzing {}".format(infile))

        coords, fields = io.ReadBinFile(infile)

        if iframe == 0:
            sumcoords = np.copy(coords)
            sumfields = np.copy(fields)
        else:
            sumcoords = sumcoords + coords
            sumfields = sumfields + fields
        count +=1

        if (iframe+1) % nblocksize == 0:
            sumcoords /= count
            sumfields /= count

            viz.writeVTK("block{}.vtk".format(iblock), sumcoords, sumfields)

            domainanalyzer = DomainAnalyzer(sumcoords,sumfields)
            ndomains[iblock], com, a, v, IQ = domainanalyzer.getDomainStats(plotMesh=True)
            if iblock == 0:
                area = np.zeros((nblocks,ndomains[iblock]+5))
                vol = np.zeros((nblocks,ndomains[iblock]+5))

            if iblock == 0:
               area[iblock][0:len(a)] = a 
               vol[iblock][0:len(v)] = v
            elif iblock and (ndomains[iblock] == ndomains[iblock-1]):
               area[iblock][0:len(a)] = a 
               vol[iblock][0:len(v)] = v
            else:
                print("Warning: Number of domains not constant throughout trajectory {} {}".format(ndomains[iblock],ndomains[iblock-1]))

            sumcoords = np.zeros(sumcoords.shape)
            sumfields = np.zeros(sumfields.shape)
            count = 0
            iblock +=1
        
            np.savetxt("ndomains.dat",ndomains)
            np.savetxt("area.dat",area)
            np.savetxt("vol.dat",vol)

