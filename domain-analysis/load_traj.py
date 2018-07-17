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
    
    import argparse as ap
    parser = ap.ArgumentParser(description='Analyze a Trajectory')
    parser.add_argument('--densityfiles',metavar='inputfiles',nargs='+',default=glob.glob('density*bin'),help='Input filename containing unformatted Field data')
    parser.add_argument('-n','--navg',default=100,type=int,help='Number of frames to average the fields over before computing domains')
    parser.add_argument('-s','--save',action='store_true',help='Input filename containing unformatted Field data')
    parser.add_argument('-l','--load',action='store_true',help='Input filename containing unformatted Field data')
    parser.add_argument('--novtk',action='store_true',help='flag to turn off saving of VTK files of averaged densities')
    # Parse the command-line arguments
    args=parser.parse_args()

    #infiles = args.infiles
    saveflag = args.save
    loadflag = args.load
    novtkflag = args.novtk
    if (saveflag == True) and (loadflag == True):
        raise ValueError ("Flags --save and --load cannot both be specified")
    elif (saveflag == False) and (loadflag == False):
        # save flag defaults to true if neither loading or saving is set
        saveflag = True
    
    nblocksize = args.navg #100 # number of frames to average together

    infiles = args.densityfiles #glob.glob('density*bin')
    # parser
    infiles.sort(key=natural_keys)
    print(infiles)

    # figure out maximum frame number in inputfiles (used for zero padding)
    # OPTION: could just set this to be some big value (i.e 20)
    maxdigit=0
    for infile in infiles:
        ndigit=sum( c.isdigit() for c in infile)  
        if ndigit > maxdigit: maxdigit = ndigit
    
    nframes = len(infiles)
    nblocks = int(nframes / nblocksize)
    iblock = 0
    ndomains = np.zeros(nblocks,dtype=np.int32)

    count = 0
    for iframe, infile in enumerate(infiles):
        
        if saveflag:
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
            if saveflag:
                sumcoords /= count
                sumfields /= count
                
                if not novtkflag:
                    viz.writeVTK("block{}.vtk".format(iblock), sumcoords, sumfields)

                domainanalyzer = DomainAnalyzer(sumcoords,sumfields)
                ndomains[iblock], com, a, v, IQ = domainanalyzer.getDomainStats(plotMesh=False,useMesh=False)
                #ndomains[iblock], com, a, v, IQ = domainanalyzer.getDomainStats(plotMesh=True)

                if (ndomains[iblock] != ndomains[iblock-1]):
                   print("Warning: Number of domains not constant throughout trajectory {} {}".format(ndomains[iblock],ndomains[iblock-1]))

                stats = np.vstack((com.T,a,v,IQ)).T
                np.savetxt("stats{}.dat".format(iblock),stats, header="comx comy comz area vol IQ")
            
                #np.savetxt("com{}.dat".format(iblock),com)

                # ---------------------
                # Eventually I'd like to remove this code but I'll leave it for now for legacy reasons
                if iblock == 0:
                    area = np.zeros((nblocks,ndomains[iblock]+20))
                    vol = np.zeros((nblocks,ndomains[iblock]+20))

                area[iblock][0:len(a)] = a 
                vol[iblock][0:len(v)] = v

                np.savetxt("ndomains.dat",ndomains)
                np.savetxt("area.dat",area)
                np.savetxt("vol.dat",vol)

                # ---------------------

                # zero quantities
                count = 0
                sumcoords = np.zeros(sumcoords.shape)
                sumfields = np.zeros(sumfields.shape)

            elif loadflag:
                print("Loading stats{}.dat".format(iblock))
                stats=np.loadtxt("stats{}.dat".format(iblock))
                com = stats[:,0:3]
                a = stats[:,3]
                v  = stats[:,4]
                IQ = stats[:,5]


            
            # here com, area, vol, IQ were either computed, or loaded
            #pdb.set_trace() #check that the variables are reasonable


            # update block 
            iblock +=1
        

