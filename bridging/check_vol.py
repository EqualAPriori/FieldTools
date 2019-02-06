#!/usr/bin/env python3

import sys
sys.path.append("/home/lequieu/Work/tools/lib/")

import numpy as np
import pdb
import matplotlib.pyplot as plt

def parse_file(ifilename):
    ''' Check BridgingOperator.dat to make sure that its
    '''
    with open(ifilename,'r') as f:
        for line in f:
            l=line.split()
            if 'Format version' in line:
                formatversion=float(l[-1])
                if formatversion != 3:
                    print("Error! Format version of %s must be 3" % ifilename)
                    exit
            elif 'NDim' in line:
                ndim=int(l[-1])
            elif 'nfields' in line:
                nfields=float(l[-1])
            elif 'k-space' in line and 'complex' in line:
                if "k-space data = 0 , complex data = 0" not in line:
                    print("Error! %s must have k-space data = 0, complex data = 0" % ifilename)
                    exit
    data = np.loadtxt(ifilename)

    return data,ndim,nfields         


if __name__ == "__main__":
    
    infile="BridgingOperator.dat"

    d,ndim, nfields = parse_file(infile)
    npw = d.shape[0]
    icol=ndim+2-1
    ndomain = int(np.max(d[:,icol]))
    vol=np.zeros(ndomain)
    for i in range(ndomain):
        vol[i] = np.count_nonzero(d[:,icol] == i+1)
    m=np.average(vol)
    std=np.std(vol)
    vol=np.sort(vol) 

    print(f"Mean: {m} +/- {std:.2f}, ( std/mean: {std/m:.2f}, std/npw: {std/npw:0.4f} )")
    print(f"Volumes: {vol}")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(vol,marker='o',ls='')
    plt.show()


    

