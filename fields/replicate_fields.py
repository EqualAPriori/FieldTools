#!/usr/bin/env python3

import sys
import os
#sys.path.append("/home/lequieu/Work/tools/lib/")
mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
sys.path.append(libpath)

import fieldtools 
import iotools as io
import numpy as np

if __name__ == "__main__":

    '''Simple example showing use of fieldstools.replicatefields()
    
    '''
    
    infile="density.bin"
    #ndim, Nx, orthorhombic, M, nfields, AllCoords, AllFields = io.ReadBinFile(infile)
    coords, fields = io.ReadBinFile(infile)
    __ndim = len(coords.shape) - 1
    
    if (__ndim != 2):
        raise  NotImplementedError('implement me!')

    newcoords,newfields = fieldtools.replicate_fields(coords,fields,(2,2))


    

    # now make a little plot
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect(1)

    xx = newcoords[:,:,0]
    yy = newcoords[:,:,1]
    im=ax.pcolormesh(xx,yy,newfields[:,:,0].T)
    fig.colorbar(im,ax=ax)
    plt.savefig('fig_replicated_fields.png')
    #plt.show()
    




