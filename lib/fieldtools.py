#!/usr/bin/env python3
'''
    Joshua Lequieu <lequieu@mrl.ucsb.edu>

    Functions to do usefule things to fields
'''

import numpy as np
import pdb


def replicate_fields(coords,fields,nreplicates):
    '''take in coords and fields and nreplicates = [nrepx,nrepy,nrepz]
       
       return new coords and fields that have been replicated
    '''
   
   #note: double underline __ are variable definitions that I copy-pasted from domaintools
   # I should really write a function that extracts all of this info from a coords,fields pair
   #    ^ but I dont know where to put it...and I dont want everything depending on it...
    __ndim = len(coords.shape) - 1
    __Nx = coords.shape[:__ndim]

    if __ndim == 1:
        raise NotImplementedError("ndim == 1 not implemented")
    elif __ndim == 2:
        __gridspacing = (coords[1,0][0], coords[0,1][1])
        __hvoxel = np.array([coords[1,0],coords[0,1]])
    elif __ndim == 3:
        __gridspacing = (coords[1,0,0][0], coords[0,1,0][1], coords[0,0,1][2])
        __hvoxel = np.array([coords[1,0,0],coords[0,1,0],coords[0,0,1]])
    __hcell = __hvoxel * __Nx

    
    assert(len(nreplicates) == __ndim)
    
    # replicating the fields is simple using numpy
    # don't want to replicate the actual values (last dimenstion) so append [1]
    reps = list(nreplicates) + [1]
    newfields= np.tile(fields,reps)

    # replicating the coords is a bit more complicated
    # first replicate in 'box' coordinates (which are just the i,j index of the voxel)
    # then get the xy position of each by using hvoxel
    newNx = np.array(__Nx) * nreplicates
    if __ndim == 2:
        
        x = np.arange(newNx[0])
        y = np.arange(newNx[1])
        yy,xx = np.meshgrid(x,y)
        # nice one-liner to rotate all of xx and yy using hvoxel
        xxrot,yyrot = np.einsum('ji, mni -> jmn', __hvoxel.T, np.dstack([xx, yy]))

        newcoords = np.dstack([xxrot,yyrot])

    elif __ndim == 3:
        raise NotImplementedError("three dimenstions not implemented yet!")
  

    return newcoords,newfields

    
