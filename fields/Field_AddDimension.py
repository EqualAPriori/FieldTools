#!/usr/bin/env python3

import sys
import os
sys.path.append("/home/lequieu/Work/tools/lib")
mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
sys.path.append(libpath)

import iotools as io
import numpy as np

import pdb

if __name__ == "__main__":

    import argparse as ap
    parser = ap.ArgumentParser(description='Add dimensions to an existing field file')
    parser.add_argument('-i','--infilename',help='Input filename containing formatted Field data')
    parser.add_argument('-o','--outfilename',help='Output filename containing formatted Field data')
    parser.add_argument('-d','--dimension',required=True, type=int,help='Dimension of output fields file')

    parser.add_argument('--Ny',default=32 ,type=int,help='New number of planewaves in y dimension')
    parser.add_argument('--Ly',default=8.0,type=float, help='New boxlength in y dimension')

    parser.add_argument('--Nz',default=32 ,type=int,help='New number of planewaves in z dimension')
    parser.add_argument('--Lz',default=8.0,type=float, help='New boxlength in z dimension')
    args = parser.parse_args()

    #infilename = "fields.dat"
    #outfilename = "fields2d.dat"
    infilename = args.infilename   #"density.dat"
    outfilename = args.outfilename #"density2d.dat"
    newDim = args.dimension
    if newDim >= 2:
        Ny = args.Ny   # number of plane waves in new dimension
        Ly = args.Ly   # length of box in new dimension
    if newDim == 3:
        Nz = args.Nz   # number of plane waves in new dimension
        Lz = args.Lz   # length of box in new dimension

    infiletype=infilename.split('.')[-1]
    outfiletype=outfilename.split('.')[-1]
    if (infiletype != "bin" and infiletype != 'dat'):
        raise RuntimeError(f"input file ({infilename}) must have file type of .bin or .dat")
    if (outfiletype != "bin" and outfiletype != 'dat' and outfiletype != 'in'):
        raise RuntimeError(f"output file ({outfilename}) must have file type of .bin, .dat or .in")

    if infiletype == 'dat':
        coords, fields = io.ReadDatFile(infilename)
    else:
        coords, fields = io.ReadBinFile(infilename)

    #cb, fb = io.ReadBinFile('density.bin')
    #cd, fd = io.ReadDatFile('density.dat')
    #pdb.set_trace()
    
    # check newDim relative to origDim
    origDim = len(coords.shape) - 1
    assert (origDim < newDim), "output dimension must be greater than input dimension"
    
    if origDim == 1:
        # convert from 1d to 2d
        
        griddim = (coords.shape[0],)
        nfields = fields.shape[origDim]
        Nx = griddim[0]
        dy=Ly/Ny
        coordsrep = np.zeros((griddim[0],Ny,2))
        fieldsrep = np.zeros((griddim[0],Ny,nfields))
        for ix in range(Nx):
            for iy in range(Ny):
                coordsrep[ix,iy,0] = coords[ix][0]
                coordsrep[ix,iy,1] = iy*dy
                fieldsrep[ix,iy,:] = fields[ix,:]
    
        # if going from 1d to 3d, then need to replace fields with fieldsrep
        if newDim == 3:
           fields = np.copy(fieldsrep)
           coords = np.copy(coordsrep)

    if newDim == 3:
        nfields = fields.shape[2] # no matter what origDim was, field dim should be 2 here
        Nx = coords.shape[0]
        Ny = coords.shape[1]
        dz=Lz/Nz
        coordsrep = np.zeros((Nx,Ny,Nz,newDim))
        fieldsrep = np.zeros((Nx,Ny,Nz,nfields))

        for ix in range(Nx):
            for iy in range(Ny):
                for iz in range(Nz):
                    coordsrep[ix,iy,iz,0] = coords[ix,iy][0]
                    coordsrep[ix,iy,iz,1] = coords[ix,iy][1]
                    coordsrep[ix,iy,iz,2] = iz*dz
                    fieldsrep[ix,iy,iz,:] = fields[ix,iy,:]
                
    
    # now write ourput
    # just assumes not kspace or complex, not a great long term solution
    if outfiletype == 'dat' or outfiletype == 'in':
        io.WriteDatFile(outfilename,coordsrep,fieldsrep,iskspace=False, iscomplex=False)
    else:
        io.WriteBinFile(outfilename,coordsrep,fieldsrep,iskspace=False, iscomplex=False)



