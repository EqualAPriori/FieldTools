#!/usr/bin/env python3

import sys
sys.path.append("/home/lequieu/Work/tools/lib")

import iotools as io
import numpy as np

import pdb

if __name__ == "__main__":

    import argparse as ap
    parser = ap.ArgumentParser(description='Analyze a Trajectory')
    parser.add_argument('-i','--infilename',default='fields.dat',help='Input filename containing formatted Field data')
    parser.add_argument('-o','--outfilename',default='fields2d.dat',help='Output filename containing formatted Field data')
    parser.add_argument('--Nz',default=32,type=int,help='New number of planewaves in y dimension')
    parser.add_argument('--Lz',default=8.0,type=float, help='New boxlength in y dimension')
    args = parser.parse_args()

    #infilename = "fields.dat"
    #outfilename = "fields2d.dat"
    infilename = args.infilename   #"density.dat"
    outfilename = args.outfilename #"density2d.dat"
    Nz = args.Nz                   #64   # number of plane waves in new dimension
    Lz = args.Lz                   #5.0  # length of box in new dimension

    infiletype=infilename.split('.')[-1]
    outfiletype=outfilename.split('.')[-1]
    if (infiletype != "bin" and infiletype != 'dat'):
        raise RuntimeError(f"input file ({infilename}) must have file type of .bin or .dat")
    if (outfiletype != "bin" and outfiletype != 'dat'):
        raise RuntimeError(f"output file ({outfilename}) must have file type of .bin or .dat")

    if infiletype == 'dat':
        coords, fields = io.ReadDatFile(infilename)
    else:
        coords, fields = io.ReadBinFile(infilename)

    #cb, fb = io.ReadBinFile('density.bin')
    #cd, fd = io.ReadDatFile('density.dat')
    #pdb.set_trace()
    
    # this is where the fields are replicated
    origDim = len(coords.shape) - 1
    assert (origDim == 2), "input fields must be 2 dimensional"
    newDim = origDim + 1
    nfields = fields.shape[origDim]
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
                
    #pdb.set_trace()
    import matplotlib.pyplot as plt
    plt.pcolormesh(coordsrep[:,:,10,0], coordsrep[:,:,10,1],fieldsrep[:,:,10,0])
    plt.show()

    #pdb.set_trace()
    #io.WriteDatFile(outfilename,coords,fields,iskspace=False, iscomplex=False)
    #io.WriteDatFile("newfields.dat",coords,fields,iskspace=False, iscomplex=True)

    # just assumes not kspace or complex, not a great long term solution
    if outfiletype == 'dat':
        io.WriteDatFile(outfilename,coordsrep,fieldsrep,iskspace=False, iscomplex=False)
    else:
        io.WriteBinFile(outfilename,coordsrep,fieldsrep,iskspace=False, iscomplex=False)



