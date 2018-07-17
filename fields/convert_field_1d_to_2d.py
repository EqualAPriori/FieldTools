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
    parser.add_argument('--Ny',default=64,type=int,help='New number of planewaves in y dimension')
    parser.add_argument('--Ly',default=10.0,type=float, help='New boxlength in y dimension')
    args = parser.parse_args()

    #infilename = "fields.dat"
    #outfilename = "fields2d.dat"
    infilename = args.infilename   #"density.dat"
    outfilename = args.outfilename #"density2d.dat"
    Ny = args.Ny                   #64   # number of plane waves in new dimension
    Ly = args.Ly                   #5.0  # length of box in new dimension

    coords, fields = io.ReadDatFile(infilename)
    
    # this is where the fields are replicated
    Dim = len(coords.shape) - 1
    assert (Dim == 1), "input fields must be 1 dimensional"
    griddim = (coords.shape[0],)
    nfields = fields.shape[Dim]
    Nx = griddim[0]
    dy=Ly/Ny
    coordsrep = np.zeros((griddim[0],Ny,2))
    fieldsrep = np.zeros((griddim[0],Ny,nfields))
    for ix in range(Nx):
        for iy in range(Ny):
            coordsrep[ix,iy,0] = coords[ix][0]
            coordsrep[ix,iy,1] = iy*dy
            fieldsrep[ix,iy,:] = fields[ix,:]

    #pdb.set_trace()
    #io.WriteDatFile(outfilename,coords,fields,iskspace=False, iscomplex=False)
    #io.WriteDatFile("newfields.dat",coords,fields,iskspace=False, iscomplex=True)

    # just assumes not kspace or complex, not a great long term solution
    io.WriteDatFile(outfilename,coordsrep,fieldsrep,iskspace=False, iscomplex=False)



