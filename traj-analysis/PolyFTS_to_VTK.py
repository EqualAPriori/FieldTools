#!/usr/bin/env python3

import sys
import numpy as np
from os import path
from FieldIO import *
import logging

def get_PolyFTS_to_VTK_IdxMap(M,Nx):
    idx=np.zeros(M,dtype=np.uint64)
    if ndim == 1:
        for ix in range(Nx[0]):
            idx[ix] = ix
    elif ndim == 2:
        #looks good
        m=0
        for iy in range(Nx[1]):
          for ix in range(Nx[0]):
            idx[m] = ix*Nx[1] + iy
            m+=1
    elif ndim == 3:
        m=0
        for iz in range(Nx[2]):
          for iy in range(Nx[1]):
            for ix in range(Nx[0]):
                idx[m] = ix*Nx[1]*Nx[2] + iy*Nx[2] + iz
                #idx[m] = iz*Nx[0]*Nx[1] + iy*Nx[0] + ix
                m+=1
    return idx

if __name__ == "__main__":
  # For command-line runs, build the relevant parser
  import argparse as ap
  parser = ap.ArgumentParser(description='Convert PolyFTS Field file to Legacy VTK format')
  parser.add_argument('infile',metavar='inputfile',nargs='?',default='./density.bin',help='Input filename containing unformatted Field data')
  parser.add_argument('-v','--verbose',metavar='inputfile',nargs='?',default='./density.bin',help='Input filename containing unformatted Field data')
  # Parse the command-line arguments
  args=parser.parse_args(sys.argv[1:])
  # Generate the output file name automatically
  outfile,ext = path.splitext(args.infile)
  outfile = outfile + ".vtk"
  
  logging.basicConfig(level=logging.INFO)
  orthorhombic = True

  # Check whether the input file exits, and whether it is a binary PolyFTS file or a formatted on.
  # Dispatch reading to relevant function.
  # Open as binary first
  print("Reading input file {}".format(args.infile))
  if TestBinFile(args.infile):
    ndim, Nx, h, M, nfields, AllFields = ReadBinFile(args.infile)
    # Check whether the cell is orthorhombic
    for i in range(ndim):
        for j in range(ndim):
            if i == j:
                continue
            if abs(h[i][j]) > 1e-8:
                orthorhombic = False
    logging.info("Orthorhombic cell? {}".format(orthorhombic))
    if orthorhombic == True:
        # Orthorhombic cells will use the VTK structured points format - need grid spacing
        spacing = [h[i][i]/Nx[i] for i in range(ndim)]
        logging.info("Mesh spacing       {}".format(spacing))
    else:
        # non orthorhombic cells with use the VTK structured grid format - need the x,y,z position of all points
        AllCoords = calcAllCoords(ndim,Nx,h,M)
  else:
    ndim, Nx, M, nfields, AllCoords, AllFields = ReadDatFile(args.infile)
    # This routine is the fallback reader and should quit if not a valid dat file
    #
    # Check whether the cell is orthorhombic
    # In future the .dat file format should be updated to have the cell tensor.
    spacing = [None]*ndim
    if ndim == 1:
        spacing[0] = AllCoords[0][1]
    if ndim == 2:
        if AllCoords[0][1] != 0.:
            orthorhombic = False
        if AllCoords[1][Nx[1]] != 0.:
            orthorhombic = False
        spacing[0] = AllCoords[0][Nx[1]]
        spacing[1] = AllCoords[1][1]
    if ndim == 3:
        # Check that the first non-zero mesh point has only a z translation
        if AllCoords[0][1] != 0. or AllCoords[1][1] != 0.:
            orthorhombic = False
        # Check that the first y translation has no x,z
        if AllCoords[0][Nx[2]] != 0 or AllCoords[2][Nx[2]] != 0:
            orthorhombic = False
        # Check that the first x translation has no y,z
        if AllCoords[1][Nx[2]*Nx[1]] != 0 or AllCoords[2][Nx[2]*Nx[1]] != 0:
            orthorhombic = False
        spacing[0] = AllCoords[0][Nx[1]*Nx[2]]
        spacing[1] = AllCoords[1][Nx[2]]
        spacing[2] = AllCoords[2][1]
    logging.info("Orthorhombic cell? ".format(orthorhombic))
    if orthorhombic:
        logging.info("Mesh spacing       ".format(spacing))

  # Generate the mapping from PolyFTS (z fastest in 3D) to VTK (x fastest)
  IdxMap = get_PolyFTS_to_VTK_IdxMap(M,Nx)
  # Remap field samples from PolyFTS order to VTK order
  AllFields = AllFields[:,IdxMap]

  # for debugging
  #AllFields[0] = np.array([float(i)/M for i in range(M)])
  #np.savetxt('idx.dat',IdxMap)
  #np.savetxt('fields.dat',AllFields.transpose())
  #pdb.set_trace()

  # Legacy VTK file.
  print ("Outputting to Legacy VTK formatted file {}".format(outfile))
  o = open(outfile,"w")

  if orthorhombic:
    o.write("# vtk DataFile Version 3.0\n")
    o.write("PolyFTS field data\n")
    o.write("ASCII\n")
    o.write("DATASET STRUCTURED_POINTS\n")
    if ndim == 1:
        o.write("DIMENSIONS {} 1 1\n".format(*Nx))
        o.write("ORIGIN 0\n")
        o.write("SPACING {} 0 0\n".format(*spacing))
    elif ndim == 2:
        o.write("DIMENSIONS {} {} 1\n".format(*Nx))
        o.write("ORIGIN 0 0 0\n")
        o.write("SPACING {} {} 0\n".format(*spacing))
    elif ndim == 3:
        o.write("DIMENSIONS {} {} {}\n".format(*Nx))
        o.write("ORIGIN 0 0 0\n")
        o.write("SPACING {} {} {}\n".format(*spacing))
    o.write("POINT_DATA {0}\n".format(M))
    for i in range(nfields):
        o.write("SCALARS field{0} float 1\n".format(i))
        o.write("LOOKUP_TABLE default\n")
        np.savetxt(o, AllFields[i], fmt="%.11f")
  else:
    o.write("# vtk DataFile Version 3.0\n")
    o.write("PolyFTS field data\n")
    o.write("ASCII\n")
    o.write("DATASET STRUCTURED_GRID\n")
    if ndim == 1:
        o.write("DIMENSIONS {} 1 1\n".format(*Nx))
    elif ndim == 2:
        o.write("DIMENSIONS {} {} 1\n".format(*Nx))
    elif ndim == 3:
        o.write("DIMENSIONS {} {} {}\n".format(*Nx))
    o.write("POINTS {} float\n".format(M))
    # Remap the mesh coordinates to VTK order
    AllCoords = AllCoords[:,IdxMap]
    #np.savetxt('coords.dat',AllCoords.transpose()) # Debug
    o.close(); o = open(outfile,'ab')
    np.savetxt(o, AllCoords.transpose(), fmt="%0.11f")
    o.close(); o = open(outfile,'a')

    o.write("\nPOINT_DATA {0}\n".format(M))
    for i in range(nfields):
    #for i in range(1):
        # write as scalar
        #o.write("SCALARS field{0} float 1\n".format(i)) #STARTED TO SEG FAULT WHEN TRYING TO READ SCALARS
        #o.write("LOOKUP_TABLE default\n")
        #np.savetxt(o, AllFields[i], fmt="%14.11f")
        #np.savetxt(o, np.vstack([AllCoords,AllFields[i]]).transpose(), fmt="%14.11f")

        # write as vector
        o.write("VECTORS field{0} float\n".format(i)) #writing as vector fixed the seg fault
        tmp=np.zeros((3,M))
        tmp[0,:] = AllFields[i]
        o.close(); o = open(outfile,'ab')
        np.savetxt(o, tmp.transpose(), fmt="%14.11f")
        o.close(); o = open(outfile,'a')
  o.close()

