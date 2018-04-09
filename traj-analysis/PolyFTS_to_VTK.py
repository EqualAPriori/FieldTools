#!/usr/bin/env python3

import sys
import numpy as np
from os import path
from FieldIO import *
import logging
import re

def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

 
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

def writeVTK(outfile, Nx, orthorhombic, M, AllCoords, AllFields):
    # Write Legacy VTK file.
    o = open(outfile,"w")

    # get mesh spacing
    if orthorhombic == True:
        spacing = [None]*ndim
        if ndim == 1:
            spacing[0] = AllCoords[0][1]
        if ndim == 2:
            spacing[0] = AllCoords[0][Nx[1]]
            spacing[1] = AllCoords[1][1]
        if ndim == 3:
            spacing[0] = AllCoords[0][Nx[1]*Nx[2]]
            spacing[1] = AllCoords[1][Nx[2]]
            spacing[2] = AllCoords[2][1]
        logging.info("Mesh spacing       {}".format(spacing))
  
  
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
            o.close();o = open(outfile,"ab")
            np.savetxt(o, AllFields[i], fmt="%.11f")
            o.close();o = open(outfile,"a")
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
  

if __name__ == "__main__":
  # For command-line runs, build the relevant parser
  import argparse as ap
  parser = ap.ArgumentParser(description='Convert PolyFTS Field file to Legacy VTK format')
  parser.add_argument('infiles',metavar='inputfiles',nargs='*',default=['./density.bin'],help='Input filename containing unformatted Field data')
  parser.add_argument('-v','--verbose',action='store_true',default=False,help='Print lots of extra info')
  # Parse the command-line arguments
  args=parser.parse_args(sys.argv[1:])

  # parser
  args.infiles.sort(key=natural_keys)
  print(args.infiles)

  # figure out maximum frame number in inputfiles (used for zero padding)
  # OPTION: could just set this to be some big value (i.e 20)
  maxdigit=0
  for infile in args.infiles:
      ndigit=sum( c.isdigit() for c in infile)  
      if ndigit > maxdigit: maxdigit = ndigit

  # for each input file, generate a VTK file
  for infile in args.infiles:

      # Generate the output file name automatically (padded with zeros)
      outfile,ext = path.splitext(infile)
      outfile_alphas=''.join([c for c in outfile if c.isalpha()])
      outfile_digits=''.join([c for c in outfile if c.isdigit()])
      zeropadded=str("{val:0>{width}}".format(val=outfile_digits,width=maxdigit))
      outfile = outfile_alphas + zeropadded  + ".vtk"
      
      # check if vtk dont already exist?
      # IMPLEMENT ME
      
      if args.verbose:
          logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
      else:
          logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.WARNING)

      orthorhombic = True

      # Check whether the input file exits, and whether it is a binary PolyFTS file or a formatted on.
      # Dispatch reading to relevant function.
      # Open as binary first
      print("Reading input file {}".format(infile))
      if TestBinFile(infile):
        ndim, Nx, orthorhombic, M, nfields, AllCoords, AllFields = ReadBinFile(infile)
      else:
        ndim, Nx, orthorhombic, M, nfields, AllCoords, AllFields = ReadDatFile(infile)

      logging.info("Orthorhombic cell? {}".format(orthorhombic))
      
       # Generate the mapping from PolyFTS (z fastest in 3D) to VTK (x fastest)
      IdxMap = get_PolyFTS_to_VTK_IdxMap(M,Nx)
      # Remap field samples from PolyFTS order to VTK order
      AllFields = AllFields[:,IdxMap]
      print ("Outputting to Legacy VTK formatted file {}".format(outfile))
      writeVTK(outfile, Nx, orthorhombic, M, AllCoords, AllFields)

