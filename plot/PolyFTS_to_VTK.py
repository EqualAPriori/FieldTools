#!/usr/bin/env python3

import sys
import os
#sys.path.append("/home/lequieu/Work/tools/lib")
mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
sys.path.append(libpath)

import numpy as np
from os import path
import logging
import re
import pdb


import viztools
import iotools 

def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]


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
      head,tail = path.split(infile)
      ndigit=sum( c.isdigit() for c in tail)  
      if ndigit > maxdigit: maxdigit = ndigit

  # for each input file, generate a VTK file
  for infile in args.infiles:

      # Generate the output file name automatically (padded with zeros)
      head,tail = path.split(infile)
      outfile,ext = path.splitext(tail)
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
      if iotools.TestBinFile(infile):
        #ndim, Nx, orthorhombic, M, nfields, AllCoords, AllFields = iotools.ReadBinFile(infile)
        AllCoords, AllFields = iotools.ReadBinFile(infile)
      else:
        #ndim, Nx, orthorhombic, M, nfields, AllCoords, AllFields = iotools.ReadDatFile(infile)
        AllCoords, AllFields = iotools.ReadDatFile(infile)

      logging.info("Orthorhombic cell? {}".format(orthorhombic))
      
      print ("Outputting to Legacy VTK formatted file {}".format(outfile))
      #viztools.writeVTK(outfile, Nx, orthorhombic, M, AllCoords, AllFields)
      viztools.writeVTK(outfile, AllCoords, AllFields)

