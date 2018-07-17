#!/usr/bin/env python3

import struct
import numpy as np
import logging
import sys
import pdb

def TestBinFile(infilename):
    ValidBinFile = True
    i = open(infilename,"rb")
    header = struct.unpack_from("@8s",i.read(8))[0] # Check first 8 bytes
    version = struct.unpack_from("@I",i.read(5),offset=1)[0] # Note the need to skip one byte due to the null termination character that got saved in the C string
    if header != b"FieldBin" or version != 51:
        logging.critical("  {} is not an unformatted field file".format(infilename))
        ValidBinFile = False
    else:
        logging.info("  {} is an unformatted field file with version #{}".format(infilename,version))
    i.close()
    return ValidBinFile

def ReadBinFile(infilename):
    i = open(infilename,"rb")
    contents = i.read()
    i.close()

    # Check header and file version number
    header = struct.unpack_from("@8s",contents)
    if header[0] != b"FieldBin":
      sys.stderr.write("\nError: Not an unformatted Field file\n")
      sys.exit(1)

    pos = 9
    version = struct.unpack_from("@I",contents,offset=pos)
    pos = pos + 4
    if version[0] != 51:
      sys.stderr.write("\nError: Only version 51 is currently supported\n")
      sys.exit(1)

    nfields = struct.unpack_from("@I",contents,offset=pos)[0]
    pos = pos + 4
    logging.info("  # fields in file = {}".format(nfields))

    ndim = struct.unpack_from("@I",contents,offset=pos)[0]
    pos = pos + 4
    logging.info("  Spatial dimensionality = {}".format(ndim))

    Nx = struct.unpack_from("@{}L".format(ndim),contents,offset=pos)
    pos = pos + 8*ndim
    M = 1
    for i in range(ndim):
      M = M * Nx[i]
    logging.info("  PW grid : {0}\tTotal # PWs = {1}".format(Nx,M))

    (kspacedata,complexdata) = struct.unpack_from("@2?",contents,offset=pos)
    pos = pos + 2
    logging.info("  k space? {}\tComplex container? {}".format(kspacedata,complexdata))
    if kspacedata:
      sys.stderr.write("\nError: k-space data is not supported\n")
      sys.exit(1)

    harray = struct.unpack_from("@{}d".format(ndim*ndim),contents,offset=pos)
    pos = pos + 8*ndim*ndim
    h = np.reshape(harray,(ndim,ndim))
    logging.info("  Cell tensor:\n{}".format(h))

    elsize = struct.unpack_from("@L",contents,offset=pos)[0]
    pos = pos + 8
    logging.info("  # bytes per element = {}".format(elsize))

    if elsize == 4 and not complexdata:
      logging.info("   * Single precision")
      fielddata = struct.unpack_from("@{}f".format(M*nfields),contents,offset=pos)
    elif elsize == 8 and complexdata:
      logging.info("   * Single precision")
      fielddata = struct.unpack_from("@{}f".format(2*M*nfields),contents,offset=pos)
    elif elsize == 8 and not complexdata:
      logging.info("   * Double precision")
      fielddata = struct.unpack_from("@{}d".format(M*nfields),contents,offset=pos)
    elif elsize == 16 and complexdata:
      logging.info("   * Double precision")
      fielddata = struct.unpack_from("@{}d".format(2*M*nfields),contents,offset=pos)
    else:
      sys.stderr.write("\nError: Unknown element size")
      sys.exit(1)

    # Now we need to return a numpy array with two indices: 
    #  1 = fieldidx (with imaginary part as a distinct field index)
    #  2 = PW indx
    # For complex fields, currently the real/imaginary part is the fastest index, 
    # so we will have to do a selective transpose
    # before reshaping the re/im sequence into the field indices
    if complexdata:
        AllFields = np.array(fielddata).reshape([nfields,M,2]).transpose((0,2,1)).reshape([nfields*2,M])
        nfields = nfields*2
    else:
        AllFields = np.array(fielddata).reshape([nfields,M])


    
    # calc all coords
    AllCoords = calcAllCoords(ndim,Nx,h,M)

    # Check whether the cell is orthorhombic
    orthorhombic = True
    for i in range(ndim):
        for j in range(ndim):
            if i == j:
                continue
            if abs(h[i][j]) > 1e-8:
                orthorhombic = False

    #return ndim, Nx, orthorhombic, M, nfields, AllCoords, AllFields

    AllCoords = np.reshape(AllCoords.T ,list(Nx) + [ndim])
    AllFields = np.reshape(AllFields.T ,list(Nx) + [nfields])
    return AllCoords, AllFields


def ReadDatFile(infilename):
    # Check whether the file is gzipped and handle that seamlessly
    ftest = open(infilename,'rb')
    twobytes=struct.unpack('2c',ftest.read(2)) # read 2 bytes and unpack as a tuple of chars
    ftest.close()
    if twobytes[0] == '\x1f' and twobytes[1] == '\x8b':
        logging.info( "  * File is gzipped")
        f = gzip.open(infilename,'r')
    else:
        logging.info( "  * File is not gzipped")
        f = open(infilename,'r')
    # Now parse the header to obtain the grid dimension and other format information
    version=int(f.readline().strip().split()[3])
    nfields=int(f.readline().strip().split()[3])
    ndim=int(f.readline().strip().split()[3])
    Nx=f.readline().strip().split()[4:4+ndim]
    Nx=[int(i) for i in Nx] # Convert all list entries to int
    kspacedata=True
    complexdata=True
    flagsline=f.readline().strip().split()
    if flagsline[4] == "0":
      kspacedata=False
    if flagsline[9] == "0":
      complexdata=False
    f.readline() # Skip the last header line
    logging.info("  * Format version = {}".format(version))
    logging.info("  * Number of fields = {}".format(nfields))
    logging.info("  * Number of spatial dimensions = {}".format(ndim))
    if ndim == 3:
      logging.info("  * Grid dimension = {}x{}x{}".format(*Nx))
    elif ndim == 2:
      logging.info("  * Grid dimension = {}x{}".format(*Nx))
    logging.info("  * K-space data? ",kspacedata)
    logging.info("  * Complex data? ",complexdata)
    if kspacedata:
      sys.stderr.write("\nError: k-space data is not supported\n")
      sys.exit(1)

    # Get the grid data but discard the grid coordinates
    AllData = np.loadtxt(f).transpose()
    AllCoords = AllData[:ndim]
    AllFields = AllData[ndim:]

    if complexdata:
        nfields = 2*nfields

    # Check whether the cell is orthorhombic
    orthorhombic = True
    if ndim == 2:
        if AllCoords[0][1] != 0.:
            orthorhombic = False
        if AllCoords[1][Nx[1]] != 0.:
            orthorhombic = False
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

    # this is old return status 
    #return ndim, Nx, np.prod(Nx), nfields, AllCoords, AllFields

    #change to return AllCoords and AllFields in a neat numpy nd array, makes ndim, Nx, ndim, nfields redundant
    AllCoords = np.reshape(AllCoords.T ,list(Nx) + [ndim])
    AllFields = np.reshape(AllFields.T ,list(Nx) + [nfields])
    return AllCoords, AllFields


def calcAllCoords(ndim,Nx,h,M):
    AllCoords = np.zeros((M,ndim))

    if ndim == 1:
      m=int(0)
      l=int(0)
      for ix in range(0,Nx[0]):
        l = ix
        xcart = h[0][0] * float(l) / Nx[0]

        AllCoords[m,:] = xcart
        m+=1

    elif ndim == 2:
      m=int(0)
      lm=[0,0]
      for iy in range(0,Nx[1]):
        lm[1] = iy
        for ix in range(0,Nx[0]):
          lm[0] = ix

          xfrac = [float(i) for i in lm]
          xfrac = np.divide(xfrac,Nx)

          xcart = [0., 0.]
          for i in range(ndim):
            for j in range(ndim):
              xcart[j] = xcart[j] + h[i][j]*xfrac[i]
          AllCoords[m,:] = xcart
          m+=1

    elif ndim == 3:
      m=int(0)
      lmn=[0,0,0]
      for ix in range(0,Nx[0]):
        lmn[0] = ix
        for iy in range(0,Nx[1]):
          lmn[1] = iy
          for iz in range(0,Nx[2]):
            lmn[2] = iz
            xfrac = [float(i) for i in lmn]
            xfrac = np.divide(xfrac,Nx)

            xcart = [0., 0., 0.]
            for i in range(ndim):
              for j in range(ndim):
                xcart[j] = xcart[j] + h[i][j]*xfrac[i]
            AllCoords[m,:] = xcart
            m+=1

    #np.savetxt("coords.dat",AllCoords)

    return AllCoords.transpose()


def WriteDatFile(outfilename, coords, fields, iskspace = False, iscomplex = False):

    outfilehndl = open(outfilename,'w')
    Dim = len(coords.shape) - 1
    nfields = fields.shape[Dim]

    assert(iscomplex==False), "fielddata isn't right if data includes complex, need to fix this"

    if iscomplex: 
        nfields = int(nfields/2)

    Nx = coords.shape[:Dim]
    if Dim == 1:
        griddim = (coords.shape[0],)
        hvoxel = np.array([coords[1]])
    elif Dim == 2:
        griddim = (coords.shape[0], coords.shape[1])
        hvoxel = np.array([coords[0,1],coords[1,0]])
    elif Dim == 3:
        griddim = (coords.shape[0], coords.shape[1],coords.shape[2])
        hvoxel = np.array([coords[0,0,1],coords[0,1,0],coords[1,0,0]])
    h = hvoxel * Nx
    kspacedata = iskspace
    complexdata = iscomplex
    M = np.prod(griddim)
    fielddata = np.ravel(fields,order='F')

    outfilehndl.write("# Format version 3\n")
    outfilehndl.write("# nfields = {0}\n".format(nfields))
    outfilehndl.write("# NDim = {0}\n".format(Dim))
    outfilehndl.write("# PW grid = ")
    for i in range(Dim):
        outfilehndl.write("{0} ".format(griddim[i]))
    outfilehndl.write("\n# k-space data = {0} , complex data = {1}\n".format(int(kspacedata),int(complexdata)))
    if Dim == 1:
      if kspacedata:
        outfilehndl.write("# Columns: kx indxx fielddata\n")
      else:
        outfilehndl.write("# Columns: x fielddata\n")

      for ix in range(0,griddim[0]):
        outfilehndl.write("%16.10g " % coords[ix][0])
        for n in range(nfields):
          outfilehndl.write("%16.10g " % fields[ix][n])
          outfilehndl.write("%16.10g " % fields[ix][n])
        outfilehndl.write("\n")

      #m=int(0)
      #l=int(0)
      #for ix in range(0,griddim[0]):
      #  l = ix
      #  if kspacedata and l > griddim[0]/2:
      #      l = l - griddim[0]
      #  if not kspacedata:
      #    xcart = h[0][0] * float(l) / griddim[0]
      #  else:
      #    xcart = h[0][0] * float(l)
      #  if kspacedata:
      #    outfilehndl.write("%7.4f %d " % (xcart,l))
      #  else:
      #    outfilehndl.write("%7.4f " % (xcart))

      #  if not complexdata and not kspacedata:
      #    for n in range(nfields):
      #      outfilehndl.write("%16.10g " % (fielddata[m + n*M]))
      #    outfilehndl.write("\n")
      #    m = m + 1
      #  else:
      #    for n in range(nfields):
      #      outfilehndl.write("%16.10g %16.10g " % (fielddata[m + n*M*2], fielddata[m + 1 + n*M*2]))
      #    outfilehndl.write("\n")
      #    m = m + 2

    elif Dim == 2:
      if kspacedata:
          outfilehndl.write("# Columns: kx ky indxx indxy fielddata\n")
      else:
          outfilehndl.write("# Columns: x y fielddata\n")

      for ix in range(0,griddim[0]):
          for iy in range(0,griddim[1]):
              outfilehndl.write("%7.4f %7.4f " % (coords[ix,iy,0], coords[ix,iy,1]))
              for n in range(nfields):
                  outfilehndl.write("%16.10g " % (fields[ix,iy,n]))
              outfilehndl.write("\n")
          outfilehndl.write("\n")

      #m=int(0)
      #lm=[0,0]
      #for ix in range(0,griddim[0]):
      #  lm[0] = ix
      #  if kspacedata and lm[0] > griddim[0]/2:
      #    lm[0] = lm[0] - griddim[0]
      #  for iy in range(0,griddim[1]):
      #    lm[1] = iy
      #    if kspacedata and lm[1] > griddim[1]/2:
      #      lm[1] = lm[1] - griddim[1]

      #    xfrac = [float(i) for i in lm]
      #    if not kspacedata:
      #      xfrac = np.divide(xfrac,griddim)

      #    xcart = [0., 0.]
      #    for i in range(Dim):
      #      for j in range(Dim):
      #        xcart[j] = xcart[j] + h[i][j]*xfrac[i]
      #    if kspacedata:
      #      outfilehndl.write("%7.4f %7.4f %d %d " % (xcart[0], xcart[1], lm[0], lm[1]))
      #    else:
      #      outfilehndl.write("%7.4f %7.4f " % (xcart[0], xcart[1]))

      #    if not complexdata and not kspacedata:
      #      for n in range(nfields):
      #        outfilehndl.write("%16.10g " % (fielddata[m + n*M]))
      #      outfilehndl.write("\n")
      #      m = m + 1
      #    else:
      #      for n in range(nfields):
      #        outfilehndl.write("%16.10g %16.10g " % (fielddata[m + n*M*2], fielddata[m + 1 + n*M*2]))
      #      outfilehndl.write("\n")
      #      m = m + 2
      #  outfilehndl.write("\n")

    elif Dim == 3:
      if kspacedata:
        outfilehndl.write("# Columns: kx ky kz indxx indxy indxz fielddata\n")
      else:
        outfilehndl.write("# Columns: x y z fielddata\n")

      for ix in range(0,griddim[0]):
        for iy in range(0,griddim[1]):
          for iz in range(0,griddim[2]):
              outfilehndl.write("%7.4f %7.4f %7.4f" % (coords[ix,iy,iz,0], coords[ix,iy,iz,1],coords[ix,iy,iz,2]))
              for n in range(nfields):
                  outfilehndl.write("%16.10g " % (fields[ix,iy,iz,n]))
              outfilehndl.write("\n")
          outfilehndl.write("\n")

      #m=int(0)
      #lmn=[0,0,0]
      #for ix in range(0,griddim[0]):
      #  lmn[0] = ix
      #  if kspacedata and lmn[0] > griddim[0]/2:
      #    lmn[0] = lmn[0] - griddim[0]
      #  for iy in range(0,griddim[1]):
      #    lmn[1] = iy
      #    if kspacedata and lmn[1] > griddim[1]/2:
      #      lmn[1] = lmn[1] - griddim[1]
      #    for iz in range(0,griddim[2]):
      #      lmn[2] = iz
      #      if kspacedata and lmn[2] > griddim[2]/2:
      #        lmn[2] = lmn[2] - griddim[2]

      #      xfrac = [float(i) for i in lmn]
      #      if not kspacedata:
      #        xfrac = np.divide(xfrac,griddim)

      #      xcart = [0., 0., 0.]
      #      for i in range(Dim):
      #        for j in range(Dim):
      #          xcart[j] = xcart[j] + h[i][j]*xfrac[i]
      #      if kspacedata:
      #        outfilehndl.write("%7.4f %7.4f %7.4f %d %d %d " % (xcart[0], xcart[1], xcart[2], lmn[0], lmn[1], lmn[2]))
      #      else:
      #        outfilehndl.write("%7.4f %7.4f %7.4f " % (xcart[0], xcart[1], xcart[2]))

      #      if not complexdata and not kspacedata:
      #        for n in range(nfields):
      #          outfilehndl.write("%16.10g " % (fielddata[m + n*M]))
      #        outfilehndl.write("\n")
      #        m = m + 1
      #      else:
      #        for n in range(nfields):
      #          outfilehndl.write("%16.10g %16.10g " % (fielddata[m + n*M*2], fielddata[m + 1 + n*M*2]))
      #        outfilehndl.write("\n")
      #        m = m + 2
      #  outfilehndl.write("\n")

    outfilehndl.close()


