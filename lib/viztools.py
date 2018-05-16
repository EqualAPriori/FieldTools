#!/usr/bin/env python3

import sys
import numpy as np
from os import path
import logging
import re
import pdb



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
    ndim = len(Nx)
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

def isOrthorhombic(AllCoords):
    ndim = len(AllCoords.shape) - 1
    h = np.zeros((ndim, ndim))
    if ndim == 2:
        if AllCoords[0,1][1] != 0 or AllCoords[1,0][0] != 0:
            return False
    if ndim == 3:
        if AllCoords[1,0,0][1] != 0 or AllCoords[1,0,0][2] != 0:
            return False
        if AllCoords[0,1,0][0] != 0 or AllCoords[0,1,0][2] != 0:
            return False
        if AllCoords[0,0,1][0] != 0 or AllCoords[0,0,1][1] != 0:
            return False
    return True



#def writeVTK(outfile, Nx, orthorhombic, M, AllCoords, AllFields):
def writeVTK(outfile, AllCoords, AllFields,binary=False):

    ndim = len(AllCoords.shape) - 1
    Nx = AllCoords.shape[:ndim]
    orthorhombic = isOrthorhombic(AllCoords)
    M = np.prod(Nx)
    nfields = AllFields.shape[-1]

    # get mesh spacing
    if orthorhombic == True:
        if ndim == 1:
            spacing = (AllCoords[1][0])
        if ndim == 2:
            spacing = (AllCoords[1,0][0], AllCoords[0,1][1])
        if ndim == 3:
            spacing = (AllCoords[1,0,0][0], AllCoords[0,1,0][1], AllCoords[0,0,1][2])
        logging.info("Mesh spacing       {}".format(spacing))
    
    #FIXME if nfields is only 1 then reshape
    #if nfields == 1:
    #    AllFields = np.reshape(AllFields,(len(AllFields),1))

    if binary:
        raise NotImplementedError("Implementation started but not tested and complete...")

    AllCoords = np.ravel(AllCoords)
    AllCoords = np.reshape(AllCoords,(M,ndim ))
    AllFields = np.ravel(AllFields)
    AllFields = np.reshape(AllFields,(M,nfields))

    # if 2D, set 3rd dimension to zero
    if ndim==2: 
        AllCoords = np.hstack((AllCoords,np.zeros((M,1))))


    # dont need to do this anymore? can just used np to ravel and "order" option?
    ## Generate the mapping from PolyFTS (z fastest in 3D) to VTK (x fastest)
    IdxMap = get_PolyFTS_to_VTK_IdxMap(M,Nx)
    # Remap field samples from PolyFTS order to VTK order
    AllCoords = AllCoords[IdxMap,:]
    AllFields = AllFields[IdxMap,:]

    AllCoords = AllCoords.T
    AllFields = AllFields.T

    # Write Legacy VTK file.
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
            o.close();o = open(outfile,"ab")
            np.savetxt(o, AllFields[i], fmt="%.11f")
            o.close();o = open(outfile,"a")
    else:
      o.write("# vtk DataFile Version 3.0\n")
      o.write("PolyFTS field data\n")
      if binary:
          o.write("BINARY\n")
      else:
          o.write("ASCII\n")
      o.write("DATASET STRUCTURED_GRID\n")
      if ndim == 1:
          o.write("DIMENSIONS {} 1 1\n".format(*Nx))
      elif ndim == 2:
          o.write("DIMENSIONS {} {} 1\n".format(*Nx))
      elif ndim == 3:
          o.write("DIMENSIONS {} {} {}\n".format(*Nx))
      o.write("POINTS {} double\n".format(M))
      # Remap the mesh coordinates to VTK order
      AllCoords = AllCoords[:,IdxMap]
      #np.savetxt('coords.dat',AllCoords.transpose()) # Debug
      o.close(); o = open(outfile,'ab')
      if binary:
          AllCoords.transpose().tofile(o)
      else:
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
          o.write("VECTORS field{0} double\n".format(i)) #writing as vector fixed the seg fault
          tmp=np.zeros((3,M))
          tmp[0,:] = AllFields[i]
          o.close(); o = open(outfile,'ab')
          if binary:
              tmp.transpose().tofile(o)
          else:
              np.savetxt(o, tmp.transpose(), fmt="%14.11f")
          o.close(); o = open(outfile,'a')
    o.close()
  
def writeCSV(outfile, AllCoords, AllFields):
    ndim = len(AllCoords.shape) - 1
    Nx = AllCoords.shape[:ndim]
    orthorhombic = isOrthorhombic(AllCoords)
    M = np.prod(Nx)
    nfields = AllFields.shape[-1]

    AllCoords = np.ravel(AllCoords)
    AllCoords = np.reshape(AllCoords,(M,ndim ))
    AllFields = np.ravel(AllFields)
    AllFields = np.reshape(AllFields,(M,nfields))

    o = open(outfile,"w")
    o.write("x coord, y coord, z coord, scalar\n")
    for i in range(M):
        o.write("{}, {}, {}, {}\n".format(AllCoords[i,0],AllCoords[i,1],AllCoords[i,2],AllFields[i,0]))
    o.close()
       
    


def writeVTK_XML(outfile, AllCoords, AllFields):

    ndim = len(AllCoords.shape) - 1
    Nx = AllCoords.shape[:ndim]
    orthorhombic = isOrthorhombic(AllCoords)
    M = np.prod(Nx)
    nfields = AllFields.shape[-1]
    x1=AllCoords[0,0,0,0]
    x2=AllCoords[-1,0,0,0] 
    y1=AllCoords[0,0,0,1] 
    y2=AllCoords[0,-1,0,1] 
    z1=AllCoords[0,0,0,2] 
    z2=AllCoords[0,0,-1,2] 


    #FIXME if nfields is only 1 then reshape
    #if nfields == 1:
    #    AllFields = np.reshape(AllFields,(len(AllFields),1))

    AllCoords = np.ravel(AllCoords)
    AllCoords = np.reshape(AllCoords,(M,ndim ))
    AllFields = np.ravel(AllFields)
    AllFields = np.reshape(AllFields,(M,nfields))

    # dont need to do this anymore? can just used np to ravel and "order" option?
    ## Generate the mapping from PolyFTS (z fastest in 3D) to VTK (x fastest)
    IdxMap = get_PolyFTS_to_VTK_IdxMap(M,Nx)
    # Remap field samples from PolyFTS order to VTK order
    AllCoords = AllCoords[IdxMap,:]
    AllFields = AllFields[IdxMap,:]

    #AllCoords = AllCoords.T
    #AllFields = AllFields.T


    o = open(outfile,"w")
    o.write("<VTKFile type='”RectilinearGrid' version='0.1' byte_order='LittleEndian'>\n")
    o.write("<RectilinearGrid WholeExtent=\" {} {} {} {} {} {}\"> \n".format(x1,x2,y1,y2,z1,z2))
    o.write("<Piece Extent=' {} {} {} {} {} {}'>\n".format(x1,x2,y1,y2,z1,z2))
    o.write("<PointData Scalars=\"field0\">\n") # fixme for multiple fields
    for i in range(M):
        o.write("<DataArray Name='field0' type='Float32' NumberOfComponents='1' >") # fixme for multiple fields
        o.write("{:f}".format(AllFields[i][0])) #fixme for multiple fields
        o.write("</DataArray>\n")
    o.write("</PointData>\n")
 
    o.write("<Coordinates>\n")
    for i in range(M):
         for j in range(3):
             o.write("<DataArray type='Float32' NumberOfComponents='1' >") # fixme for multiple fields
             o.write("{}".format(AllCoords[i][j]))
             o.write("</DataArray>\n")
    o.write("</Coordinates>\n")
 
    o.write("</Piece>\n")
    o.write("</RectilinearGrid>\n")
    o.write("</VTKFile>\n")
    o.close() 

