#!/usr/bin/env python3

import numpy as np
import pdb

'''
    Function to AnalyzeDomains that come out of a PolyFTS simulation
    
    - identify domains using the "burning alrgorithm"
        Based off of OperatorBridging in Josh's PolyFTS bridging code

    - mesh using marching cubes
    - calculate volume and area of domains

'''
import sys
sys.setrecursionlimit(100000)

class DomainAnalyzer:
    def __init__(self,coords,fields):
        self.__coords = coords
        self.__fields = fields
        self.__ndim = len(coords.shape) - 1
        self.__Nx = coords.shape[:3]
        self.__nfields = len(fields.shape) - self.__ndim
        self.__M = np.prod(self.__Nx)
        # assume box starts at (0,0,0) and ends at (lx,ly,lz)
        if not np.all(self.__coords.ravel()[0:self.__ndim] == np.array([0,0,0])):
            raise ValueError("coords[0,0,0] != (0,0,0)")
        self.__boxl = tuple(self.__coords[-1,-1,-1]) 
        self.__boxh = tuple(self.__coords[-1,-1,-1]*0.5) 


        density_threshold = 0.5
        # create boolean selector from density fields for region definition
        isdomain_array = (fields[:,:,:,0] > density_threshold)

        # FIXME, things break for non-cubic boxes. It must have to do with the vtk vs numpy indexing

        import PolyFTS_to_VTK
        AllCoords = np.reshape(coords,(self.__M, self.__ndim))
        AllCoords = AllCoords.T

        # identify domains
        self.__regionID = None # initially empty, created in computeRegionIDs
        self.__nregions = self.identifyAndIndexDomains(isdomain_array)
        

        #tmp = np.ravel(self.__image_flags)
        #tmp = np.resize(tmp,(len(tmp),3))
        #tmp = tmp.T
        #PolyFTS_to_VTK.writeVTK("image_flags.vtk", self.__Nx, True, self.__M, AllCoords, tmp)


        
        
        tmp = np.ravel(self.__isborder)
        tmp = np.resize(tmp,(1,len(tmp)))
        PolyFTS_to_VTK.writeVTK("isborder.vtk", self.__Nx, True, self.__M, AllCoords,tmp)

        domain_ids = np.ravel(self.__regionID)
        domain_ids = np.resize(domain_ids,(1,len(domain_ids)))
        PolyFTS_to_VTK.writeVTK("domains.vtk", self.__Nx, True, self.__M, AllCoords, domain_ids)

        pdb.set_trace()
    
        self.calcDomainCOM(1)
        self.encloseDomainInBox(1)

        #from skimage import measure
        #verts, faces, normals, values = measure.marching_cubes_lewiner(isdomain_array, 0)

        #import matplotlib.pyplot as plt
        #from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        #fig = plt.figure(figsize=(10, 10))
        #ax = fig.add_subplot(111, projection='3d')

        ### Fancy indexing: `verts[faces]` to generate a collection of triangles
        #mesh = Poly3DCollection(verts[faces])
        #mesh.set_edgecolor('k')
        #ax.add_collection3d(mesh)

        #ax.set_xlim(0, 32)  # a = 6 (times two for 2nd ellipsoid)
        #ax.set_ylim(0, 32)  # b = 10
        #ax.set_zlim(0, 64)  # c = 16

        #plt.tight_layout()
        #plt.show()





    def getRegionIDs(self):
        return self.__regionID, self.__numRegionNeighbors, self.__nregions

    
    def calcDomainCOM(self,idomain, units='box'):
        ''' given a domain index, apply PBC and return the center of mass
            Can return result in 'box' units (0 to Nx) or in 'coord' units (0 to boxl)
        '''

        isdomain = (self.__regionID == idomain)
        N = np.sum(isdomain)
        indicies = np.transpose(np.nonzero(isdomain))
        coords = np.zeros((N,self.__ndim))

        
        #TODO could I do this without for loop? (will be faster)
        for i in range(N):
            index = tuple(indicies[i])

            if units == "box":
               coord = index + self.__image_flags[index] * self.__Nx
            elif units == "coord":
               coord = self.__coords[index] + self.__image_flags[index] * self.__boxl
            else:
                raise ValueError("Invalid units entry of \'%s\'" % units)
            
            coords[i] = coord
        
        # now average in order to get center of the domain (each point weighted evenly)
        return np.average(coords,axis=0)

    def meshDomain(self,idomain):
        '''
        Function to:
        1) apply PBC to the domains so that an entire domain is continuous (ie not split across boundaries)
        2) Grab a little margin around each domain so that marching cubes can interpolate
        3) Mesh the domain using marching cubes
        '''
        

        com = calcDomainCOM(idomain,units='box')

        # get 



    def encloseDomainInBox(self,idomain):
        ''' 
        determine minimal orthorhombic box that encloses a domain 
        return value: a boolean mask of shape self.__Nx
        '''

        is_target_domain = (self.__regionID == idomain)



    def identifyAndIndexDomains(self, isdomain_array):
        ''' This function populates the regionID member variable
            if regionID == 0, it is the continuous domain
            points with regionID == i, correspond to the ith domain
            Also sets - image_flags (which PBC a domain belongs to) and 
                      - isborder (whether a grid is adjacent to a domain)
        '''
        # if regionID == -1, it has not been visited
        #M = np.size(isdomain_array)
        self.__regionID = np.full(self.__Nx,-1,dtype=np.int32);
        self.__isborder = np.full(self.__Nx,False,dtype=np.bool);
        self.__image_flags = np.zeros(list(self.__Nx) + [self.__ndim])
        

        region_number = 1;

        #this is where the recursive magic happens
        for i in np.ndindex(self.__Nx):
          if (self.__regionID[i] == -1):
            if (isdomain_array[i]):
              current_image_flag = np.zeros(self.__ndim)
              self.spread_region(i, region_number, isdomain_array,current_image_flag);
              region_number += 1;
            else:
              self.__regionID[i] = 0;
              self.__image_flags[i]= np.zeros(self.__ndim)

              # get neighbors to determine if this is a border
              neighbors, neigh_image_flags = self.getNeighbors(i)
              nneigh_in_domain = sum([isdomain_array[n] for n in neighbors])
              if nneigh_in_domain >= 1:
                  self.__isborder[i] = True

        nregions = region_number-1;

        return nregions
        

    def spread_region(self, coord_center, region_number, isdomain_array,current_image_flag):
        ''' recursive function:
            given a point, find the neighbors of that point, 
            for each neighbor, send back into function
        '''

        self.__regionID[coord_center] = region_number;
        self.__image_flags[coord_center] = current_image_flag

        neighbors,neigh_image_flags = self.getNeighbors(coord_center, current_image_flag);

        for i in range(len(neighbors)):
            neighbor = neighbors[i]
            image_flag = neigh_image_flags[i]
            if (self.__regionID[neighbor] == -1):
                if (isdomain_array[neighbor]):
                  self.spread_region(neighbor, region_number, isdomain_array, image_flag);
                else:
                  self.__regionID[neighbor] = 0;
                  # must have neighbors that are domain (since spread region is only called 
                  #   if coord_center is a domain). Therefore, it's a border
                  self.__isborder[neighbor] = True
                  

            if self.__regionID[neighbor] == 0:
              # set image flags of non-domain adjacent to domain according to the domain
              # basically, I need the border to have the correct image flags

              # this will override whatever image_flag was set before, generally this shouldnt be a 
              #   problem unless two domains are very close to each other and their borders overlap
              #   Raise and error just in case...
              #   make sure, not (0,0,0) and not new == old
              if not np.all(self.__image_flags[neighbor]==0) and not np.all(self.__image_flags[neighbor] == image_flag):
                raise RuntimeError("Tried to set image_flag of {0} to {2} but it was already set to {1} and not (0,0,0)".format(neighbor,self.__image_flags[neighbor], image_flag))
              
              self.__image_flags[neighbor] = image_flag
        
    def getNeighbors(self,coord_center,center_image_flag=[]):
       ''' given a coord (tuple), return 
            1) the neighbors of that coord (also tuple) AND 
            2) the image_flag (which PBC) that neighbor corresponds to
       '''

       # set default
       if center_image_flag == []:
            center_image_flag = np.zeros(self.__ndim)

       neighbors = [];
       neigh_image_flags = np.tile(center_image_flag, (2*self.__ndim,1))
       for i in range(self.__ndim):
          coord_neigh = np.copy(coord_center)
          coord_neigh[i] -= 1;
          self.applyPBC(coord_neigh, neigh_image_flags[2*i]);
          neighbors.append(tuple(coord_neigh))

          coord_neigh = np.copy(coord_center)
          coord_neigh[i] += 1
          self.applyPBC(coord_neigh,neigh_image_flags[2*i+1])
          neighbors.append(tuple(coord_neigh))

       return neighbors, neigh_image_flags

    def applyPBC(self,coord,image_flag):
      for i in range(self.__ndim):
         if coord[i] >= self.__Nx[i]: 
             coord[i] = 0
             image_flag[i] += 1
         if coord[i] < 0:             
             coord[i] = self.__Nx[i] - 1
             image_flag[i] -= 1
            

# testing code if executed directly
if __name__ == "__main__":
    pass
