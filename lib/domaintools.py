#!/usr/bin/env python3

import numpy as np
from skimage import measure
import pdb

'''
    Function to AnalyzeDomains that come out of a PolyFTS simulation
    
    - identify domains using the "burning alrgorithm"
        Based off of OperatorBridging in Josh's PolyFTS bridging code

    - mesh using marching cubes
    - calculate volume and area of domains using mesh

'''

# for burning algorithm
import sys
sys.setrecursionlimit(100000)

class DomainAnalyzer:
    def __init__(self,coords,fields, density_field_index=0, density_threshold = 0.5):
        ''' Define and calculate useful variables for DomainAnalysis routines
        '''

        self.__coords = coords
        self.__fields = fields

        self.__density_field_index = density_field_index
        self.__density_threshold = density_threshold

        self.__ndim = len(coords.shape) - 1
        self.__Nx = coords.shape[:3]
        self.__nfields = len(fields.shape) - self.__ndim
        self.__M = np.prod(self.__Nx)
        # assume box starts at (0,0,0) and ends at (lx,ly,lz)
        if not np.all(self.__coords.ravel()[0:self.__ndim] == np.array([0,0,0])):
            raise ValueError("coords[0,0,0] != (0,0,0)")
        self.__boxl = tuple(self.__coords[-1,-1,-1]) 
        self.__boxh = tuple(self.__coords[-1,-1,-1]*0.5) 
        self.__gridspacing = (self.__coords[1,0,0][0], self.__coords[0,1,0][1], self.__coords[0,0,1][2])

    def setDensityThreshold(density_threshold):
        self.__density_threshold = density_threshold
   
    def getDomainStats(self, plotMesh=False):
        ''' Calculate properties of each of the domains
            return com, surface_area, volume, IQ
        '''

        # create boolean selector from density fields for region definition
        isdomain_array = (self.__fields[:,:,:,self.__density_field_index] > self.__density_threshold)

        # FIXME, things break for non-cubic boxes. It must have to do with the vtk vs numpy indexing

        # identify domains
        self.__regionID = None # initially empty, created in computeRegionIDs
        self.__ndomains = self.identifyAndIndexDomains(isdomain_array)

        #nstats = 1+ 3*getCenter + getArea + getVol + getIQ
        #stats = np.zeros((self.__ndomains,nstats))
        com = np.zeros((self.__ndomains, self.__ndim))
        surface_area = np.zeros(self.__ndomains)
        volume = np.zeros(self.__ndomains)
        IQ = np.zeros(self.__ndomains)

        #for each domain
        for idomain in range(0,self.__ndomains):
            # calc center of domain
            com[idomain,:] = self.calcDomainCOM(idomain,units='coord')
            
            # mesh domain
            verts, faces, normals, values = self.meshDomain(idomain+1)

            # get surface area, volume and isoperimetric quotient
            surface_area[idomain] = measure.mesh_surface_area(verts, faces)

            volume[idomain] = self.mesh_volume(verts,faces)

            IQ[idomain] = self.calcIQ(surface_area[idomain], volume[idomain])
            
            if plotMesh: 
                self.plotMesh(verts,faces, "mesh.{}.png".format(idomain+1))

        return self.__ndomains, com, surface_area, volume, IQ
    
    def calcIQ(self, area, vol):
        return 36.0*np.pi * vol*vol / (area * area * area)

    def meshDomain(self,idomain):
        '''
        Function to:
        1) apply PBC to the domains so that an entire domain is continuous (ie not split across boundaries)
        2) Grab a little margin around each domain (the domain's "border") so that marching cubes can interpolate. The border is computed in identifyAndIndexDomains().
        3) Mesh the domain using marching cubes
        '''
        
        isdomain = (self.__regionID == idomain)
        isborder = (self.__borderID == idomain)

        com = self.calcDomainCOM(idomain,units='box')
        
        # center box and properties around center of mass (so that domains don't cross pbc)
        # np.roll is the key function here
        alldensity = self.__fields[:,:,:, self.__density_field_index]
        #coords_tmp = np.copy(self.__coords)
        for i in range(self.__ndim):
            shift = int(0.5*self.__Nx[i] - com[i])
            isdomain = np.roll(isdomain,shift,axis=i)
            isborder = np.roll(isborder,shift,axis=i)
            #coords_tmp = np.roll(coords_tmp,shift,axis=i)
            alldensity = np.roll(alldensity,shift,axis=i)

               
        # isolate the domain of interest
        isdomain_or_isborder = isdomain + isborder # since both bool, sum is the union of the two fields
        mydensity = np.zeros(self.__Nx)
        mydensity[isdomain_or_isborder] = alldensity[isdomain_or_isborder]

        # plot for debugging
        #import PolyFTS_to_VTK
        #AllCoords = np.reshape(coords_tmp,(self.__M, self.__ndim))
        #AllCoords = AllCoords.T
        #tmp = np.ravel(isdomain)
        #tmp = np.resize(tmp,(1,len(tmp)))
        #PolyFTS_to_VTK.writeVTK("isdomain.vtk", self.__Nx, True, self.__M, AllCoords,tmp)
        #tmp = np.ravel(isborder)
        #tmp = np.resize(tmp,(1,len(tmp)))
        #PolyFTS_to_VTK.writeVTK("isborder.vtk", self.__Nx, True, self.__M, AllCoords,tmp)
        #tmp = np.ravel(mydensity)
        #tmp = np.resize(tmp,(1,len(tmp)))
        #PolyFTS_to_VTK.writeVTK("mydensity.vtk", self.__Nx, True, self.__M, AllCoords,tmp)

        # mesh! (using scikit-image)
        if self.__ndim == 2:
            raise NotImplementedError("Meshing in 2 dimensions is in development")
            contours = skimage.measure.find_contours(mydensity, self.__density_threshold) 
        elif self.__ndim == 3:
            #from skimage import measure
            verts, faces, normals, values = measure.marching_cubes_lewiner(mydensity, self.__density_threshold, spacing = self.__gridspacing)
            return verts, faces, normals, values
        else:
            raise ValueError("Meshing makes no sense in 1 dimension!")


    def mesh_volume(self, verts, faces):
        actual_verts = verts[faces]
        v0 = actual_verts[:,0,:]
        v1 = actual_verts[:,1,:]
        v2 = actual_verts[:,2,:]
        # 1/6 \sum v0 \cdot (v1 x v2)
        return 1.0/6.0 * np.abs( (v0*np.cross(v1,v2)).sum(axis=1).sum() )
 
    def plotMesh(self, verts, faces, filename=None):
        ''' Plot a mesh from marching cubes
        '''
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        # Display resulting triangular mesh using Matplotlib. This can also be done
        # with mayavi (see skimage.measure.marching_cubes_lewiner docstring).
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Fancy indexing: `verts[faces]` to generate a collection of triangles
        mesh = Poly3DCollection(verts[faces])
        mesh.set_edgecolor('k')
        ax.add_collection3d(mesh)

        ax.set_xlim(0, self.__boxl[0])  
        ax.set_ylim(0, self.__boxl[1])  
        ax.set_zlim(0, self.__boxl[2])  

        plt.tight_layout()
        if not filename:
            plt.show()
        else:
            plt.savefig(filename)

        plt.close()

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
        self.__borderID = np.full(self.__Nx,0,dtype=np.int32);
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
              # note - dont assign borders here, this is acomplished inside of spread_region()
              self.__regionID[i] = 0;
              self.__image_flags[i]= np.zeros(self.__ndim)

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
                  

            if self.__regionID[neighbor] == 0:
              # must have neighbors that are domain (since spread region is only called 
              #   if coord_center is a domain). Therefore, it's a border
              if self.__borderID[neighbor] != 0 and self.__borderID[neighbor] != region_number:
                #raise RuntimeError("Trying to set borderID[{0}]={1} but is it is already set to {2}".format(neighbor, region_number, self.__borderID[neighbor]))
                print("Trying to set borderID[{0}]={1} but is it is already set to {2}".format(neighbor, region_number, self.__borderID[neighbor]))
              self.__borderID[neighbor] = region_number
 
              # set image flags of non-domain adjacent to domain according to the domain
              # basically, I need the border to have the correct image flags

              # this will override whatever image_flag was set before, generally this shouldnt be a 
              #   problem unless two domains are very close to each other and their borders overlap
              #   Raise and error just in case...
              #   make sure, not (0,0,0) and not new == old
              if not np.all(self.__image_flags[neighbor]==0) and not np.all(self.__image_flags[neighbor] == image_flag):
                #raise RuntimeError("Tried to set image_flag of {0} to {2} but it was already set \
                #    to {1} and not (0,0,0)".format(neighbor,self.__image_flags[neighbor], image_flag))
                print("Tried to set image_flag of {0} to {2} but it was already set \
                    to {1} and not (0,0,0)".format(neighbor,self.__image_flags[neighbor], image_flag))
              
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
            
