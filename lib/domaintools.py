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
        self.__Nx = coords.shape[:self.__ndim]
        self.__nfields = len(fields.shape) - self.__ndim
        self.__M = np.prod(self.__Nx)


        # assume box starts at (0,0,0) and ends at (lx,ly,lz)
        if not np.all(self.__coords.ravel()[0:self.__ndim] == np.zeros(self.__ndim)):
            raise ValueError("coords[0,0,0] != (0,0,0)")

        if self.__ndim == 2:
            self.__boxl = tuple(self.__coords[-1,-1]) 
            self.__boxh = tuple(self.__coords[-1,-1]*0.5) 
            self.__gridspacing = (self.__coords[1,0][0], self.__coords[0,1][1])
            self.__hvoxel = np.array([coords[1,0],coords[0,1]])
        elif self.__ndim == 3:
            self.__boxl = tuple(self.__coords[-1,-1,-1]) 
            self.__boxh = tuple(self.__coords[-1,-1,-1]*0.5) 
            self.__gridspacing = (self.__coords[1,0,0][0], self.__coords[0,1,0][1], self.__coords[0,0,1][2])
            self.__hvoxel = np.array([coords[1,0,0],coords[0,1,0],coords[0,0,1]])
        self.__hcell = self.__hvoxel * self.__Nx
        self.__volvoxel = np.linalg.det(self.__hvoxel)
        assert (np.abs(self.__volvoxel - np.linalg.det(self.__hcell) / self.__M) < 1e-5), "Volume of voxel != (Volume of cell / n voxels). This should be true!"

        # check if orthorhombic
        self.__orthorhombic = True
        hnorm = self.__hcell /  np.linalg.norm(self.__hcell, axis=0)
        if self.__ndim == 2 and np.dot(hnorm[0],hnorm[1]) != 0:
            self.__orthorhombic = False
        elif self.__ndim == 3 :
            if np.dot(hnorm[0],[1,0,0]) != 0 or np.dot(hnorm[1],[0,1,0]) != 0 or np.dot(hnorm[2],[0,0,1]) != 0:
                self.__orthorhombic = False

        # check if density field is reasonable between 0-1, if not throw warning
        if self.__ndim == 2:
            mindensity= np.min(self.__fields[:,:,self.__density_field_index])
            maxdensity= np.max(self.__fields[:,:,self.__density_field_index])
        elif self.__ndim == 3:
            mindensity= np.min(self.__fields[:,:,:,self.__density_field_index])
            maxdensity= np.max(self.__fields[:,:,:,self.__density_field_index])
        if maxdensity > 1.0 or mindensity < 0.0:
            print("Warning: The density field is not between 0-1 (min: {}, max: {}). The specified threshold of {} might be inappropriate.".format(mindensity,maxdensity,self.__density_threshold))

        self.__needToIndexDomains = True

    def setDensityThreshold(density_threshold):
        self.__density_threshold = density_threshold
        # if changing the Density threshold, will need to index domains again
        self.__needToIndexDomains = True
  
    def getNdim():
        return self.__ndim

    def getDomainStats(self, useMesh=True, plotMesh=False):
        ''' Calculate properties of each of the domains
            return com, surface_area, volume, IQ

            if useMesh == True, calculate a isosurface mesh to calculate the volumes and areas. 
                This is very accurate, but can have issues creating a good mesh if domains are poorly defined (as in certain CL systems)

                (Specifically the issue is if two domains are only separated by a single grid point.  When this happens, 
                 the border around the domain belongs to two domains simultaneously and my current burning algorithm throws 
                 an error. I use the border around a domain when applying PBC's to make sure a domain is continuous. 
                 Eventually I might think of a better algorithm that will be robust to this edge case...
                )

            useMesh == False uses the less accurate approach of summing over the voxels to get the volume and area
               the volume is still pretty accurate, the area...well, I'm not even going to implement it since in CL I only want volume

        '''

        if useMesh and not self.__orthorhombic: 
            print("Warning: computing volume using mesh, but cell is not orthorhombic. This will lead to errors in the surface areas calculation of the domains")


        # create boolean selector from density fields for region definition
        if self.__ndim == 2:
            isdomain_array = (self.__fields[:,:,self.__density_field_index] > self.__density_threshold)
        elif self.__ndim == 3:
            isdomain_array = (self.__fields[:,:,:,self.__density_field_index] > self.__density_threshold)

        # FIXME, things break for non-cubic boxes. It must have to do with the vtk vs numpy indexing

        # identify domains
        if self.__needToIndexDomains:
            self.__regionID = None # initially empty, created in computeRegionIDs
            self.__ndomains = self.identifyAndIndexDomains(isdomain_array)
        else:
            print("Note: Using cached domain ID's")

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
            
            if useMesh:

                if self.__ndim == 2:
                    # mesh domain
                    contours = self.meshSingleDomain(idomain+1)
                    assert (len(contours) == 1), "The contour should only be one curve, if not the area and volume calculations will be completely wrong!"

                    # get surface area (perimeter) and volume (area)
                    surface_area[idomain] = self.contour_perimeter(contours[0])
                    volume[idomain] = self.contour_area(contours[0])

                    if plotMesh: 
                        self.plotContours2D(contours,filename="mesh.{}.png".format(idomain+1))

                if self.__ndim == 3: 
                    # mesh domain
                    verts, faces, normals, values = self.meshSingleDomain(idomain+1)

                    # get surface area, volume and isoperimetric quotient
                    surface_area[idomain] = measure.mesh_surface_area(verts, faces)
                    volume[idomain] = self.mesh_volume(verts,faces)

                    if plotMesh: 
                        self.plotMesh3D(verts,faces, filename="mesh.{}.png".format(idomain+1))

                IQ[idomain] = self.calcIQ(surface_area[idomain], volume[idomain])

            else:
                surface_area[idomain] = -1.0 #FIXME surface_area is currently not calculated if no mesh
                volume[idomain] = self.voxel_volume(idomain+1) # get volume from voxels
                IQ[idomain] = 0.0

        return self.__ndomains, com, surface_area, volume, IQ
    
    def calcIQ(self, area, vol):
        '''returns isoperimetric coefficient. 1 for perfect circle or sphere, less for other shapes
           note that in 2d "area" is actually perimeter, and "vol" is actually area
           This difference didn't seem to warrant a completely different method though
        '''
        if self.__ndim == 2:
            return 4.0*np.pi*vol / (area * area)
        elif self.__ndim == 3:
            return 36.0*np.pi * vol*vol / (area * area * area)

    def meshAllDomains(self,datafile=None,plotfile=None):
        ''' Mesh all domains using marching cubes or marching squares
            Options:
            - Save plot of mesh to plotfile if specified
            - save mesh data to file if specified


        '''
        if self.__ndim == 2:
            mydensity = self.__fields[:,:, self.__density_field_index]
            contours = measure.find_contours(mydensity, self.__density_threshold) 
            if datafile:
                self.writeContours(contours,datafile)
            if plotfile:
                self.plotContours2D(contours,surface=mydensity,filename=plotfile)

        elif self.__ndim == 3:
            mydensity = self.__fields[:,:,:, self.__density_field_index]
            verts, faces, normals, values = measure.marching_cubes_lewiner(mydensity, self.__density_threshold, spacing = self.__gridspacing)
            if datafile:
                raise NotImplementedError("Support for writing 3D mesh not implemented")
            if filename:
                self.plotMesh3D(verts,faces,filename=filename)
            return verts,faces, normals, values


    def meshSingleDomain(self,idomain):
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
        if self.__ndim == 2:
            alldensity = self.__fields[:,:, self.__density_field_index]
        elif self.__ndim == 3:
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
            #raise NotImplementedError("Meshing in 2 dimensions is in development")
            contours = measure.find_contours(mydensity, self.__density_threshold) 
            return contours
        elif self.__ndim == 3:
            #from skimage import measure
            verts, faces, normals, values = measure.marching_cubes_lewiner(mydensity, self.__density_threshold, spacing = self.__gridspacing)
            return verts, faces, normals, values
        else:
            raise ValueError("Meshing makes no sense in 1 dimension!")

    def contour_perimeter(self,contour):
        '''calculate perimeter of contour by suming up the line-segment lengths
        '''
        #TODO vectorize this for loop
        p = 0.0
        n=contour.shape[0]
        for i in range(n-1):
           v = contour[i+1] - contour[i] 
           p += np.square(v).sum()
        return p

    def contour_area(self,contour):
        ''' Calculate area of shape enclosed in contour
            similar to calculating mesh volume
            use trick from http://geomalgorithms.com/a01-_area.html
        '''
        assert (np.all(contour[0] == contour[-1])), "Contour must be closed! (1st point == last point)"
        
        #TODO vectorize this for loop
        area = 0.0
        n=contour.shape[0]
        for i in range(n-1):
            area += np.cross(contour[i],contour[i+1])
        return 0.5*np.abs(area)


    def mesh_volume(self, verts, faces):
        '''calculate volume of a mesh, using cross product trick
        '''
        actual_verts = verts[faces]
        v0 = actual_verts[:,0,:]
        v1 = actual_verts[:,1,:]
        v2 = actual_verts[:,2,:]
       
        # TODO: dont do the volume rescaling here, instead change the actual position of "verts" in getDomainStats my scaling each vert position by h (or something along these lines)

        # introduce factor to scale the volume if non-orthorhombic box
        # this is because the mesh is generated assuming a
        if self.__orthorhombic:
            factor=1.0
        else:
            factor = self.__volvoxel / np.prod(self.__gridspacing)

        # 1/6 \sum v0 \cdot (v1 x v2)
        return factor * 1.0/6.0 * np.abs( (v0*np.cross(v1,v2)).sum(axis=1).sum() )

    def voxel_volume(self,idomain):
        ''' Get volume of idomain using voxels
        '''
        #v_voxel = np.prod(self.__gridspacing) # volume of single voxel
        v_voxel = self.__volvoxel
        n_voxel = np.sum(self.__regionID == idomain) # number of voxels in ith domain
        return v_voxel*n_voxel
 
    def writeContours(self, contours,filename):
        ''' write contours to data files
            The format is built for using the gnuplot command "plot 'file' index 0 u 1:2"
            Each individual contor is plotted in two x,y columns
            Each contour is separated by two new lines (see gnuplot "index" for explanation)
        '''
        with open(filename,'wb') as f:
            f.write(b"# NContours = %d\n" % len(contours))
            for contour in contours:
                #np.savetxt(f,contour,footer='\n',comments='')
                np.savetxt(f,contour)
                f.write(b"\n\n")

    def plotContours2D(self, contours, surface=None, filename=None):
        ''' Plot a mesh from marching squares
        '''
        import matplotlib.pyplot as plt

        # Display the image and plot all contours found
        fig, ax = plt.subplots()
        if surface is not None:
            #ax.imshow(surface.T, interpolation='nearest', cmap=plt.cm.gray)
            ax.imshow(surface.T, interpolation='nearest')

        for n, contour in enumerate(contours):
            ax.plot(contour[:, 0], contour[:, 1], linewidth=2)

        #ax.axis('image')
        #ax.set_xticks([])
        #ax.set_yticks([])
        if not filename:
            plt.show()
        else:
            plt.savefig(filename)
        plt.close()

    def plotMesh3D(self, verts, faces, filename=None):
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

        self.__needToIndexDomains = False
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
                raise RuntimeError("Trying to set borderID[{0}]={1} but is it is already set to {2}".format(neighbor, region_number, self.__borderID[neighbor]))
                #print("Trying to set borderID[{0}]={1} but is it is already set to {2}".format(neighbor, region_number, self.__borderID[neighbor]))
              self.__borderID[neighbor] = region_number
 
              # set image flags of non-domain adjacent to domain according to the domain
              # basically, I need the border to have the correct image flags

              # this will override whatever image_flag was set before, generally this shouldnt be a 
              #   problem unless two domains are very close to each other and their borders overlap
              #   Raise and error just in case...
              #   make sure, not (0,0,0) and not new == old
              if not np.all(self.__image_flags[neighbor]==0) and not np.all(self.__image_flags[neighbor] == image_flag):
                raise RuntimeError("Tried to set image_flag of {0} to {2} but it was already set \
                    to {1} and not (0,0,0)".format(neighbor,self.__image_flags[neighbor], image_flag))
                #print("Tried to set image_flag of {0} to {2} but it was already set \
                #    to {1} and not (0,0,0)".format(neighbor,self.__image_flags[neighbor], image_flag))
              
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
            
