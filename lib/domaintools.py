#!/usr/bin/env python3
'''
    Joshua Lequieu <lequieu@mrl.ucsb.edu>
'''

import numpy as np
from skimage import measure
import pdb
#import viztools as viz

# Need to increase recursion limit for burning algorithm
import sys
sys.setrecursionlimit(800000)
# 

class DomainAnalyzer:
    '''
        Class to AnalyzeDomains that come out of a PolyFTS simulation
        
        - identify domains using the "burning alrgorithm"
            Based off of OperatorBridging in Josh's PolyFTS bridging code

        - mesh using marching cubes
        - calculate volume and area of domains using mesh

    '''


    def __init__(self,coords,fields, density_field_index=0, density_threshold = 0.5):
        ''' Define and calculate useful variables for DomainAnalysis routines
        '''

        self.__coords = coords
        self.__fields = fields

        self.__density_field_index = density_field_index
        self.__density_threshold = density_threshold

        self.__ndim = len(coords.shape) - 1
        self.__Nx = coords.shape[:self.__ndim]
        self.__nfields = fields.shape[self.__ndim]
        self.__M = np.prod(self.__Nx)


        # assume box starts at (0,0,0) and ends at (lx,ly,lz)
        if not np.all(self.__coords.ravel()[0:self.__ndim] == np.zeros(self.__ndim)):
            raise ValueError("coords[0,0,0] != (0,0,0)")

        if self.__ndim == 2:
            #self.__boxl = tuple(self.__coords[-1,-1] + self.__coords[1,1]) 
            #self.__boxh = tuple(np.array(self.__boxl)*0.5) 
            self.__gridspacing = (self.__coords[1,0][0], self.__coords[0,1][1])
            self.__hvoxel = np.array([coords[1,0],coords[0,1]])
        elif self.__ndim == 3:
            #self.__boxl = tuple(self.__coords[-1,-1,-1] + self.__coords[1,1,1]) 
            #self.__boxh = tuple(np.array(self.__boxl)*0.5) 
            self.__gridspacing = (self.__coords[1,0,0][0], self.__coords[0,1,0][1], self.__coords[0,0,1][2])
            self.__hvoxel = np.array([coords[1,0,0],coords[0,1,0],coords[0,0,1]])
        self.__hcell = self.__hvoxel * self.__Nx
        self.__volvoxel = np.linalg.det(self.__hvoxel)
        assert (np.abs(self.__volvoxel - np.linalg.det(self.__hcell) / self.__M) < 1e-5), "Volume of voxel != (Volume of cell / n voxels). This should be true!"

        self.__boxl = tuple(np.sqrt(np.sum(np.square(self.__hcell),axis=1)))
        self.__boxh = tuple(np.array(self.__boxl)*0.5) 


        # check if orthorhombic
        self.__orthorhombic = True
        hnorm = self.__hcell /  np.linalg.norm(self.__hcell, axis=0)
        if self.__ndim == 2 and np.dot(hnorm[0],hnorm[1]) != 0:
            self.__orthorhombic = False
        elif self.__ndim == 3 :
            if np.dot(hnorm[0],[1,0,0]) != 0 or np.dot(hnorm[1],[0,1,0]) != 0 or np.dot(hnorm[2],[0,0,1]) != 0:
                self.__orthorhombic = False
                print("Warning! Cell is not orthorhombic. This code was written for orthorhombic cells and non-orthorhombic support is in progress. So be careful, and check that the code is doing what you think it should!")


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
  
    def getNdim(self):
        return self.__ndim
    def getBoxl(self):
        return self.__boxl
    def getVolVoxel(self):
        return self.__volvoxel

    def getDomainStats(self, useMesh=True, plotMesh=False,outputMesh=False,add_periodic_domains=False, applyPBC=True):
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

            add periodic domains = true adds a center for mass at each of the locations for each periodic domain

        '''
        if useMesh and not self.__orthorhombic: 
            print("Warning: computing volume/area using mesh, but cell is not orthorhombic. This will lead to errors in the surface areas calculation of the domains")


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
            com[idomain,:] = self.calcDomainCOM(idomain+1,units='coord')
            
            if useMesh:

                if self.__ndim == 2:
                    # mesh domain
                    contours = self.meshSingleDomain(idomain+1,wrap_before_mesh=applyPBC)
                    assert (len(contours) == 1), "The contour should only be one curve, if not the area and volume calculations will be completely wrong!"

                    # get surface area (perimeter) and volume (area)
                    surface_area[idomain] = self.contour_perimeter(contours[0])
                    volume[idomain] = self.contour_area(contours[0])

                    if plotMesh: 
                        self.plotContours2D(contours,filename="mesh.{}.png".format(idomain+1))

                if self.__ndim == 3: 
                    # mesh domain
                    verts, faces, normals, values = self.meshSingleDomain(idomain+1,wrap_before_mesh=applyPBC)

                    # get surface area, volume and isoperimetric quotient
                    surface_area[idomain] = measure.mesh_surface_area(verts, faces)
                    volume[idomain] = self.mesh_volume(verts,faces)
                    if plotMesh: 
                        self.plotMesh3D(verts,faces, filename="mesh.{}.png".format(idomain+1))
                    if outputMesh:
                        self.writeMesh(verts,faces,fileprefix="mesh.{}.".format(idomain+1))

                IQ[idomain] = self.calcIQ(surface_area[idomain], volume[idomain])

            else:
                surface_area[idomain] = -1.0 #FIXME surface_area is currently not calculated if no mesh
                volume[idomain] = self.voxel_volume(idomain+1) # get volume from voxels
                IQ[idomain] = 0.0
        if add_periodic_domains:
            for idomain in range(1,self.__ndomains+1):  
                extracom = self.pbc_domain_locs(idomain,com[idomain-1])
                if extracom:
                    com = np.concatenate((com,extracom))
                    extra_num = len(extracom)
                    IQ = np.concatenate((IQ,[IQ[idomain-1]]*extra_num)) 
                    surface_area = np.concatenate((surface_area,[surface_area[idomain-1]]*extra_num)) 
                    volume = np.concatenate((volume,[volume[idomain-1]]*extra_num)) 
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

            # need to scale contours to be in terms of box dimensions
            for c in contours:
                c /= self.__Nx
                c *= self.__boxl

            if datafile:
                self.writeContours(contours,datafile)
            if plotfile:
                self.plotContours2D(contours,surface=mydensity,filename=plotfile)

            return contours

        elif self.__ndim == 3:
            mydensity = self.__fields[:,:,:, self.__density_field_index]
            verts, faces, normals, values = measure.marching_cubes_lewiner(mydensity, self.__density_threshold, spacing = self.__gridspacing)
            if datafile:
                raise NotImplementedError("Support for writing 3D mesh not implemented")
            if plotfile:
                self.plotMesh3D(verts,faces,filename=plotfile)
            return verts,faces, normals, values


    def meshSingleDomain(self,idomain, wrap_before_mesh=True):
        '''
        Function to:
        1) apply PBC to the domains so that an entire domain is continuous (ie not split across boundaries)
        2) Grab a little margin around each domain (the domain's "border") so that marching cubes can interpolate. The border is computed in identifyAndIndexDomains().
        3) Mesh the domain using marching cubes
        '''
        
        isdomain = (self.__regionID == idomain)
        #isborder = (self.__borderID == idomain)
        isborder = np.zeros(self.__Nx,dtype=np.bool)
        # convert to tuple to correctly set indicies of isborder
        isborder[tuple(self.__regionBorder[idomain-1])] = True

        if self.__ndim == 2:
            alldensity = self.__fields[:,:, self.__density_field_index]
        elif self.__ndim == 3:
            alldensity = self.__fields[:,:,:, self.__density_field_index]

        # center box and properties around center of mass (so that domains don't cross pbc)
        # np.roll is the key function here
        # if domains percolate then this will break 
        com_box = self.calcDomainCOM(idomain,units='box')
        com_coord = self.calcDomainCOM(idomain,units='coord')
        #coords_tmp = np.copy(self.__coords)
        for i in range(self.__ndim):
            shift = int(0.5*self.__Nx[i] - com_box[i])
            isdomain = np.roll(isdomain,shift,axis=i)
            isborder = np.roll(isborder,shift,axis=i)
            #coords_tmp = np.roll(coords_tmp,shift,axis=i)
            alldensity = np.roll(alldensity,shift,axis=i)

               
        # isolate the domain of interest
        isdomain_or_isborder = isdomain + isborder # since both bool, sum is the union of the two fields
        mydensity = np.zeros(self.__Nx)
        mydensity[isdomain_or_isborder] = alldensity[isdomain_or_isborder]
        
        #
        #tmp =mydensity[:,:,:,np.newaxis] 
        #viz.writeVTK('test.vtk',self.__coords,tmp)

        # plot for debugging
      #  import sys
      #  sys.path.append('/home/trent/Documents/college/polymers/ResearchTools/plot/')
      #  sys.path.append('../../')
      #  import PolyFTS_to_VTK
      #  AllCoords = np.reshape(coords_tmp,(self.__M, self.__ndim))
      #  AllCoords = AllCoords.T
      #  tmp = np.ravel(isdomain)
      #  tmp = np.resize(tmp,(1,len(tmp)))
      #  PolyFTS_to_VTK.writeVTK("isdomain.vtk", self.__Nx, True, self.__M, AllCoords,tmp)
      #  tmp = np.ravel(isborder)
      #  tmp = np.resize(tmp,(1,len(tmp)))
      #  PolyFTS_to_VTK.writeVTK("isborder.vtk", self.__Nx, True, self.__M, AllCoords,tmp)
      #  tmp = np.ravel(mydensity)
      #  tmp = np.resize(tmp,(1,len(tmp)))
      #  PolyFTS_to_VTK.writeVTK("mydensity.vtk", self.__Nx, True, self.__M, AllCoords,tmp)

        # mesh! (using scikit-image)
        if self.__ndim == 2:
            #raise NotImplementedError("Meshing in 2 dimensions is in development")
            contours = measure.find_contours(mydensity, self.__density_threshold) 
            #fixme! need to re-add center of mass to contours
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

    def writeMesh(self,verts,faces,fileprefix="mesh."):
       '''save mesh to a file'''
       np.savetxt(fileprefix + "verts.dat",verts,header='Autogenerated mesh file. Contains x y z positions of each vertex' )
       np.savetxt(fileprefix + "faces.dat",faces, header='Autogenerated mesh file. Contains vertex indicies of each triangle in mesh')
        


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

        self.__regionID = np.full(self.__Nx,-1,dtype=np.int32);
        # image_flags are only for the domains themselves, the image flags of the border are not needed
        self.__image_flags = np.zeros(list(self.__Nx) + [self.__ndim])
        ###self.__borderID = np.full(self.__Nx,0,dtype=np.int32);
        
        self.__regionBorder = [[]]

        region_number = 1;

        #this is where the recursive magic happens
        for i in np.ndindex(self.__Nx):
          if (self.__regionID[i] == -1):
            if (isdomain_array[i]):
              current_image_flag = np.zeros(self.__ndim)
              self.spread_region(i, region_number, isdomain_array,current_image_flag);
              self.__regionBorder.append([])
              region_number += 1;
            else:
              # note - dont assign borders here, this is acomplished inside of spread_region()
              self.__regionID[i] = 0;
              self.__image_flags[i]= np.zeros(self.__ndim)
        
        # now cleaning up
        nregions = region_number-1;
        
        # remove last element from lists (should be empty)
        assert (self.__regionBorder[-1] == [])
        del self.__regionBorder[-1] 
        
        # check that lengths of region structs are correct
        assert (len(self.__regionBorder) == nregions)

        # convert border and imageflag lists to numpy arrays
        for i in range(nregions):
            self.__regionBorder[i] = np.array(self.__regionBorder[i]).transpose()
       
        # change caching flag
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
            image_flag = tuple(neigh_image_flags[i])
            if (self.__regionID[neighbor] == -1):
                if (isdomain_array[neighbor]):
                  self.spread_region(neighbor, region_number, isdomain_array, image_flag);
                else:
                  self.__regionID[neighbor] = 0;
                  

            if self.__regionID[neighbor] == 0:
              
              # only append to list if neighbor isn't in there already
              if neighbor not in self.__regionBorder[region_number-1]:
                  # must have neighbors that are domain (since spread region is only called 
                  #   if coord_center is a domain). Therefore, it's a border
                  self.__regionBorder[region_number-1].append(neighbor)

                  # set image flags of non-domain adjacent to domain according to the domain
                  # basically, I need the border to have the correct image flags
                  # NOTE: image flags of borders aren't used anymore
                  #self.__regionBorderImageFlags[region_number-1].append(image_flag)

       
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

    def pbc_domain_locs(self,idomain,local_com):
        '''This function returns the locations of the other domains on the periodic boundary.  
        for example for a domain with its center on the corner of the box, it would return all
        the other box corners'''
        extra_com = []
        domain = (self.__regionID == idomain)
        local_flags = self.__image_flags[domain]
        unique_flags = set([])
        for i in range(np.shape(local_flags)[0]):
            unique_flags.add(tuple(local_flags[i]))
        unique_flags.remove((0,0,0))#remove duplicate com
        for flag in unique_flags:
            flag = np.array(flag)
            new_com = -1*flag*self.__boxl+local_com
            #find the location of the extra periodic com by adding the box length times the flag to the current com 
            extra_com.append(new_com)
            num_extra = len(extra_com)
        return extra_com
            

class DomainTracker:
    def __init__(self, boxl, vol_threshold=0.2):
        self.__boxl = boxl # stores max box position, (lower corner is at 0,0,0)
        self.__boxh = 0.5*boxl
        self.__ndim = len(boxl) 
        self.__vol_threshold = vol_threshold # volume threshold below which to ignore domains, percentage
        self.__is_init_pos = False
        #self.__msd = # stores average squared displacement (averaged over all micelles)

    def setInitialPositions(self,ndomains, com):
        ''' Set initial positions of domains 
        '''
        self.__ndomains = ndomains
        self.__pos0 = np.copy(com)  # initial position of each domain
        self.__pos_prev = np.copy(com)
        self.__imageflags = np.zeros((self.__ndomains,self.__ndim)) # which PBC image is the domain in (so that MSD can exceed the size of box)
        self.__sqdisp = np.zeros(self.__ndomains) # stores squared displacement of each micelle

        self.__is_init_pos = True

    def getMSD(self):
        ''' Returns mean squared displacement (averaged over all micelles)
        '''
        assert(self.__is_init_pos)
        return np.average(self.__sqdisp)

    def updateDisp(self,pos_curr, vol_curr):
        ''' Given current position of all micelles, update the squared displacement

            However, there's several things to be careful of:

            - The number of domains needs to be the same, the threshold will be needed to ignore domains with small volumes
                ^ frankly I'm not sure if this will work
            - the domain ordering may have changed so I need to figure out which domain_index in pos_curr corresponds to the domain_idx in pos_prev
        '''
        #check vol_curr
        pdb.set_trace()

        self.__pos_prev = np.copy(pos_curr)

        ndomains_curr = pos_curr.shape[0]
        assert(self.__ndomains == ndomains_curr)
       
        # compute index_map between pos_curr to pos_prev
        index_map = np.zeros(self.__ndomains, dtype=np.int32 )
        for i in range(self.__ndomains):
            min_dist, min_idx = 1e20, -1 
            for j in range(i+1,self.__ndomains):
                dist=getDistMinImage(pos_curr[i],self.__pos_prev[j])
                if dist < min_dist:
                    min_dist = dist
                    min_idx = j
            index_map[i] = min_idx
        #check index map
        pdb.set_trace()
        
        # now compute difference between pos_curr and pos_prev (without PBC)
        for i in range(self.__ndomains):
            dist = pos_curr[i] - pos_prev[index_map[i]] #FIXME check if this index_map usage is right
            for j in range(self.__ndim):
                # if domain has moved more than boxh then its crossed the PBC, so we update image flag accordingly
                if dist[j] > self.__boxh[j]:
                    self.__imageflags[i][j] += 1
                if dist[j] < -self.__boxh[j]:
                    self.__imageflags[i][j] -= 1

    def getDistMinImage(self,posA, posB):
        '''get the distance between posA and posB using the minimum image convention
        '''
        dist = posB - posA 
        for i in range(self.__ndim):
            if dist[i] >= self.__boxh[i]: 
                dist[i] -= self.__boxl[i]
            if dist[i] <= -self.__boxh[i]: 
                dist[i] += self.__boxl[i]
        return dist





