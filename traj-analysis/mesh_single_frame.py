#!/usr/bin/env python3
import sys
sys.path.append("/home/lequieu/Work/tools/lib/")

import iotools as io
import pdb
import viztools as viz
import numpy as np
from skimage import measure

if __name__ == "__main__":
    
    infile="density.bin"
    #ndim, Nx, orthorhombic, M, nfields, AllCoords, AllFields = io.ReadBinFile(infile)
    coords, fields = io.ReadBinFile(infile)

    mydensity = fields[:,:,:,0]
    density_threshold = 0.5
    gridspacing = (coords[1,0,0][0], coords[0,1,0][1], coords[0,0,1][2])
    boxl = tuple(coords[-1,-1,-1]) 

    verts, faces, normals, values = measure.marching_cubes_lewiner(mydensity, density_threshold, spacing = gridspacing)

    filename = "mesh_all.png"
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

    ax.set_xlim(0, boxl[0])  
    ax.set_ylim(0, boxl[1])  
    ax.set_zlim(0, boxl[2])  

    plt.tight_layout()
    if not filename:
        plt.show()
    else:
        plt.savefig(filename)

    plt.close()

