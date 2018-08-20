#!/usr/bin/env python

import numpy as np
import scipy.interpolate 
import re, glob
import matplotlib.pyplot as plt
import argparse
import itertools
import pdb


class PhaseBoundary:
    def __init__(self, phaseA, phaseB):
        self.phaseA = phaseA
        self.phaseB = phaseB
        self.boundary_points = []
        self.npoints = 0
        self.dist_threshold = [1e30, 1e30]

    def add_point(self,pos):
        self.boundary_points.append(pos)
        self.npoints += 1
    def get_pos_column(self,dim): # get entire column of pos as column, helpful for plotting
        out = [0]*self.npoints
        for i in range(self.npoints):
            out[i] = self.boundary_points[i][dim]
        return out
    def set_dist_threshold(self,dist_threshold):
        dim = len(self.boundary_points[0])
        if len(dist_threshold) != dim:
            raise ValueError("Cannot set dist_threshold, input dimension is %d but dimension of points are %d" % (len(dist_threshold), dim))
        self.dist_threshold = []
        for i in range(len(dist_threshold)):
            self.dist_threshold.append(float(dist_threshold[i]))
    def get_linesegments(self):
        ''' get a list of numpy 2d arrays where each element is a  line segment to be plotted'''
        dim = len(self.boundary_points[0])
        #dist_threshold = [0]*dim
        #dist_threshold[0] = 0.055 # along fA
        #dist_threshold[1] = 3     # along chiN

        # based on distances between adjacent points, 
        linesegments = [None]*1
        nline = 1

        pos0 = self.boundary_points[0]
        linesegments[nline-1] = np.array([pos0])
        for i in range(1,self.npoints):
            pos1 = self.boundary_points[i]
             
            mindist = 1e30
            jmin = -1
            append_flag = False  # if true append to previous line segment, if false create new line segment
            for j in range(nline):
                pos2 = linesegments[j][-1,:] # get last point on this line segment

                totaldist = 0
                local_append_flag = True
                for k in range(dim):
                    dist = abs(pos1[k]-pos2[k])
                    if (dist < self.dist_threshold[k]):
                        totaldist += dist*dist
                    else:
                        # if any distances are above threshold do not append locally
                        local_append_flag = False
                if local_append_flag and totaldist < mindist: 
                    mindist,jmin = totaldist, j # only compute min dist if within thresholds

                if local_append_flag: # if a single local_append_flag = True then append_flag = True
                    append_flag = True

            if append_flag: #append pos1 to closest linesegment
                linesegments[jmin] = np.append(linesegments[jmin],np.array([pos1]),axis=0)
            else: # create new linesegment
                linesegments.append(np.array([pos1]))
                nline += 1
        return linesegments

class PhaseBoundaryHolder:
    def __init__(self):
        self.boundaries = []
    def get_boundary(self,phaseA, phaseB):
        '''Return boundary number if it exists, otherwise return -1'''
        n = len(self.boundaries)
        for i in range(n):
            if (self.boundaries[i].phaseA == phaseA) and (self.boundaries[i].phaseB == phaseB):
                return i
            elif (self.boundaries[i].phaseA == phaseB) and (self.boundaries[i].phaseB == phaseA):
                return i
        return -1
    def set_dist_thresholds(self,dist_threshold):
        for boundary in self.boundaries:
            boundary.set_dist_threshold(dist_threshold)
        
    def new_boundary(self,phaseA, phaseB, pos):
        self.boundaries.append(PhaseBoundary(phaseA,phaseB))
        self.boundaries[-1].add_point(pos)

    def add_point(self, phaseA, phaseB, pos):
        iboundary = self.get_boundary(phaseA,phaseB)
        if iboundary == -1:
            self.new_boundary(phaseA,phaseB,pos)
        else:
            self.boundaries[iboundary].add_point(pos)

    def write(self,filename):
        f = open (filename, 'w')
        n = len(self.boundaries)
        for i in range(n):
            m = self.boundaries[i].npoints
            for j in range(m):
                f.write("%f %f # %s-%s\n" % (self.boundaries[i].boundary_points[j][0],self.boundaries[i].boundary_points[j][1], self.boundaries[i].phaseA, self.boundaries[i].phaseB))
        f.close()    

    def ax_plot_nodes(self, ax, nodes):
        '''utility to add nodes to the provided axis'''
        for node in nodes:
            x,y=node.pos[0], node.pos[1]
            ax.plot(x,y,marker='.',color='black',markersize=0.5)

    def ax_plot_boundaries(self,ax,marker=None, color=None, showlabels=False,ignorephases=[]):

        n = len(self.boundaries)
        for i in range(n):
            mymarker = marker.next()
            mycolor = color.next()
            if self.boundaries[i].phaseA in ignorephases or  self.boundaries[i].phaseB in ignorephases:
                continue

            lines = self.boundaries[i].get_linesegments()
            for j in range(len(lines)):
                if showlabels and j==0:
                    mylabel= "%s-%s" % (self.boundaries[i].phaseA, self.boundaries[i].phaseB)
                else:      
                    mylabel = ""
                x = lines[j][:,0]
                y = lines[j][:,1]
                #ax.plot(x,y,color=mycolor,marker=mymarker,linewidth=3,markersize=4)
                ax.plot(x,y,linewidth=2,color=mycolor,marker=mymarker,label=mylabel)

    def write(self,filename):
        
        f = open(filename,'wb')
        f.write("# raw data for phase boundaries. Formatted for gnuplot, though pretty self explanatory\n")

        n = len(self.boundaries)
        for i in range(n):
            lines = self.boundaries[i].get_linesegments()
            for j in range(len(lines)):
                if j==0:
                    f.write(b"#%s-%s\n" % (self.boundaries[i].phaseA, self.boundaries[i].phaseB))
                else:      
                    mylabel = ""
                x = lines[j][:,0]
                y = lines[j][:,1]
                data =np.vstack((x,y)).T
                np.savetxt(f,data,fmt="%0.8f")
            f.write(b"\n\n")
                #ax.plot(x,y,color=mycolor,marker=mymarker,linewidth=3,markersize=4)

        f.close()


    def plot(self,filename,plottype,nodes=None,xlabel='',ylabel='',axisrange=[None,None,None,None]):

        plt.rcParams['axes.labelsize'] = 16
        plt.rcParams['axes.labelsize'] = 16
        plt.rcParams['xtick.labelsize'] = 14
        plt.rcParams['ytick.labelsize'] = 14

        if plottype == 'plain':
            marker = itertools.cycle(('o')) 
            color  = itertools.cycle(('k')) 
            ignorephases=['sDIS']

        elif plottype == 'nodes':
            marker = itertools.cycle(('+', 'o', '*','v','^','<','>','s','p','h','H','x')) 
            color  = itertools.cycle(('tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan')) 
            #self.plot_with_nodes(filename,nodes,xlabel=args.xlabel, ylabel=args.ylabel)
            ignorephases=[]
        else:
            raise ValueError("Invalid plottype %s" % args.plottype)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if plottype == 'plain':
            self.ax_plot_boundaries(ax,marker=marker,color=color, showlabels=False,ignorephases=ignorephases)
        elif plottype == 'nodes':
            self.ax_plot_nodes(ax,nodes)
            self.ax_plot_boundaries(ax,marker=marker,color=color, showlabels=True)

            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        #ax.axis([None,None, 10, None])
            
        ax.axis(axisrange)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        #ax.autoscale()
        #ax.set_xlim(left=0.06,right=0.9)

        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        if filename == '':
            plt.show()
        else:
            if filename.split('.')[-1] == 'eps':
                plt.savefig(filename, bbox_inches='tight', format='eps')
            else:
                plt.savefig(filename, bbox_inches='tight')

#    def plot_plain(self,filename,xlabel='',ylabel=''):
#        
#        #fig = plt.figure(dpi=300)
#        fig = plt.figure()
#        ax = fig.add_subplot(111)
#        self.ax_plot_boundaries_plain(ax)
#
#                #ax.axis([None,None, 10, None])
#        ax.set_xlabel(xlabel)
#        ax.set_ylabel(ylabel)
#        #ax.autoscale()
#        ax.set_xlim(left=0.06,right=0.9)
#
#        #box = ax.get_position()
#        #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#        #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#
#    def plot_with_nodes(self,filename,nodes,xlabel='', ylabel=''):
#
#        fig = plt.figure()
#        ax = fig.add_subplot(111)
#
#        self.ax_plot_nodes(ax,nodes)
#        self.ax_plot_boundaries_color(ax)
#       
#        #ax.axis([None,None, 10, None])
#        ax.set_xlabel(xlabel)
#        ax.set_ylabel(ylabel)
#        box = ax.get_position()
#        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    def plot_heatmaps(self,nodes):
        fig =plt.figure(figsize=(5,8))
        #TODO clear figure
        npts = 100 #resolution of interpolation
        nnodes = len(nodes) 
        phases = nodes[0].phases
        nphases = len(phases)
        nrow, ncol = 4,2


        filename="fig_heatmaps.png" 
        for j in range(nphases):
            targetphase=phases[j]
            ax = fig.add_subplot(nrow, ncol,j+1)

            ax.set_title(targetphase)

            x = np.zeros((nnodes,))
            y = np.zeros((nnodes,))
            z = np.zeros((nnodes,))
            
            for i in range(nnodes):
                Ftarget = nodes[i].get_phase_F(targetphase)
                Fmin = nodes[i].get_min_F()
                x[i] = nodes[i].pos[0]
                y[i] = nodes[i].pos[1]
                z[i] = Ftarget - Fmin
            xi = np.linspace(min(x),max(x),npts)
            yi = np.linspace(min(y),max(y),npts)
            zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
            z_min, z_max = zi.min(), np.abs(zi).max()
            p = ax.pcolormesh(xi,yi,zi,vmin=0,vmax=1e-3)

            n = len(self.boundaries)
            for i in range(n):
                mylabel="%s-%s" % (self.boundaries[i].phaseA, self.boundaries[i].phaseB)
                x = self.boundaries[i].get_pos_column(0)
                y = self.boundaries[i].get_pos_column(1)
                ax.plot(x,y,linewidth=0,marker='o',markersize=1)
            ax.axis([None,None, 10, None])


        #fig.colorbar(p)
        plt.subplots_adjust(wspace=0.2, hspace=0.5)
        plt.savefig(filename)


        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        #plt.savefig(filename)
     
        

def linesegment(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x,y
    else:
        #raise ValueError('Lines L1 and L2 do not intersect')
        return L1[0],L1[1]

class Node:
    __slots__ = 'pos', 'dim', 'phases', 'F', 'nextnode', 'visited', 'F_THRESHOLD'
    def __init__(self,pos, phases,F):
        self.pos = pos
        dim = len(self.pos)
        if len(F) != len(phases):
            raise ValueError('length of phases and F are not equal!',phases, F)
        self.phases = phases
        self.F = F
        self.nextnode = [None]*2*dim
        self.visited = False
        self.F_THRESHOLD = 1e-4
        self.check_all_phases_DIS()
    def get_phase_index(self,myphase):
        for i,phase in enumerate(self.phases):
            if myphase in phase:
                return i
        raise ValueError('%s phase not found in node at %s' % (myphase,self.pos))

    def check_all_phases_DIS(self):
        DIS_idx = self.get_phase_index('DIS')
        flag=True
        for i in range(len(self.F)):
            if abs(self.F[i] - self.F[DIS_idx]) > self.F_THRESHOLD:
                flag = False
        self.is_all_DIS = flag; 
		#self.phases[DIS_idx] = 'sDIS' # change name to spinodal DIS (sDIS)
        
    def get_phase_F(self,phase):
        return self.F[self.get_phase_index(phase)]
    def get_min_F(self):
        idx = self.get_min_idx()
        if idx > len(self.phases) or idx < 0:
            raise ValueError('Invalid minimum index %d > %d' % (idx,len(self.phases)))
        return self.F[idx]

    def get_min_phase(self):
		idx = self.get_min_idx()
		return self.phases[idx]
		#if not self.is_all_DIS:
		#	idx = self.get_min_idx()
		#	return self.phases[idx]
		#else:
		#	return "sDIS" # if all are dis, then return a different DIS phase

    def get_min_idx(self):
        val, idx = min((val, idx) for (idx, val) in enumerate(self.F)) 

        # custom check for DIS boundaries
        DIS_idx = self.get_phase_index('DIS')
        DIS_val = self.F[DIS_idx]
        if ( abs (val - DIS_val) < self.F_THRESHOLD ):
            idx = DIS_idx
        
        return idx

    def set_node_up(self,node,dim=None): self.nextnode[2*dim+1] = node
    def set_node_dn(self,node,dim=None): self.nextnode[2*dim+0] = node
    def get_node_up(self,dim=None): return self.nextnode[2*dim+1]
    def get_node_dn(self,dim=None): return self.nextnode[2*dim+0]
    def get_dim(self) : return len(self.pos)

class Iterator:
    def __init__(self,mynode):
        self.mynode = mynode 
    def move_up(self,dim=None): self.mynode = self.mynode.get_node_up(dim=dim)
    def move_dn(self,dim=None): self.mynode = self.mynode.get_node_dn(dim=dim)
    def peek_up(self,dim=None): return self.mynode.get_node_up(dim=dim)
    def peek_dn(self,dim=None): return self.mynode.get_node_dn(dim=dim)


def compare_node_pos(nodeA, nodeB,dim_to_ignore=-1):
    '''compare the position of two nodes, and return true if equal, false otherwise. Can ignore a dimension if necessary'''
    if nodeA.get_dim() != nodeB.get_dim():
        raise ValueError('Dimension of nodeA and nodeB are not equal!')
    dim = nodeA.get_dim()
    match = True
    for i in range(dim):
        if ((dim_to_ignore != -1) and (i == dim_to_ignore)): continue;
        elif nodeA.pos[i] != nodeB.pos[i]: match = False
    return match 


def interpolate_nodes(nodeA, nodeB, dim):

    # linear interpolate to get phase boundary
    min_phaseA = nodeA.get_min_phase()
    min_phaseB = nodeB.get_min_phase()

    # whats going on here can be confusing so heres a picture
    #
    #        ^
    #        |
    #        |    minphaseB     minphaseA
    #        |    (point 3)    (point 2)
    #        |              \  /                               
    #      F |               \/                                
    #        |               /\                                
    #        |              /  \                               
    #        |     (point 1)    (point 4)
    #        |     minphaseA    minphaseB
    #        |
    #        |      NodeA         Node B
    #        ----------------------------------->
    #              fA,chiN or something else

    p1 = (nodeA.pos[dim], nodeA.get_phase_F(min_phaseA))
    p2 = (nodeB.pos[dim], nodeB.get_phase_F(min_phaseA))
    p3 = (nodeA.pos[dim], nodeA.get_phase_F(min_phaseB))
    p4 = (nodeB.pos[dim], nodeB.get_phase_F(min_phaseB))
    
    lineA = linesegment(p1,p2)
    lineB = linesegment(p3,p4)
    pos_boundary_dim,F_boundary = intersection(lineA,lineB)
    
    pos_boundary = list(nodeA.pos)
    pos_boundary[dim] = pos_boundary_dim

    return tuple(pos_boundary)

def plot_node_connectivity(nodes):
    #plot it
    x_head_width=0.3
    x_head_length=0.01
    y_head_width=0.01
    y_head_length=0.5

    plt.axis([0,1,10,45])
    scale=0.7
    for node in nodes:
        x = node.pos[0]
        y = node.pos[1]
        plt.plot(x,y,'.')
        for i in [0,1]:
            nextnode = node.get_node_up(dim=i)
            if nextnode != None:
                xnext,ynext = nextnode.pos[0], nextnode.pos[1]
                dx,dy = xnext - x, ynext - y
                plt.plot(xnext,ynext,'.')
                if i == 0:
                    plt.arrow(x,y,scale*dx,scale*dy, head_width=x_head_width, head_length=x_head_length, fc='k', ec='k')
                elif i == 1:
                    plt.arrow(x,y,scale*dx,scale*dy, head_width=y_head_width, head_length=y_head_length, fc='k', ec='k')

            nextnode = node.get_node_dn(dim=i)
            if nextnode != None:
                xnext = nextnode.pos[0]
                ynext = nextnode.pos[1]
                dx = xnext - x
                dy = ynext - y
                plt.plot(xnext,ynext,'.')
                if i == 0:
                    plt.arrow(x,y,scale*dx,scale*dy, head_width=x_head_width, head_length=x_head_length, fc='k', ec='k')
                elif i == 1:
                    plt.arrow(x,y,scale*dx,scale*dy, head_width=y_head_width, head_length=y_head_length, fc='k', ec='k')
    plt.savefig("fig_node_connectivity.png")
    #plt.show()


#this assumes that F0_phases.dat has been generated in each chiN*/phiA* directory
def initialize_nodes(dirs,fnmeIn):
    ''' 
        Go through directories (dirs) and look for filename (fnmeIn) and load the phase information 
        the file is expected to contain:

        <Phase name> <Free energy of that phase>
        HEXPhase 3.2
        LAMPhase 2.4
        DISPhase 1.2
    '''
    dirs.sort(reverse=True)
    nphases0 = None;
    nodes = []
    for mydir in dirs:
        
        dir1=mydir.split('/')[0]
        dir2=mydir.split('/')[1]

        chiN = float(re.sub("[a-z,A-Z]","",dir1.split('_')[-1]))
        phiA = float(re.sub("[a-z,A-Z]","",dir2))

        fnme = "%s/%s" % (mydir,fnmeIn)
        
        # read phases
        phases = []
        F=[]
        with open(fnme,'r') as f:
            for line in f:
                l=line.split()
                phase = re.sub("Phase","",l[0]) #remove Phase
                phases.append(phase)
                F.append(float(l[1]))

        pos = (phiA,chiN)
        nodes.append(Node(pos,phases,F))

        f.close()

    populate_nodes_neighbors(nodes)
    return nodes
   
def populate_nodes_neighbors(nodes):
    # now populate each node's neighbors, this is N^2 complexity, but shouldnt be noticable 
    for inode in range(len(nodes)):
        dim = nodes[inode].get_dim()
        for idim in range(dim):
            xinode = nodes[inode].pos[idim]

            min_dist_up, min_dist_dn = 1e30,1e30
            node_dn,node_up = None,None
            for jnode in range(len(nodes)):
                
                # skip if nodes are the same, or if they're not at the same position along all dim other than idim
                if inode == jnode: continue
                elif (not compare_node_pos(nodes[inode],nodes[jnode],dim_to_ignore=idim)): continue

                xjnode = nodes[jnode].pos[idim]
                dist = abs(xjnode - xinode)
                if (xjnode < xinode) and (dist < min_dist_dn): min_dist_dn = dist; node_dn = jnode
                if (xjnode > xinode) and (dist < min_dist_up): min_dist_up = dist; node_up = jnode

            #if nodes[inode].pos == (0.3,19.0):
            #    pdb.set_trace()        
            if node_dn != None:
                nodes[inode].set_node_dn(nodes[node_dn],dim=idim)
            if node_up != None:
                nodes[inode].set_node_up(nodes[node_up],dim=idim)


def calc_phase_boundaries(nodes):
    # now we have a data structure where each phase point is stored in a node, and each node contains pointers to its nearest neighbors
    # populate the boundaryholder object
    boundaryholder = PhaseBoundaryHolder();

    nvisited = 0
    nnodes = len(nodes)
    #it = Iterator(nodes[0])
    #while (nvisited < nnodes):
    for inode in range(nnodes):
        my_node = nodes[inode]
        #my_node = it.mynode
        my_phase = my_node.get_min_phase()
        #if my_node.pos == (0.15,25.0):
        #    pdb.set_trace()

        for idim in range(my_node.get_dim()):
            # dont do vertical interpolation if DIS is involved (this cleans up the phase diagram)
            #if idim == 1 and "DIS" in my_phase: continue;

            #up_node  = it.peek_up(dim=idim)        
            up_node  = my_node.get_node_up(dim=idim)        
            #if up_node != None and up_node.visited == False and not up_node.is_all_DIS:
            if up_node != None and up_node.visited == False:
                up_phase = up_node.get_min_phase()
                #if idim == 1 and "DIS" in up_phase: continue;

                if (my_phase != up_phase): 
                    pos_boundary = interpolate_nodes(my_node,up_node,idim)
                    #if ("DIS" in my_phase) or ("DIS" in up_phase):
                    #    my_phase = up_phase = "DIS"
                    if (my_node.is_all_DIS): my_phase = 'sDIS' # change name to spinodal DIS (sDIS)
                    if (up_node.is_all_DIS): up_phase = 'sDIS'
                    if 'DIS' in my_phase and 'DIS' in up_phase:  continue
                    boundaryholder.add_point(my_phase,up_phase,pos_boundary)    

            #dn_node  = it.peek_dn(dim=idim)        
            dn_node  = my_node.get_node_dn(dim=idim)        
            #if dn_node != None and dn_node.visited == False and not dn_node.is_all_DIS:
            if dn_node != None and dn_node.visited == False:
                dn_phase = dn_node.get_min_phase()
                #if idim == 1 and "DIS" in dn_phase: continue;

                if (my_phase != dn_phase): 
                    pos_boundary = interpolate_nodes(my_node,dn_node,idim)
                    #if ("DIS" in my_phase) or ("DIS" in dn_phase):
                    #    my_phase = dn_phase = "DIS"
                    if (my_node.is_all_DIS): my_phase = 'sDIS' # change name to spinodal DIS (sDIS)
                    if (dn_node.is_all_DIS): dn_phase = 'sDIS'
                    if 'DIS' in my_phase and 'DIS' in dn_phase:  continue
                    boundaryholder.add_point(my_phase,dn_phase,pos_boundary)    

        my_node.visited = True     
        nvisited += 1;
    
    return boundaryholder

# ==============================================================================
#   Begin Main
# ==============================================================================

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Tool to compute phase boundaries')
    parser.add_argument('-f', '--filename', action='store', default='F0_phases.dat',help='file that contains the phases and free energies at each phase point')
    parser.add_argument('-d', '--dirs', action='store', nargs='+', default=glob.glob("tau*/phiA*"),help='list of directories that contain each phase point')
    parser.add_argument('-o', '--outfig', action='store', default='',help='name of output figure file')
    parser.add_argument('--raw', action='store', default='',help='name of raw output file (for plotting in another program')
    parser.add_argument('-t', '--plottype', action='store', default='nodes',help='type of plot to generate')
    parser.add_argument('--xlabel', action='store', default="$f_A$",help='label for xaxis, can use \'$\' to write latex')
    parser.add_argument('--ylabel', action='store', default="$\chi N$",help='')
    parser.add_argument('--axisrange', action='store', nargs=4, default=[None,None,None,None],help='')
    parser.add_argument('--linecutoff', action='store', nargs=2, default=[1e30,1e30],help='maximum length of lines to draw in phase diagrams, useful to clean them up')
    print "IMPLEMENT CUSTOM AXIS RANGES AND LABELS FROM COMMAND LINE"
    args = parser.parse_args()
    
    #fnmeIn="F0_phases.dat"
    #dirs=glob.glob("tau*/phiA*");
    if args.axisrange != [None, None, None, None]:
        args.axisrange = [float(x) for x in args.axisrange]
    print args.axisrange

    nodes = initialize_nodes(args.dirs, args.filename)
    boundaryholder = calc_phase_boundaries(nodes)
    dist_threshold = args.linecutoff
    boundaryholder.set_dist_thresholds(dist_threshold)

    #boundaryholder.write("boundaries.dat")
    #boundaryholder.plot_plain("fig_boundaries.png")
    #boundaryholder.plot_heatmaps(nodes)
    #plot_node_connectivity(nodes)
    if args.raw != '':
        print "Saving raw phase boundary data to \'%s\'" % args.raw
        boundaryholder.write(args.raw)

    print "Plotting with type \'%s\'" % args.plottype
    boundaryholder.plot(args.outfig,args.plottype, nodes=nodes,xlabel=args.xlabel, ylabel=args.ylabel,axisrange=args.axisrange)



