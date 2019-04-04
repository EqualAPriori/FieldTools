#!/usr/bin/env python3

import numpy as np
import scipy.interpolate 
import re, glob
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import argparse
import itertools
import pdb
import copy
import os,sys

class PhaseBoundary:
    def __init__(self, phaseA, phaseB):
        self.phaseA = phaseA
        self.phaseB = phaseB
        self.boundary_points = []
        self.boundary_F = []
        self.npoints = 0
        self.dist_threshold = [1e30, 1e30]
    def add_point(self,pos,F=None):
        self.boundary_points.append(pos)
        self.boundary_F.append(F)
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

        # now I need to sort the points within each line segment
        for il,l in enumerate(linesegments):
            if l.shape[0] > 2:
                # start with the top-left value (min x, max y) and then work the way down the curve
                maxy = np.max(l[:,1])
                ismaxy = l[:,1] == maxy
                if sum(ismaxy) != 1: # if more than one maxy, then use the one with the minx
                    ltmp = np.copy(l)
                    ltmp[~ismaxy] = 1e30 # set values that aren't maxy to a big value
                    minx = np.min(ltmp[:,0])
                    isminx = (ltmp[:,0] == minx) 

                    isboth = isminx * ismaxy
                    
                    index = isboth.argmax() # contains the row in l that is at the top-left
                    #print("Error! more than one point matches minx miny")
                    #pdb.set_trace()
                else:
                    index = ismaxy.argmax() # contains the row in l that is at the top-left
            else:
                index = 0
               
            # create two arrays
            # 1) keep track of whether the index has been used
            is_index_used = np.zeros(l.shape[0],dtype=np.bool)
            # 2) store the sorted line
            l_sorted = np.zeros(l.shape)
            
            # store 1st point
            l_sorted[0,:] = l[index,:]
            is_index_used[index] = True
            
            # compute a scale factor that weights the differences in x and y axis ranges
            sizeofdim = [0]*dim
            for k in range(dim):
                sizeofdim[k] = np.max(l[:,k]) - np.min(l[:,k])
            
            if l.shape[0] > 2:
                i=0 # stores index of l_sorted
                while i < l.shape[0]-1:
                    # find closest (unused) element of "l" and add it to "l_sorted"
                    jmin=0
                    mindist = 1e30
                    for j in range(l.shape[0]):
                        if is_index_used[j] != True:
                            dist = 0
                            for k in range(dim):
                                dist += (l_sorted[i][k]-l[j][k]) * (l_sorted[i][k]-l[j][k]) / sizeofdim[k] / sizeofdim[k]
                            dist = np.sqrt(dist) 
                            if dist < mindist:
                                jmin = j
                                mindist = dist
                    try:
                        l_sorted[i+1,:] = l[jmin,:]
                    except:
                        pdb.set_trace()
            
                    is_index_used[jmin] = True
                    i += 1
            else:
                l_sorted = np.copy(l)
                

            linesegments[il] = np.copy(l_sorted)
        
        #finally return
        return linesegments

class PhaseBoundaryHolder:
    def __init__(self):
        self.boundaries = []
        self.phase_F_dict = {}
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
        
    def new_boundary(self,phaseA, phaseB, pos,F=None):
        self.boundaries.append(PhaseBoundary(phaseA,phaseB))
        self.boundaries[-1].add_point(pos,F)

    def add_point(self, phaseA, phaseB, pos,F=None):
        iboundary = self.get_boundary(phaseA,phaseB)
        if iboundary == -1:
            self.new_boundary(phaseA,phaseB,pos,F)
        else:
            self.boundaries[iboundary].add_point(pos,F)


    def ax_plot_nodes(self, ax, nodes):
        '''utility to add nodes to the provided axis'''
        for node in nodes:
            x,y=node.pos[0], node.pos[1]
            ax.plot(x,y,marker='.',color='black',markersize=0.5)

    def ax_plot_boundaries(self,ax,marker=None, color=None, showlabels=False):

        n = len(self.boundaries)
        for i in range(n):
            mymarker = next(marker)
            mycolor = next(color)
            lines = self.boundaries[i].get_linesegments()
            for j in range(len(lines)):
                if showlabels and j==0:
                    mylabel= "%s-%s" % (self.boundaries[i].phaseA, self.boundaries[i].phaseB)
                else:      
                    mylabel = ""
                x = lines[j][:,0]
                y = lines[j][:,1]
                #ax.plot(x,y,color=mycolor,marker=mymarker,linewidth=3,markersize=4)

                ax.plot(x,y,color=mycolor,marker=mymarker,label=mylabel)

    def write(self,filename,dim=2,z=None):
        '''
        this method writes phase boundary locations and free energies to files for plotting in another program
        '''
        if dim ==2:
           with open(filename,'w') as f:
                f.write("# raw data for phase boundaries. Formatted for gnuplot, though pretty self explanatory\n")

                n = len(self.boundaries)
                for i in range(n):
                    lines = self.boundaries[i].get_linesegments()
                    for j in range(len(lines)):
                        if j==0:
                            f.write("#%s-%s\n" % (self.boundaries[i].phaseA, self.boundaries[i].phaseB))
                        else:      
                            mylabel = ""
                        x = lines[j][:,0]
                        y = lines[j][:,1]
                        data =np.vstack((x,y)).T
                        np.savetxt(f,data,fmt="%0.8f")
                    f.write("\n\n")
        elif dim == 1:
            with open(filename,'w+') as f:
                f.write("#Raw free energy data pairs for each microphase\n")
                for p in self.phase_F_dict:
                    f.write("#"+p+"\n")
                    data = np.vstack(self.phase_F_dict[p]).T
                    np.savetxt(f,data,fmt="%0.8f")
                    f.write("\n\n")
        elif dim ==3:
            with open(filename,'a') as f:
                n = len(self.boundaries)
                for i in range(n):
                    lines = self.boundaries[i].get_linesegments()
                    for j in range(len(lines)):
                        if j==0:
                            f.write("#%s-%s\n" % (self.boundaries[i].phaseA, self.boundaries[i].phaseB))
                        else:      
                            mylabel = ""
                        x = lines[j][:,0]
                        y = lines[j][:,1]
                        my_z = [z]*len(x)
                        data =np.vstack((x,y,my_z)).T
                        np.savetxt(f,data,fmt="%0.8f")
                        f.write("\n\n")
                f.write("\n\n")



                  

    def plot(self,filename,plottype,nodes=None,xlabel='',ylabel='',axisrange=[None,None,None,None],n=2, z=None,colors=None,symbols=None,xticks=None,yticks=None,aspect=None,refPhase=None):
        if plottype == 'plain':
            marker = itertools.cycle(('o')) 
            color  = itertools.cycle(('k')) 

        elif plottype == 'nodes':
            marker = itertools.cycle(('+', 'o', '*','v','^','<','>','s','p','h','H','x')) 
            color  = itertools.cycle(('tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan')) 
            #self.plot_with_nodes(filename,nodes,xlabel=args.xlabel, ylabel=args.ylabel)
        elif plottype == 'nodesblack':
            marker = itertools.cycle(('+', 'o', '*','v','^','<','>','s','p','h','H','x'))
            color  = itertools.cycle(('k'))
        elif plottype == 'simplecolors' or 'nolegend':
            marker = itertools.cycle(('+', 'o', '*','v','^','<','>','s','p','h','H','x'))
            color  = itertools.cycle(('r','b','g','k','m','c','y'))
        else:
            raise ValueError("Invalid plottype %s" % plottype)
        if colors:
           color = itertools.cycle(colors)
        if symbols:
           marker = itertools.cycle(symbols)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if plottype == 'plain':
            if n == 1:
                self.plot1d(ax,nodes,refPhase=refPhase)
                self.plotvert(ax)
            elif n == 2:
                self.ax_plot_boundaries(ax,marker=marker,color=color, showlabels=False)
        elif plottype == 'nodes' or plottype == 'nodesblack' or plottype == 'simplecolors' or plottype == 'nolegend':
            if n==1:
                self.plot1d(ax,nodes,color=color,marker=marker,showlabels=True,refPhase=refPhase)
                self.plotvert(ax)
            elif n==2:
                if plottype != 'nolegend':
                    self.ax_plot_nodes(ax,nodes)
                self.ax_plot_boundaries(ax,marker=marker,color=color, showlabels=True)

            if plottype != 'nolegend':
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


        ax.axis(axisrange)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        #ax.autoscale()
        #ax.set_xlim(left=0.06,right=0.9)

        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        if xticks:
           ax.set_xticks(xticks)
        if yticks:
           ax.set_yticks(yticks)
        if filename == '':
            plt.gcf().set_tight_layout(True) #makes the graph look better by removing whitespace
            plt.show()
        else:
            if aspect:
                 ratio = float(aspect)
                 xleft, xright = ax.get_xlim()
                 ybottom, ytop = ax.get_ylim()
                 # the abs method is used to make sure that all numbers are positive
                 # because x and y axis of an axes maybe inversed. 
                 # convert to figure coords from data coords
                 ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)
            if filename.split('.')[-1] == 'eps':
                plt.savefig(filename, bbox_inches='tight', format='eps')
            else:
                plt.savefig(filename, bbox_inches='tight')
    
    def plot1d(self,ax,nodes,marker='.', color='k', showlabels=False,refPhase=None):
        '''
        This method loops over each node in the boundary holder and 
        finds what phases are present, and then plots the free energies
        '''
        
        # assume all nodes have same dimension...should be reasonable...
        dim2plot = -1
        if nodes[0].dim == 1:
            dim2plot = 0
        elif nodes[0].dim == 2:
            # determine which dimension is changing, that should be the x axis
            first = True 
            for i,n in enumerate(nodes):
                pos = n.pos
                if i != 0:
                    isdimequal = list(np.array(pos) == np.array(posprev) )
                    if sum(isdimequal) != 1:
                       raise RuntimeError("More than one dimension is changing within nodes. Dont know which dimension to plot as x-axis in 1d plot")
                    if dim2plot == -1:
                        dim2plot = isdimequal.index(False) # the index that is not equal is the one we want to plot
                    elif dim2plot != isdimequal.index(False):
                        raise RuntimeError("dim2plot is changing throughout nodes. Can't determine which dim to plot as x-axis is 1d plot")
                posprev = pos


        phases = copy.copy(nodes[1].phases) #initialize possible phase list
        for n in nodes:
            for p in n.phases:
                if p not in phases:
                    phases.append(p)
        for p in phases:
            mymarker = next(marker)
            mycolor = next(color)
            x,y = [],[]
            for n in nodes:
                if n.has_phase(p):
                    x.append(n.pos[dim2plot])
                    F=n.get_phase_F(p)
                    if refPhase != None:
                        Fref = n.get_phase_F(refPhase)
                        F -= Fref
                    y.append(F)
                    #pdb.set_trace()
            if showlabels:
                ax.plot(x,y,marker=mymarker,color=mycolor,label = p)
            else:
                ax.plot(x,y,marker=mymarker,color=mycolor)
            self.phase_F_dict[p] = (x,y)

    def plotvert(self,ax):
        '''
        this method adds a vertical line for each free energy curve intersection on the 
        1 dimensional free energy curve plot
        it also adds the words for each phase
        '''
        bottom = ax.get_ylim()[0]#get the bottom area to be plotted
        left = ax.get_xlim()[0]#find the left side of the graph for centering text
        right = ax.get_xlim()[1]
        xlist = [left,right]#make a list of the x values that the text can be centered inbetween
        for i in range(len(self.boundaries)):
           x = self.boundaries[i].get_pos_column(0)#intersection locations
           top = self.boundaries[i].boundary_F
           for j in range(len(top)):
               y = np.linspace(bottom,top[j])
               xx = np.array([x[j]]*len(y))#make an array the size of the linearly spaced y
               ax.plot(xx,y,linestyle=':',color='k')
               xlist.append(x[j])    
        #print(xlist)
        xlist.sort()
        for i in range(1,len(xlist)):
            x = (xlist[i]+xlist[i-1])/2 #center point
            mindiff = 1e30
            center_node = None
            for n in nodes:
                if abs(n.pos[0]-x) < mindiff:
                    mindiff = abs(n.pos[0]-x)
                    centernode = n
            minphase = centernode.get_min_phase()
            y = (centernode.get_min_F()-bottom)/2.0+bottom #vertical center point
            t= ax.text(x,y,centernode.get_min_phase(),fontsize='large',color ='k',horizontalalignment='center')
   
    def plot3d(self,ax,filename,marker=None,color=None,nodes=None,axisrange=[None,None,None,None],z=None):
        '''
        this method plots the phase boundaries on a given dictionary of axes, but does not show the plot
        this allows multiple boundary holders to be plotted, to generate a layered plot
        '''
        for boundary in reversed(self.boundaries):
        #   if not z:
               lines = boundary.get_linesegments()
               for j in range(len(lines)):
                   showlabels=True
                   if showlabels:
                       mylabel= boundary.phaseA+'-'+boundary.phaseB#+' at '+args.zlabel+' '+str(z)
                   else:
                       mylabel = ""
                   x = lines[j][:,0]
                   y = lines[j][:,1]
                   #ax.plot(x,y,color=mycolor,marker=mymarker,linewidth=3,markersize=4)
                   my_ax =ax[mylabel] #plot on the axis with the appropriate phase boundary by accessing a dictionary
                   my_ax.plot(x,y,linewidth=3,label=mylabel,color= color,marker=marker,)
                   my_ax.text(x[-1]+.01,y[-1]-.01,args.zlabel + '=' + str(z),fontsize='small',fontweight='semibold',verticalalignment='center',horizontalalignment='left') #graph text at the toward the end of the line 

            #pdb.set_trace()
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
    __slots__ = 'pos', 'dim', 'phases', 'F', 'nextnode', 'visited', 'F_THRESHOLD','is_all_DIS'
    def __init__(self,pos, phases,F):
        self.pos = pos
        self.dim = len(pos)
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
        #    idx = self.get_min_idx()
        #    return self.phases[idx]
        #else:
        #    return "sDIS" # if all are dis, then return a different DIS phase

    def get_min_idx(self):
        val, idx = min((val, idx) for (idx, val) in enumerate(self.F)) 

        # custom check for DIS boundaries
        DIS_idx = self.get_phase_index('DIS')
        DIS_val = self.F[DIS_idx]
        if ( abs (val - DIS_val) < self.F_THRESHOLD ):
            idx = DIS_idx
        
        return idx
    def has_phase(self,phase):
        return phase in self.phases
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

    return tuple(pos_boundary),F_boundary

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

        <Phase name> <Free energy of that phase>   [exit status]
        HEXPhase 3.2 2 
        LAMPhase 2.4 2
        DISPhase 1.2 0
    '''
    dirs.sort(key=numerical_sort)
    nphases0 = None;
    nodes = []
    for mydir in dirs:
        
        dir1=mydir.split('/')[-2]
        dir2=mydir.split('/')[-1]

        chiN = float(re.sub("[a-z,A-Z]","",dir1.split('_')[-1]))
        phiA = float(re.sub("[a-z,A-Z]","",dir2))
#        split = re.split('/',mydir)
#        pos = []
#        for word in split:
#             pos.append(float(re.sub('[^0-9.]','',re.split('_',word)[-1]))) #add the last number from each directory level to the list of position
        fnme = "%s/%s" % (mydir,fnmeIn)
        #pos = tuple(reversed(pos)) 
        # read phases
        phases = []
        F=[]
        status=[]
        try:
            with open(fnme,'r') as f:
                for line in f:
                    l=line.split()
                    phase = re.sub("Phase","",l[0]) #remove Phase
                    phases.append(phase)
                    F.append(float(l[1]))
            pos = (phiA,chiN)
            nodes.append(Node(pos,phases,F))
        except FileNotFoundError:
           print('WARNING: Free energy data not found at {} try using extractF0 to generate it'.format(fnme))

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


def calc_phase_boundaries(nodes,n=2):
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
            #Skip dimension if it is not one of the dimensions to interpolate on
            if idim not in args.interp_dimension:
                continue
            # dont do vertical interpolation if DIS is involved (this cleans up the phase diagram)
            #if idim == 1 and "DIS" in my_phase: continue;

            #up_node  = it.peek_up(dim=idim)        
            up_node  = my_node.get_node_up(dim=idim)        
            #if up_node != None and up_node.visited == False and not up_node.is_all_DIS:
            if up_node != None and up_node.visited == False:
                up_phase = up_node.get_min_phase()
                #if idim == 1 and "DIS" in up_phase: continue;

                if (my_phase != up_phase): 
                    try:
                       pos_boundary,F_boundary = interpolate_nodes(my_node,up_node,idim)
                       #if ("DIS" in my_phase) or ("DIS" in up_phase):
                       #    my_phase = up_phase = "DIS"
                       if (my_node.is_all_DIS): my_phase = 'sDIS' # change name to spinodal DIS (sDIS)
                       if (up_node.is_all_DIS): up_phase = 'sDIS'
                       if 'DIS' in my_phase and 'DIS' in up_phase:  continue
                       boundaryholder.add_point(my_phase,up_phase,pos_boundary,F_boundary)    
                    except ValueError as e:
                       print('WARNING: ' +str(e)+" skipping this possible phase boundary")

            #dn_node  = it.peek_dn(dim=idim)        
            dn_node  = my_node.get_node_dn(dim=idim)        
            #if dn_node != None and dn_node.visited == False and not dn_node.is_all_DIS:
            if dn_node != None and dn_node.visited == False:
                dn_phase = dn_node.get_min_phase()
                #if idim == 1 and "DIS" in dn_phase: continue;

                if (my_phase != dn_phase): 
                    try:
                        pos_boundary,F_boundary = interpolate_nodes(my_node,dn_node,idim)
                        #if ("DIS" in my_phase) or ("DIS" in dn_phase):
                        #    my_phase = dn_phase = "DIS"
                        if (my_node.is_all_DIS): my_phase = 'sDIS' # change name to spinodal DIS (sDIS)
                        if (dn_node.is_all_DIS): dn_phase = 'sDIS'
                        if 'DIS' in my_phase and 'DIS' in dn_phase:  continue
                        boundaryholder.add_point(my_phase,dn_phase,pos_boundary,F_boundary)    
                    except ValueError as e:
                         print('WARNING '+str(e)+" skipping this possible phase boundary")

        my_node.visited = True     
        nvisited += 1;
    
    return boundaryholder



def numerical_sort(mydir):  
    '''
    This function weights a directory, so that it can be sorted
    based on the numbers found in each portion of the directory
    a greater weight is applied to numbers on the left
    '''
    weight = 0
    words = re.split('/',mydir)
    for i, word in enumerate(reversed(words)):
        #get the last part of the word seperated by an underscore, because this will be used for indexing
        last_word = re.split('_',word)[-1]
        num = re.sub('[^0-9.]','',last_word)
        if num:
           num = float(num)
           #apply a significant difference in weight, because the numbers can differ in order of magnitude
           #ie phi = 0.1 chiN = 40
           weight += num*(10**(5*(i+1)))
    return weight     



def slice2d(dirs):
    '''This function takes a 3d set of directories, and divides it into sets of 2d slices and boundary holders'''
    dirs.sort()
    dir3d,dir2d = [],[]
    for mydir in dirs:  #split the directories into 2d slices
         split = re.split('/',mydir)
         dir3d.append(split[0])
         dir2d.append(split[1]+'/'+split[2])
    unique3d = list(set(dir3d)) #find all the unique first directories
    unique3d.sort()
    sets2d = []
    for i in range(len(unique3d)): #match the 2d slice to the its coresponding value
       myset2d = []
       for j in range(len(dir3d)):
           if dir3d[j] == unique3d[i]:
                myset2d.append(dir2d[j])
       sets2d.append(myset2d)
    nodes = []
    pwd = os.getcwd()
    z = []
    for i in range(len(sets2d)):
        os.chdir(pwd+'/'+unique3d[i])
        z.append(float(re.sub('[^0-9.]','',re.split('_',unique3d[i])[-1])))#append the last set of numbers from the highest directory to the z list
        nodes.append(initialize_nodes(sets2d[i],args.filename))
        os.chdir(pwd)
    boundaryholders = []
    z = itertools.cycle(z)
    for mynodes in nodes:
        boundaryholders.append(calc_phase_boundaries(mynodes))
    present_boundaries = []
    for boundaryholder in boundaryholders:
        for boundary in boundaryholder.boundaries:
            present_boundaries.append(boundary.phaseA+'-'+boundary.phaseB)
    present_boundaries = list(set(present_boundaries)) #eliminate non unique entries 
    return boundaryholders,present_boundaries,z
        
# ==============================================================================
#   Begin Main
# ==============================================================================

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Tool to compute phase boundaries')
    parser.add_argument('-f', '--filename', action='store', default='F0_phases.dat',help='file that contains the phases and free energies at each phase point')
    parser.add_argument('-d', '--dirs', action='store', nargs='+', default=glob.glob("tau*/phiA*"),help='list of directories that contain each phase point')
    parser.add_argument('-o', '--outfig', action='store', default='',help='name of output figure file')
    parser.add_argument('--raw', action='store', default='',help='name of raw output file (for plotting in another program')
    parser.add_argument('-t', '--plottype', action='store', default='simplecolors',help='type of plot to generate')
    parser.add_argument('--xlabel', action='store', default=r"$f_A$",help='label for xaxis, can use \'$\' to write latex')
    parser.add_argument('--zlabel', action='store', default=r"$\tau$")
    parser.add_argument('--ylabel', action='store', default=r"$\chi N$",help='')
    parser.add_argument('--axisrange', action='store', nargs=4, default=[None,None,None,None],help='')
    parser.add_argument('--linecutoff', action='store', nargs='+', default=[1e30,1e30],help='maximum length of lines to draw in phase diagrams, useful to clean them up')
    parser.add_argument('-n','--dim',action='store',default=None,help='Number of dimensions to plot phase data in \n   1 => Free energy curves\n   2 => Phase Diagram \n   3 => 3d phase diagram  (guesses by default)')
    parser.add_argument('-i','--interp_dimension',action='store',default=[0],nargs='+',help='Dimensions to interpolate the phase diagram along ex: [0,1] would interpolate in 2 dimensions')
    parser.add_argument('-p','--plotstyle3d',action='store',default='flat',help='This argument changes the 3d plot style. Flat => multiple graphs with different linestyles on top of each other')
    parser.add_argument('--stylesheet',action= 'store',default=os.path.dirname(os.path.realpath(sys.argv[0]))+'/better_style.mplstyle',help='This argument is the Matplotlib stylesheet that will be used for graphing') #the default is located in the directory this script is located at
    parser.add_argument('--aspect',action='store',default=None,help='The aspect ratio for the outputted figure use 1 for a square fig, works for a 2d graph right now')
    parser.add_argument('-r', '--refphase', action='store', default=None,help='name of phase to reference to, only matters if 1d')
    print("IMPLEMENT CUSTOM AXIS RANGES AND LABELS FROM COMMAND LINE")
    args = parser.parse_args()
    
    #fnmeIn="F0_phases.dat"
    #dirs=glob.glob("tau*/phiA*");
    if os.path.isfile(args.stylesheet):
        plt.style.use(args.stylesheet)
        print('Graphing using the {} stylesheet'.format(args.stylesheet))
    else:
        print('WARNING: No stylesheet found at {} graphing with default Matplotlib settings'.format(args.stylesheet))
        print(os.path.dirname(os.path.realpath(sys.argv[0])))

    # Guess how many dimensions you're plotting
    if not args.dim:
        #Here it finds the where there are numbers in the first directory string, and then checks if those numbers are differnt in subsequent directory strings
        dirs=args.dirs


        folders = re.split('/',dirs[0])
        print('Trying to guess how many dimensions to plot in')
        args.dim = 0
        locs = [] # stores the indicies of directories in tree that contain numbers
        nums1 = [] # stores the actual numbers in each subdirectory of tree
        names=[]
        index = 0
        for folder in folders:
           if '*' in folder:
                raise ValueError("There is a '*' in you directories. Looks like a wildcard didn't get expanded. Check your path, something is likely wrong!")
           if re.search('[0-9]',folder):
               locs.append(index)
               first_num = re.sub('[^0-9.]',"",re.split('_',folder)[-1])#get the last number in the folder, sometimes there are multiple numbers
               nums1.append(float(first_num))
           index +=1
        for i in range(len(nums1)):
           #in every path if the number at the current location is different add one to the dimension
           for mydir in dirs[1:]:
               if mydir[-1] == '/':
                 raise RuntimeError(f'Directory "{mydir}" contains a trailing slash, this messes up parsing. Please remove and try again.')

               mysubdir = re.split('/',mydir)[locs[i]]
               if float(re.sub('[^0-9.]',"",re.split('_',mysubdir)[-1])) != nums1[i]:#grab the last number following the underscore since there may be multiple numbers in a directory
                    args.dim += 1
                    names.append(re.sub('[0-9,.]','',mysubdir))
                    break
        if args.dim == 0 or args.dim >3:
            raise ValueError('could not guess how many dimensions display in please specify with the -n flag')

        print('Graphing in {0} dimensions according to guess. If this is incorrect specify the number manually with the -n flag'.format(str(args.dim)))

        if args.dim == 1:
            args.xlabel = names[0]
            args.ylabel = 'Free Energy'
        elif args.dim == 2:
            args.ylabel = names[0]
            args.xlabel = names[1]
    else:
        args.dim = int(args.dim)


    if args.axisrange != [None, None, None, None]:
        args.axisrange = [float(x) for x in args.axisrange]
    print('axis ranges:',args.axisrange)
    print('xlabel:',args.xlabel)
    print('ylabel:',args.ylabel)
    print('zlabel:',args.zlabel)

    if args.interp_dimension != [0]:
        args.interp_dimension = [int(i) for i in args.interp_dimension]


    if args.dim == 1:
        nodes = initialize_nodes(args.dirs, args.filename)
        boundaryholder = calc_phase_boundaries(nodes)
        boundaryholder.plot(args.outfig,args.plottype, nodes=nodes,xlabel=args.xlabel, ylabel=args.ylabel,axisrange=args.axisrange,n=args.dim,aspect=args.aspect,refPhase=args.refphase)
        if args.raw != '':
            print("Saving free energy curve data to \'%s\'" % args.raw)
            boundaryholder.write(args.raw,dim=1)
    elif args.dim == 2:
        nodes = initialize_nodes(args.dirs, args.filename)
        boundaryholder = calc_phase_boundaries(nodes)
        dist_threshold = args.linecutoff
        boundaryholder.set_dist_thresholds(dist_threshold)
        os.getcwd()
        if args.raw != '':
            print("Saving raw phase boundary data to \'%s\'" % args.raw)
            boundaryholder.write(args.raw)
        print("Plotting with type \'%s\'" % args.plottype)
        boundaryholder.plot(args.outfig,args.plottype, nodes=nodes,xlabel=args.xlabel, ylabel=args.ylabel,axisrange=args.axisrange,aspect=args.aspect)

    elif args.dim == 3:
        if args.plotstyle3d == 'flat':
            boundaryholders,present_boundaries,z = slice2d(args.dirs)
            marker = itertools.cycle(('+', 'o', 's','v','^','<','>','x','p','h','H','*'))
            color  = itertools.cycle(('b','r','g','k','m','c'))
            #this dictionary used to differentiate when they were all grouped together
#            boundary_line_types = {}
#            for bound in present_boundaries:
#                boundary_line_types[bound] = {'color':next(color)}
#            boundary_line_styles = itertools.cycle(('-','--','-.',':'))
            fig,axlist = plt.subplots(nrows=len(present_boundaries),ncols=1)
            axdict = dict(zip(present_boundaries,axlist)) #make a dictionary relating a phase boundary to an axis
            bounddict = dict(zip(axlist,present_boundaries)) #make a dictionary with the opposite set of keys and values swapped for other uses
            #boundary
            for i,boundaryholder in enumerate(boundaryholders):
                my_z = next(z)
                boundaryholder.plot3d(axdict,'',marker=next(marker),color=next(color),z= my_z)
                
                if args.raw != '':
                    if i == 0:
                       with open(args.raw,'w+') as f:
                           f.write("#Raw phase boundary data for a layered phase boundary plot\n")
                           print("Saving raw phase boundary data to \'%s\'" % args.raw)
                    boundaryholder.write(args.raw,dim = 3,z=my_z)
               
            max_ax = [1e30,0,1e30,0]
            for ax in axlist:   #get furthest posible axis ranges and set them all equal
                   for i in [0,2]:
                       if (ax.axis()[i]) < (max_ax[i]):
                           max_ax[i] = ax.axis()[i]
                   for i in [1,3]:
                       if (ax.axis()[i]) > (max_ax[i]):
                           max_ax[i] = ax.axis()[i]
            if not all(args.axisrange):
                for ax in axlist:
                   #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                   ax.axis(max_ax)
                   ax.text(0.02,0.98,bounddict[ax],transform=ax.transAxes,fontsize='large',color ='k',horizontalalignment='left',verticalalignment='top',weight='semibold')#place text with the boundary type in the upper left corner
            for ax in axlist[0:-1]: ax.xaxis.set_ticklabels([]) #remove numbers but not axis ticks
            axlist[-1].set_xlabel(args.xlabel)
            axlist[int(len(axlist)/2)].set_ylabel(args.ylabel) #put the y label in the center of the plots
            for ax in axlist: ax.yaxis.set_ticks([2,4,6])#FORCE Y AXIS TICKS REMOVE IF YOU WANT DEFAULT NUMBERS
            filename= args.outfig
            if args.outfig == '':
                plt.gcf().set_tight_layout(True) #makes the graph look better by removing whitespace
                plt.show()
            else:
                if args.aspect:
                     ratio = float(args.aspect)
                     xleft, xright = ax.get_xlim()
                     ybottom, ytop = ax.get_ylim()
                     # the abs method is used to make sure that all numbers are positive
                     # because x and y axis of an axes maybe inversed. 
                     # convert to figure coords from data coords
                     for ax in axlist:
                         ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)
                if filename.split('.')[-1] == 'eps':
                    plt.savefig(filename, bbox_inches='tight', format='eps')
                else:
                    plt.savefig(filename, bbox_inches='tight')
        else:
             print("ERROR: unrecognized 3d plot type arguement")
    else:
        print("Error dimension argument must be between 1 and 3")




