#!/usr/bin/env python

import matplotlib.pyplot as plt
import argparse
import re
import numpy as np
import glob
import pdb



parser = argparse.ArgumentParser()
#parser.add_argument('chiN_dir',type=str)
parser.add_argument('-d', '--dirs', action='store',required=True, nargs='+',help='list of directories that contain each phase point')
parser.add_argument('-r', '--refphase', action='store', default='',help='name of phase to reference to')
parser.add_argument('-f', '--filename', action='store', default='F0_phases.dat',help='file that contains the phases and free energies at each phase point')
args = parser.parse_args()

if args.refphase == '':
    subtract_ref = False
else:
    subtract_ref = True
    refphase=args.refphase
    if "Phase" not in refphase:
        refphase += "Phase"

filename=args.filename

dirs=args.dirs
dirs.sort()
print dirs

#chiN_dir=args.chiN_dir
#phiA_dirs=glob.glob("%s/phiA*" % chiN_dir)
#phiA_dirs.sort()
#print phiA_dirs



phases=[]
npoints=0
# get max phases
dir1s = []
dir2s = []
for mydir in dirs:
    #mydir=re.sub("./","",mydir) # remove front ./ if exists
    fnme="%s/%s" % (mydir,filename)
    dir1s.append(mydir.split("/")[0])
    dir2s.append(mydir.split("/")[1])
    
    with open(fnme,'r') as f:
        for line in f:
            l = line.split()
            phase_tmp= l[0]
            F_tmp = float(l[1])
            if phase_tmp not in phases:
                phases.append(phase_tmp)
        npoints += 1

# check to see if its dir1 or dir2 thats varrying
if dir1s.count(dir1s[0]) == len(dir1s) and  dir2s.count(dir2s[0]) != len(dir2s):
    dir_that_varies=1   #dir2 varries, use this on x axis
elif dir1s.count(dir1s[0]) != len(dir1s) and  dir2s.count(dir2s[0]) == len(dir2s):
    dir_that_varies=0   #dir2 varries, use this on x axis
else:
    print "Error! It loosk like both dir1s and dir2s vary! the arguement of -d is wrong"
    print "dir1s: ",dir1s
    print "dir2s: ",dir2s


xpos=np.zeros(npoints)
data=np.full((npoints,len(phases)),np.nan)
# now populate data
i=0
for mydir in dirs:
    #mydir=re.sub("./","",mydir) # remove front ./ if exists
    fnme="%s/%s" % (mydir,filename)
    with open(fnme,'r') as f:
        
        xpos_tmp=re.sub("[a-z,A-Z]","", mydir.split("/")[dir_that_varies])
        if '_' in xpos_tmp: # bit of a hack for tau0.9_chiN13.0 directories
            xpos_tmp=xpos_tmp.split('_')[-1]
        xpos[i] = float(xpos_tmp)
        for line in f:
            l = line.split()
            phase_tmp=l[0]
            F_tmp = float(l[1])
            phase_idx = phases.index(phase_tmp)
            data[i,phase_idx] = F_tmp
        i += 1


# custom check for DIS boundaries
if subtract_ref:
    refphase_idx = phases.index(refphase)

fig = plt.figure()
ax = fig.add_subplot(111)
for i in range(len(phases)):
    x = xpos
    if subtract_ref:
        y = data[:,i] - data[:,refphase_idx]
    else:
        y = data[:,i]
    mask=np.isfinite(y)
    ax.plot(x[mask],y[mask],label=phases[i],ls='-',marker='.')
#ax.locator_params(nticks=10)
start, end = ax.get_xlim()
#ax.xaxis.set_ticks(np.arange(0.1, 0.9,0.1 ))
#ax.axis([None,None,None,0.1])
ax.legend()
ax.grid()
plt.show()
