#!/usr/bin/env python3

import numpy as np
import pdb

'''
This is a utility script that postprocesses a bridging calculation
 - Its main job is to integrate BridgingOperator.dat to get Pmtilda and other useful quantities
 - It also outputs the bridging fraction computed from Pm_tilda and Pm_hat
 - note that BridgingOperator.dat must be in ascii and 'Format verson 3'

Questions? Contact Joshua Lequieu (lequieu@mrl.ucsb.edu)
'''

def parse_file(filename):
    ''' Check BridgingOperator.dat to make sure that its
    '''
    with open(ifilename,'r') as f:
        for line in f:
            l=line.split()
            if 'Format version' in line:
                formatversion=float(l[-1])
                if formatversion != 3:
                    print("Error! Format version of %s must be 3" % ifilename)
                    exit
            elif 'NDim' in line:
                ndim=int(l[-1])
            elif 'nfields' in line:
                nfields=float(l[-1])
            elif 'k-space' in line and 'complex' in line:
                if "k-space data = 0 , complex data = 0" not in line:
                    print("Error! %s must have k-space data = 0, complex data = 0" % ifilename)
                    exit
    data = np.loadtxt(ifilename)

    return data,ndim,nfields         

ifilename = "BridgingOperator.dat"
ofilename = "Pm.dat"
d,ndim, nfields = parse_file(ifilename)

npw = d.shape[0]
idxtargetregion = d[:,ndim+1][d[:,ndim+0]==1][0]


# define boolean arrays to select different regions (using dimension to grab correct column)
istargetregion = d[:,ndim+0]
isnottargetregion = 1-istargetregion
isanyregion = (d[:,ndim+1] != 0)
isanyregionbuttarget = (isanyregion != istargetregion )

# ----------------------------------------------
# first calc Pbar and Pbar_prime
# ----------------------------------------------

col_offset = ndim + 2 # +2 is for isregionid and then regionidx
# M is total number of arms
M = d.shape[1] - col_offset #d.shape[1] - col_offset

comment=''' 
Pm_storage stores the integrated (or marginal) distributions of Pm(r)
There are 5 different columns, columns 2-5 only differ by how the integral is performed
1) The value of m (the number of arms) constrained in region R, runs to M (the total number of di-block arms)
2) P = \int_{R} dr P_m(r) 
3) P^\prime = \int_{not R} dr P_m(r) 
###4) P + P^\prime = \int_{all space} P_m(r) **This is Spencer2017 Equation 12**
4) Ptilda (see Josh's CFDC presentation) - roughly gives the probability of m arms bridging the two domains
    - constrains the star center to be outside R, if an arm is in R
    in the limit that domains fill space, Ptilda becomes Phat
5) Phat (see Josh's CFDC presentation)
'''
Pm_storage = np.zeros((M+1,5))

total = 0
for m in range(1,M+1):
    col = m + col_offset - 1
    P       = np.sum(istargetregion * d[:,col]) / npw       
    Pprime  = np.sum(isnottargetregion * d[:,col]) / npw
    total += P + Pprime

    Pm_storage[m,0] = m
    Pm_storage[m,1] = P         
    Pm_storage[m,2] = Pprime
assert(np.allclose(total,1.0)), "Integral over P and Pprime != 1"

# ----------------------------------------------
# then calc Ptilda and Phat
# ----------------------------------------------
sumPtilda=0
sumPhat=0
# b here gives the number of strands that 'bridge' the B domain
for b in range(M+1):
    colb = b + col_offset - 1
    colMminusb = (M-b-1) + col_offset
    if b==0:
        Ptilda = np.sum(istargetregion * d[:,colMminusb])/npw 
        Phat = np.sum(istargetregion * d[:,colMminusb])/npw 
    elif b > 0 and b < M:
        colMminusb = (M-b-1) + col_offset
        Ptilda = np.sum(istargetregion * d[:,colMminusb])/npw + np.sum(isnottargetregion * d[:,colb])/npw
        Phat = np.sum(istargetregion * d[:,colMminusb])/npw + np.sum(isanyregionbuttarget * d[:,colb])/npw
    else:
        Ptilda = np.sum(isnottargetregion * d[:,colb])/npw
        Phat = np.sum(isanyregionbuttarget * d[:,colb])/npw
    if b != 0:
        sumPtilda += Ptilda
        sumPhat += Phat
    Pm_storage[b,3] = Ptilda
    Pm_storage[b,4] = Phat

print("Bridging Fraction (1-PM): %f" % (1-(Pm_storage[-1,2]+Pm_storage[-1,1])))
print("Bridging Fraction (\sum Ptilda): %f" % sumPtilda)
print("Bridging Fraction (\sum Phat): %f" % sumPhat)
#header="# m Pmbar Pmbar_prime Ptilda Phat"
np.savetxt(ofilename,Pm_storage,header=comment, fmt='%.10f')
