#!/usr/bin/env python3

import argparse
import os,glob,re

def checkStatus(wdir,status2ignore):
    filename="{}/STATUS".format(wdir)
    try:
        with open(filename,'r') as f:
            status=int(f.readline())
    except FileNotFoundError as e:
        #print (e)
        print("{} not found...skipping!".format(filename))
        return False

    if status == 0 and 0 in status2ignore:
        print("{} is not converged (killed while running) ...skipping!".format(wdir))
        return False
    if status == 1 and 1 in status2ignore:
        print("{} is divergent ...skipping!".format(wdir))
        return False
    if status == 3 and 3 in status2ignore:
        print("{} is not converged (reached max steps) ...skipping!".format(wdir))
        return False

    return True


def printStatusWarning(status,status2ignore,widr):
    if status == 0 and 0 in status2ignore:
        print("{} is not converged (killed while running) ...skipping!".format(wdir))
    if status == 1 and 1 in status2ignore:
        print("{} is divergent ...skipping!".format(wdir))
    if status == 3 and 3 in status2ignore:
        print("{} is not converged (reached max steps) ...skipping!".format(wdir))

def parseLogForF0(filename):
    with open(filename,"r") as f:
        found=False
        for line in f:
            if "Intensive Hamiltonian" in line: 
                found=True
                F0=float(line.split()[3])
        if found==True: 
            return F0
        else:
            print("Warning! No Intensive Hamiltonian in {}".format(filename))
            return 0.0
                 
        

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dirs', action='store',required=True, nargs='+',help='list of directories that contain each phase point')
parser.add_argument('--ignorephase', action='store',default=[''], nargs='+',help='phases to ignore')
parser.add_argument('--ignorestatus', action='store', default=[1,3], nargs='+',help='status to ignore')
parser.add_argument('-f', '--filename', action='store', default='F0_phases.dat',help='file that contains the phases and free energies at each phase point')
args = parser.parse_args()

filename=args.filename
phases2ignore = [a+"Phase" for a in args.ignorephase ]
status2ignore = args.ignorestatus
dirs=args.dirs
dirs.sort()
print("Dirs: {}".format(dirs))
print ("Status To Ignore: {}".format(status2ignore))
print ("Phases To Ignore: {}".format(phases2ignore))
idir=os.getcwd()

for mydir in dirs:
    os.chdir(mydir)
    phases=glob.glob("*Phase")
    with open(filename,"w") as f:
        for phase in phases:
            wdir=idir+'/'+mydir+"/"+phase
            shortphase=re.sub("Phase","",phase)
            if phase not in phases2ignore:
                if checkStatus(wdir,status2ignore):
                    F0 = parseLogForF0("{}/{}.out".format(wdir,shortphase))
                    if F0 != 0.0:
                        f.write("{} {}\n".format(phase,F0))
    os.chdir(idir)
        
