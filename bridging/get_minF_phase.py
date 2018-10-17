#!/usr/bin/python

import argparse
import glob,re
import os.path

parser = argparse.ArgumentParser(description='PolyFTS Bridging Operator')
parser.add_argument('-d', '--dir',   action='store', nargs='?', type=str,required=True,help='directory to calculate phases in')
parser.add_argument('-f', '--file',   action='store', nargs='?', type=str,default='minphase.dat',help='file to write minphase in')
parser.add_argument('-i', '--ignore',action='store', nargs='+', type=str, default=['SIGMA', 'A15', 'C14', 'C15', 'GYR'],help='Phases to ignore when calculating minF phase')
parser.add_argument('-o', '--output',action='store_true', default=False,help='save result to text file in the directory of interest')
args = parser.parse_args()

#command line args
phases_to_ignore=args.ignore
wdir=args.dir
ofile="%s/%s" % (wdir,args.file)

dirs=glob.glob("%s/*Phase" % wdir );

phases=[]
F=[]
for mydir in dirs:
    phase=re.sub('Phase','',mydir.split('/')[-1])
    if phase not in phases_to_ignore:
        

        phases.append(phase)

        #find min intensive hamiltonian
        logfile="%s/%s.out" % (mydir,phase)
        if not os.path.isfile(logfile):
            print "Skipping %s, %s does not exist" % (phase,logfile)
            continue
        
        Ffinal = False
        with open(logfile,'r') as f:
            for line in f:
                if "Intensive Hamiltonian" in line:
                    Ffinal=float(line.split()[3])

        if Ffinal == False:
            raise ValueError('Intensive Hamiltonian not found in %s' % logfile)

        F.append(Ffinal)


# now find min phase
val, idx = min((val, idx) for (idx, val) in enumerate(F)) 

# check if all are DIS
THRESH=1e-4
alldis=True
for myF in F:
    if (abs(myF-val) > THRESH):
       alldis=False 
if alldis:
   idx=phases.index('DIS') 

# and output to file
if args.output:
    with open(ofile,'w') as f:
        f.write(phases[idx])
else:
    print phases[idx]

