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


def printStatusWarning(status,status2ignore,wdir):
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



def checkPhi(filename,wdir):
	'''checks if the first component has the proper mole fraction, and that errors are not caused by rounding due to the contour_ds. This assumes that the .out file structure follows the same format that I am working with.  It should work for most star polymers'''        
	phistr = re.split('/',wdir)
	for i in phistr:    #this splits the directory and finds where phis is followed by numbers then extracts that number
		if re.search('phiA[0-9].*',i):
			phi = float(re.sub("[a-z,A-Z]","",i))
	if not phistr:
		print('WARNING: no mole fraction found in {} not checking if contour_ds is small enough'.format(wdir))
		return True
	with open(filename,"r") as f:
		isstar = False
		arm_lengths = []
		arm_index = -1	
		block_frac = []
		multiplicity = []
		species = []
		for line in f:
			if "Initializing propagator for a STAR (multi-block) polymer" in line: 	#check for a star polymer	
				isstar = True
			if not isstar: #when there is a star polymer check weight percents
				continue	
			arm = False
			if "WARNING: Arm type" in line and "length rounded to" in line:
				words = re.split(" ",line) #Get each seperate word from the line
				for word in words:
					if re.search('[0-9].*',word):# if the word is actually a number
						if not arm:#when arm hasnt been defined the first number is the arm
							arm = True
						else:#second number is the arm length add it to the arm length list
							arm_lengths.append(float(word))
			if not arm_lengths: #check if there were warnings at all
				continue

			if "Arm index" in line:
				arm_index +=1
			if arm_index != -1:
				if 'Multiplicity' in line: #add the multiplicity to the list	
					multiplicity.append(int(re.sub("[^0-9]","",line)))#cut out the interger multiplicity from the line
				elif 'Species of blocks' in line:
					words = re.split(' ',line)
					spec = []
					for word in words:
						if re.search('[0-9]',word):#check if the word is a number
							spec.append(int(word)-1)#account for indexing from 1 in POLYFTS
					species.append(spec)
				elif 'Block fractions' in line:
					words = re.split(' ',line)
					fracs = []
					for word in words:
						if re.search('[0-1]\.[0-9].*',word):#check if the word is a number
							fracs.append(float(word))
					block_frac.append(fracs)	
		if not arm_lengths:
			if verbose:
				print('In the file {} there were no warnings so phi is correct'.format(filename))
			return True #no warnings in whole file so values must be good
		species_amounts = [0]*(max(max(species))+1)#initialize the species amounts to the number of species there are
		
		for armnum in range(arm_index+1):
			for specnum in range(len(species_amounts)):#iterate over every possible species
				for myspec in range(len(species[armnum])):
					if species[armnum][myspec] == specnum:
						try:
							species_amounts[specnum] += multiplicity[armnum]*arm_lengths[armnum]*block_frac[armnum][myspec]#if the species is the one looped over add the ammount of it to the list
						except IndexError:
							print('ERROR: in the file {} one of the fields was not found. What was found is listed below: \n Arm lengths: {} \n Multiplicity: {} \n Species: {} \n Block Fractions {} \n Currently iterating over arm {} and species {}'.format(filename,arm_lengths,multiplicity,species,block_frac,armnum,myspec))
		real_phi = species_amounts[args.phi_species_num]/sum(species_amounts)	
		intol =  abs(phi-real_phi) <= args.fraction_tolerance#check if phi is within tolerance
		if not intol:
			print('WARNING: In the file {} countour_ds was not small enough to resolve the polymer to the desired mole fraction.  The desired phi was {}, and the actual phi was {}. The tolerance is was {}. Skipping!'.format(filename,phi,real_phi,args.fraction_tolerance))
		if verbose:
			print("In the file {} the actual phi was found to be {}, and the desired phi was {}".format(filename,real_phi,phi))
		return intol	


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dirs', action='store',required=True, nargs='+',help='list of directories that contain each phase point')
parser.add_argument('--ignorephase', action='store',default=[''], nargs='+',help='phases to ignore')
parser.add_argument('--ignorestatus', action='store', default=[1,3], nargs='+',help='status to ignore')
parser.add_argument('-f', '--filename', action='store', default='F0_phases.dat',help='file that contains the phases and free energies at each phase point')
parser.add_argument('-t','--fraction_tolerance',action='store',default=0.005,help='Tolerance on rounding of the weight fraction, values where the rounding is greater than this value will be excluded')
parser.add_argument('-p','--phi_species_num',action='store',default=0,help='The species number used to check that the number fraction (phiA) is correct. The default is 0 (the first species) enter -1 to skip this check')
parser.add_argument('-v','--verbose',action='count',help='Make the ouput of the script verbose, will show found phi values')
parser.add_argument('--writemin', action='store_true', help='write minimum phase to file')
args = parser.parse_args()

filename=args.filename
phases2ignore = [a+"Phase" for a in args.ignorephase ]
status2ignore = [int(i) for i in args.ignorestatus]
dirs=args.dirs
dirs.sort()
if args.verbose:
	verbose = True
else:
	verbose = False
print("Dirs: {}".format(dirs))
print ("Status To Ignore: {}".format(status2ignore))
print ("Phases To Ignore: {}".format(phases2ignore))
#need to convert the given strings into integers, this was missing
status2ignore = [int(s)for s in status2ignore]
args.fraction_tolerance = float(args.fraction_tolerance)
args.phi_species_num=int(args.phi_species_num)
idir=os.getcwd()

for mydir in dirs:
    os.chdir(mydir)
    phases=glob.glob("*Phase")
    with open(filename,"w") as f:
        phaseslist=[]
        F0list=[]
        for phase in phases:
            wdir=idir+'/'+mydir+"/"+phase
            shortphase=re.sub("Phase","",phase)
            if phase not in phases2ignore:
                if checkStatus(wdir,status2ignore):
                    if args.phi_species_num == -1 or checkPhi("{}/{}.out".format(wdir,shortphase),wdir):
                        F0 = parseLogForF0("{}/{}.out".format(wdir,shortphase))
                        if F0 != 0.0:
                            f.write("{} {}\n".format(phase,F0))
                            phaseslist.append(phase)
                            F0list.append(float(F0))

    #write minphase.dat if flag is set
    if args.writemin:
        # now find min phase
        val, idx = min((val, idx) for (idx, val) in enumerate(F0list))

        # check if all are DIS
        THRESH=1e-4
        alldis=True
        for myF in F0list:
            if (abs(myF-val) > THRESH):
               alldis=False
        if alldis:
           idx=phaseslist.index('DIS')

        # and output to file
        ofile="minphase.dat"
        with open(ofile,'w') as f:
            phase=phaseslist[idx]
            shortphase=re.sub("Phase","",phase)
            f.write("{}".format(shortphase))

    os.chdir(idir)
