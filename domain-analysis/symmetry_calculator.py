#!/usr/bin/env python3
import subprocess
import pdb
import itertools
import numpy as np
#import pymatgen as mg #this is a materials science library
import os,sys,argparse
import re
script_dir = os.path.realpath(sys.argv[0])#the absolute path to this script
#allow importing from the lib folder in a parallel directory
sys.path.append('/'.join(re.split('/',script_dir)[0:-2])+'/lib')
script_dir = '/'.join(re.split('/',script_dir)[0:-1])#remove the script from the script dir
from domaintools import DomainAnalyzer
import iotools as io
#from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
#from pymatgen.io.cifio
#structure = mg.Structure.from_file('structure.cif')
#spg = SpacegroupAnalyzer(structure,symprec=0.05)
#print(spg)
#print(spg.get_space_group_symbol())


# ==============================================================================
#   Begin Main
# ==============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tool to guess crystallographic symmetry based on density.bin files. Requires a density.bin or dat file in the directory this script is run from. Also requires Platon to be compiled two directory levels higher than this script' )
    parser.add_argument('-f', '--filename', action='store', default='density.bin',help='file that contains the density information: default = density.bin')
    parser.add_argument('--round',action='store',default=4,help='The amount of decimal places to round the volumes and center of mass locations tothis is important for uniqueness calculations. default =4')
    parser.add_argument('-t','--tolerance',action='store',default=0.25,help='How tolerant platon will be to domains that are slightly out of place Default = 0.25 radii of gyration')
    parser.add_argument('-i','--ignore',action='store',default=0,help=' minimum volume for domains, domains with smaller volumes will be ignored')
    #parser.add_argument('-o','--overwrite',action='count',help='Use with this option to overwrite .cif files generated by previous meshing')
    args = parser.parse_args()
    args.round = int(args.round)
    args.tolerance = float(args.tolerance)
    #the default tolerance in platon is 0.5 angstroms (radii of gyration)
    #in order to change this expand the box by a factor, so the relative 
    #tolerance is tighter
    box_multiplier = 0.5/args.tolerance
#    if not os.path.isfile('structure.cif') and not args.overwrite:
#        make_domains(args.filename)
#    else:
#        print('Found an existing structure.cif skipping domain calculation, disable this feature with -o')

    #======================================================================
    #this section finds the domains, using domaintools with meshing enabled
    #======================================================================
    filetype = re.split('\.',args.filename)[-1]
    if filetype =='bin':
        AllCoords, AllFields = io.ReadBinFile(args.filename)
    elif filetype =='dat':
        AllCoords, AllFields = io.ReadDatFile(args.filename)
    else:
        raise NotImplementedError('This only supports .bin and .dat files')
    #normalize by box length
    #extra AllCoords[1,1,1] is because boxes are mesured by lower left corner
    box_length = AllCoords[-1,-1,-1] + AllCoords[1,1,1]
    AllCoords = AllCoords/(box_length)
    domainanalyzer = DomainAnalyzer(AllCoords,AllFields)
    ndomains, com,area,vol,IQ = domainanalyzer.getDomainStats(plotMesh=False,add_periodic_domains=True)
    domainanalyzer.meshAllDomains(plotfile="mesh_all.png")
    stats = np.vstack((com.T,area,vol,IQ)).T
    np.savetxt("stats_mesh.dat",stats, header="comx comy comz area vol IQ")
    
    #======================================================================
    #this section makes a crystalography file for platon (.cif) 
    #======================================================================
    args.ignore = float(args.ignore)
    vol = np.around(vol,decimals=args.round)#numbers wont be exact so round them
    not_too_small = vol >= args.ignore
    com = np.around(com,decimals=args.round)
    vol = vol[not_too_small]
    com = com[not_too_small]
    unique_vol = list(set(vol))#create a list of the unique volumes
    unique_vol.sort()
    volume = np.prod(box_length)
    #expand the box to change tolerances
    box = itertools.cycle(box_length*box_multiplier)
    #a list of elements to differntiate different sizes of domains for symmetry
    #platon can do som weird bonding guesses, so start by adding unreactive elements, then transition metals that form frankcastper phases
    elements = ['He','Ne','Ar','Kr','Xe','Ag','Au','Pt','Sc','Ti','V','Cr','Mn','Fe','Ni','Cu','Zn','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Cd']
    #elements = ['Ag','Au','Cu','Pt','Ni']
    print('There are {} unique volumes. If this seems too high or low change the rounding parameter with --round'.format(len(unique_vol)))
    if len(elements)>=len(unique_vol):
        elementdict = dict(zip(unique_vol,elements))
    else:
        raise ValueError("There are more unique volumes than elelments in this file's element list: {}  Either increase the rounding on volume with --round(likely the problem, there should not be too many unique domains) or add more to the element list".format(script_dir))
    #this string describes a basic crystalography file
    base = '''
data_global
_cell_length_a 
_cell_length_b 
_cell_length_c 
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_cell_volume 

loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
'''
   # with open(script_dir+'/base.cif') as base:#open a base cif file to edit
    with open('structure.cif','w+') as cif:
        for line in base.splitlines():
            if '_cell_length_' in line:
                cif.write(line.rstrip()+' '+str(next(box))+'\n')#add cell lengths if that is in the line rstrip removes newline
            elif '_cell_volume' in line:
                cif.write(line.rstrip()+' '+str(volume)+'\n')
            else:
                cif.write(line+'\n')
        for line,vol in zip(com[:],vol):
            cif.write(elementdict[vol]+"   " + "   ".join([str(word) for word in line])+'\n') #add the locations of the centers of mass to the cif file, note that I am adding them as helium because platon doesnt seem to try weird bonding stuff with helium            
        cif_path = cif.name #the absolute path to the cif file to ensure no errors from platon 
    #======================================================================
    #this section runs platon and parses the file
    #======================================================================
    #run platon and discard the output, because it saves an output file automatically
    #PLATON MUST BE INSTALLED 2 DIRECTORIES HIGHER FOR THIS CALL TO IT TO RUN
    subprocess.call(['/'.join(re.split('/',script_dir)[0:-2])+'/platon/platon','-o','-m',cif_path],stdout=subprocess.DEVNULL)
    if not os.path.isfile('structure.lis'):
        raise FileNotFoundError("structure.lis in the local direcory, maybe platon did not create it?")
    spacegroup,precision = [],[]
    #only check for space groups after a it says it has found a new symmetry
    check4spacegroups = False
    with open('structure.lis') as f:
        for line in f:
            if 'Conventional, New or Pseudo Symmetry' in line:
                check4spacegroups = True
            if check4spacegroups and 'Space Group  H-M:' in line:
                spacegroup = re.split(' ',line)[5]
                #for i,word in enumerate(re.split(' ',line)):
                 #   print(i,word)
            if ':: Change of Crystal System indicated.' in line:
                for word in re.split(' ', line):
                    if re.search('[0-9]',word):
                        precision = float(word)/box_multiplier
                        if float(word) != 0.5:
                            print("WARNING: PLATON'S MAXIMUM ACCEPTABLE DEVIATION WAS NOT FOUND TO BE 0.5 ANGSTROMS, SO THE TOLERANCE SPECIFIED IN RADII OF GYRATION IS WRONG. IN ALL TEST CASES THE MAXIMUM ACCEPTABLE DEVIATION WAS 0.5 RADII OF GYRATION SO THIS SEEMS TO BE PLATON'S DEFAULT, BUT IN CASE ITS NOT I ADDED THIS WARNING.  READ THE structure.lis FILE TO SEE PLATON'S OUTPUT")
    if spacegroup and precision:
        print('Platon found the spacegroup to be {} with a maximum acceptable deviation of {} times the radius of gyration \n For more detailed output open the structure.lis file'.format(spacegroup,precision))
    else:
        print("ERROR: platon did not find a space group, check the structure.lis file for Platon's raw output")

