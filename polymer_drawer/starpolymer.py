#!/usr/bin/env python3
import numpy as np
import pdb
import matplotlib.pyplot as plt
import argparse,os,sys,glob,re
pi = np.pi



class star():
    def __init__(self,filename,cwd):
        with open(filename,"r") as f:
            isstar = False
            self.armlist = []
            arm_index = -1	
            for line in f:
                if "Initializing propagator for a STAR (multi-block) polymer" in line: 	#check for a star polymer	
                    isstar = True
                if not isstar: #when there is a star polymer check weight percents
                    continue
                
                if "Arm index" in line:
                    arm_index +=1
                    self.armlist.append(arm(self))
                if arm_index != -1:
                    if 'Multiplicity' in line: #add the self.multiplicity to the list	
                        self.armlist[arm_index].multiplicity = int(re.sub("[^0-9]","",line))#cut out the interger self.multiplicity from the line
                    elif 'Arm length,' in line:
                        self.armlist[arm_index].length = float(re.sub("[^.0-9]","",line))#cut out the decimal armlength from the line
                    elif 'Species of blocks' in line:
                        words = re.split(' ',line)
                        spec = []
                        for word in words:
                            if re.search('[0-9]',word):#check if the word is a number
                                spec.append(int(word)-1)#account for indexing from 1 in POLYFTS
                        self.armlist[arm_index].species = spec
                    elif 'Block fractions' in line:
                        words = re.split(' ',line)
                        fracs = []
                        for word in words:
                            if re.search('[0-1]\.[0-9].*',word):#check if the word is a number
                                fracs.append(float(word))
                        self.armlist[arm_index].block_frac = fracs	
            if not isstar:
                raise TypeError('The file at {} is not a star polymer'.format(cwd+'/'+filename))
            #number of diblock arms for graphing with respect to narms and phi
            self.n_diblock=self.armlist[1].multiplicity
            #initilize the species amounts to all 0
            species_amounts = [0]*(max([max(a.species) for a in self.armlist])+1)
            for myarm in self.armlist:
                for spec,bf in zip(myarm.species,myarm.block_frac):
                    species_amounts[spec] = bf*myarm.length
            phiA = species_amounts[0]/sum(species_amounts)
            self.phi = round(phiA,2)#round for aesthetics

class arm:
    def __init__(self,mystar):
        self.length = None
        self.block_frac = None
        self.multiplicity = None
        self.species = None
        self.x = []
        self.y = []
        self.star = mystar
        self.angle = []
        self.sectionsize = []
        self.arclength = []
    def draw_arm(self,ax,angles,sectionsize,colors):
        ''' a function that draws arms, takes an arm as input, as well as the angle to rotate the arm
            a list of colors and the section size that is the angle that the section of the unit circle
            that the arm can take up'''
        numsins = 20
        ax.axis('equal')#need to be equal for proper visual scaling
        n = self.length
        slope = np.tan(sectionsize/2)
        #make the polymers extra long, so that the section that is the right arclength can be cut off
        #this is not the most robust solution but it is quick and probably wont fail
        for angle in angles:
            attempts = 0
            intersecting =True
            #loop while checking for intersections, due to randomness
            #and the generation this cant be avoided
            while intersecting:
                x = np.linspace(0,n*5,10000)
                w = 1.0/(np.random.rand(numsins)*3*n/10 +n/10 )
                phase_shift = np.random.rand(numsins)*np.pi
                y = np.zeros(np.shape(x))
                for w, phase_shift in zip(w,phase_shift):
                    sign = [-1,1][np.random.randint(2)]#Add a random sign for more random appearance
                    y += np.cos(x*w+phase_shift)*sign
                y = y/np.max(np.abs(y))*x
                arclength = [0]
                i=1
                while arclength[-1] < n:
                     arclength.append(arclength[-1]+np.sqrt((y[i]-y[i-1])**2+(x[i]-x[i-1])**2))
                     i+=1
                x,y = x[0:i],y[0:i]
                xp = x*np.cos(angle)-y*np.sin(angle) #apply a rotation matrix to each coordinate of the polymer
                yp = x*np.sin(angle)+y*np.cos(angle)
                intersecting = self.check_intersect(xp,yp)
                if attempts > 15:
                    self.redraw(ax,angles,sectionsize,color)
                    return
                attempts += 1 
            self.x.append(xp)
            self.y.append(yp)
            self.arclength.append(arclength)
        sections = (1-np.cumsum(self.block_frac))*self.length#use the block fractions to divide into arclength sections the 1- is  because they are reversed in polyfts
        sections = np.flip(sections,0)
        self.species.reverse()
        for x,y,arclength in zip(self.x,self.y,self.arclength):
            arclength = np.array(arclength)
            for i in range(len(sections)):
                if i != len(sections)-1:
                    insection = np.logical_and(arclength>=sections[i],arclength<sections[i+1])
                else:
                    insection = np.logical_and(arclength>=sections[i],arclength<=self.length)
                ax.plot(x[insection],y[insection],color=colors[self.species[i]],lw=3)


    def check_intersect(self,xp,yp):
        front_cutoff = int(len(xp)/3) #the beginning of the arms wont intersect so ignore them in intersection checks
        for myarm in self.star.armlist:
            for x,y in zip(myarm.x,myarm.y):
                for myx, myy in zip(xp[front_cutoff:],yp[front_cutoff:]):
                    #check intersection between point and other polymer chain
                    #if intersecting make a new random arm
                    if np.any(np.logical_and(np.abs(x-myx) < .005,np.abs(y-myy) < .005)):
                        return True 
        return False
    def redraw(self,ax,angles,sectionsize,color):
        self.x,self.y,self.arclength= [],[],[]
        self.draw_arm(angles,sectionsize,color)















parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dirs', action='store',required=True, nargs='+',help='list of directories that contain each .out file to be graphed')
args = parser.parse_args()

dirs=args.dirs
dirs.sort()
print("Dirs: {}".format(dirs))
idir=os.getcwd()
starlist = []
for mydir in dirs:
    os.chdir(mydir)
    starlist.append(star(glob.glob('*.out')[0],idir))#create an arm for each .out file (you should only have one in each directory)
        

    os.chdir(idir)
colorlist = ['r','b','c','m','k','g','y']
n_arms = []
phi = []
for mystar in starlist:
    n_arms.append(mystar.n_diblock)
    phi.append(mystar.phi)
n_arms,phi = list(set(n_arms)),list(set(phi))#find unique n and phi
n_arms.sort()#need to resort since set removed order
phi.sort()
f,axlist = plt.subplots(len(n_arms),len(phi))
#make dictionaries that relate phi to a specific subplot
phi_dict = dict([(phi,i) for i, phi in enumerate(phi)])
n_dict = dict([(n_arms,i) for i,n_arms in enumerate(n_arms)])
for mystar in starlist:
    #get the subplot to graph on
    ax_index =  [n_dict[mystar.n_diblock],phi_dict[mystar.phi]]
    if len(np.shape(axlist)) == 1:#in case of a 1 dimensional plot
        ax_index = sum(ax_index)
        ax = axlist[ax_index]
    else:
        ax = ax_list[ax_index[0],ax_index[1]]
    
    #FUTURE IMPROVEMENT:
    #this is a hardcoded style for graphing miktoarms, this can be generalized for any star polymer
    #but that is not worthwile at the moment
    mystar.armlist[0].draw_arm(ax,[-pi],pi/2,colorlist)
    diblock_spacing = np.linspace(-pi/2,pi/2,mystar.armlist[1].multiplicity+2)[1:-1] #divide the starting angle into multiplicity plus two sections so endpoints can be removed
    sectionsize = pi/mystar.armlist[1].multiplicity
    mystar.armlist[1].draw_arm(ax,diblock_spacing,sectionsize,colorlist)
plt.axis('off')
plt.show()

