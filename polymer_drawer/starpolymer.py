#!/usr/bin/env python3
import numpy as np
import pdb
import matplotlib.pyplot as plt
import argparse,os,sys,glob,re
import pickle
pi = np.pi



class star():
    def __init__(self,filename,cwd):
        self.lines = []
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
           # species_amounts = [0]*(max([max(a.species) for a in self.armlist])+1)
           # for myarm in self.armlist:
           #     for spec,bf in zip(myarm.species,myarm.block_frac):
           #         species_amounts[spec] = bf*myarm.length
           # phiA = species_amounts[0]/sum(species_amounts)
           # self.phi = round(phiA,2)#round for aesthetics
           #figure out what phi A is for graphing
            words =  re.split('/',cwd)
            for word in words:
                if 'phiA0' in word:
                    self.phi = float(re.sub('[^0-9.]','',word))
            if not self.phi:
               raise ValueError('Could not find phi in directory name')
      
    def draw_linear(self,ax,colors):
        '''
         a method that draws linear polymers in the case that n = 1 
        '''
        numsins = 20
        ax.axis('equal')#need to be equal for proper visual scaling
        ax.set_xticks([])
        ax.set_yticks([])

        n =1 
        #make the polymers extra long, so that the section that is the right arclength can be cut off
        #this is not the most robust solution but it is quick and probably wont fail
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
        x,y = x[:i],y[:i]
        plotted = 0
        arclength = np.array(arclength)
        for myarm in self.armlist:
       	     sections = plotted+(1-np.cumsum(myarm.block_frac))*myarm.length#use the block fractions to divide into arclength sections the 1- is  because they are reversed in polyfts
       	     sections = np.flip(sections,0)
       	     myarm.species.reverse()
       	     for i in range(len(sections)):
       	         if i != len(sections)-1:
       	             insection = np.logical_and(arclength>=sections[i],arclength<sections[i+1])
       	         else:
       	             insection = np.logical_and(arclength>=sections[i],arclength<=plotted+myarm.length)
       	         self.add_line(x[insection],y[insection],color=colors[myarm.species[i]])
             plotted+= myarm.length

    def add_line(self,x,y,color):
        self.lines.append(line(x,y,color))
    def draw_lines(self,ax):
       '''
       a method to draw all the lines in the star at once for normalization
       '''
       x_max,y_max = 0,0
       for myline in self.lines:
           if max(abs(myline.x)) > x_max: xmax = max(abs(myline.x))
           if max(abs(myline.y)) > y_max: ymax = max(abs(myline.y))
       for myline in self.lines:
           ax.plot(myline.x,myline.y,color=myline.color,lw=4)
           ax.axis('square')
           #ax.set_ylim((-ymax,ymax))
           #ax.set_xlim((-xmax,xmax))
           ax.set_xticks([])
           ax.set_yticks([])



class line:
    '''
    a class that holds lines to be drawn later after normalization
    '''
    def __init__(self,x,y,color):
       self.x,self.y,self.color = x,y,color



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
    def draw_arm(self,angles,sectionsize,colors):
        ''' a function that draws arms, takes an arm as input, as well as the angle to rotate the arm
            a list of colors and the section size that is the angle that the section of the unit circle
            that the arm can take up'''
        numsins = 20

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
                    self.redraw(angles,sectionsize,colors)
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
                self.star.add_line(x[insection],y[insection],colors[self.species[i]])


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
    def redraw(self,angles,sectionsize,colors):
        self.x,self.y,self.arclength= [],[],[]
        self.draw_arm(angles,sectionsize,colors)















parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dirs', action='store',required=True, nargs='+',help='list of directories that contain each .out file to be graphed')
parser.add_argument('-p','--pickle',action='store',default='',help='the name of the filename to save a python pickle of the figure to for editing somwhere else, or sharing interactively')
parser.add_argument('--xlabel',action='store',default='f_A',help='the label for the x axis use $"label"$ for latex')
parser.add_argument('--ylabel',action='store',default='n',help='the label for the y axis use $"label"$ for latex')
args = parser.parse_args()

dirs=args.dirs
dirs.sort()
print("Dirs: {}".format(dirs))
idir=os.getcwd()
starlist = []
for mydir in dirs:
    os.chdir(mydir)
    starlist.append(star(glob.glob('*.out')[0],mydir))#create an arm for each .out file (you should only have one in each directory)
    os.chdir(idir)
colorlist = ['r','b','c','m','k','g','y']
n_arms = []
phi = []
for mystar in starlist:
    n_arms.append(mystar.n_diblock)
    phi.append(mystar.phi)
n_arms,phi = list(set(n_arms)),list(set(phi))#find unique n and phi
n_arms.sort(reverse=True)#need to resort since set removed orde
phi.sort()
fig,axlist = plt.subplots(len(n_arms),len(phi))
#make dictionaries that relate phi to a specific subplot
phi_dict = dict([(phi,i) for i, phi in enumerate(phi)])
n_dict = dict([(n_arms,i) for i,n_arms in enumerate(n_arms)])
for mystar in starlist:
    #get the subplot to graph on
    ax_index =  [n_dict[mystar.n_diblock],phi_dict[mystar.phi]]
    if len(np.shape(axlist)) == 1:#in case of a 1 dimensional plot
        ax = axlist[sum(ax_index)]
        multiplot = True
    elif len(np.shape(axlist)) == 0: #for a 0 dim plot
        ax = axlist
        multiplot=False
    else:
        ax = axlist[ax_index[0],ax_index[1]]
        multiplot=True
    
    #label edge axes
    if multiplot:
        if ax_index[1] == 0:
           ax.set_ylabel(str(mystar.n_diblock),fontsize=22)
        if ax_index[0] == len(axlist)-1:
           ax.set_xlabel(str(mystar.phi),fontsize = 22)
    else:
       ax.set_xlabel(str(mystar.phi),fontsize = 22)
       ax.set_ylabel(str(mystar.n_diblock),fontsize=22)
    #FUTURE IMPROVEMENT:
    #this is a hardcoded style for graphing miktoarms, this can be generalized for any star polymer
    #but that is not worthwile at the moment
    if mystar.n_diblock != 1:
        mystar.armlist[0].draw_arm([-pi],pi/2,colorlist)
        diblock_spacing = np.linspace(-pi/2,pi/2,mystar.armlist[1].multiplicity+2)[1:-1] #divide the starting angle into multiplicity plus two sections so endpoints can be removed
        sectionsize = pi/mystar.armlist[1].multiplicity
        mystar.armlist[1].draw_arm(diblock_spacing,sectionsize,colorlist)
    else:
        mystar.draw_linear(ax,colorlist)
    mystar.draw_lines(ax)
if args.xlabel != '':
    plt.text(0.5, 0.01, r'$\phi_A$', fontsize=30, transform=plt.gcf().transFigure)
if args.ylabel != '':
    plt.text(0.001, 0.5, r'$n$', fontsize=30, transform=plt.gcf().transFigure,rotation='vertical')
plt.tight_layout()
if args.pickle != '':
    with open(args.pickle,'wb') as pic:
        pickle.dump(fig,pic)
plt.show()

