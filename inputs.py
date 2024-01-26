# -*- coding: utf-8 -*-
"""
Last Modified 25 Jan 2024

@author: Chiara
"""

# -*- coding: utf-8 -*-
"""
Inputs definitions for simulation.py
Define parameters:
time step, total time,
surface dimensions, surface boundary,
numbers of species (for extendibility), 
diffusion coefficient,
blinking, bleaching parameters,
initial spots 
It includes some controls on the variables.

"""
import numpy as np
#from sys import exit
import sys

#%%simulation's parameters
dt=40 #time step (ms)
t=40000 #total time (ms)
size=16 #surface side in um.
#BOUNDARY conditions. Admitted: 'periodic', 'reflective' or None
boundary = 'reflective'
  
#%%Species' features
numSpecies=1 #allows extensibility to multiple species

D_species=np.array([0.3],dtype=float) #um^2/s

blinkbleach=True
Tbleach_species=[np.inf] #ms NoBleach: Ton=np.inf
Ton_species=[np.inf] #ms NoBlink: Ton=np.inf
Toff_species=[np.inf] #ms


#%%INITIAL SPOTS
#InitialSpots: types and number of initial spots
#Dictionary. key=tuple of stoichiometry, value=initial number of spots
InitialSpots= {             (1,) :110,
                    
              }

#%% controls
        
# check initial spots

if type(InitialSpots) != dict:
    print('InitialSpots must be a dictionary!')
    sys.exit()
       
#values=number of initial species must be integers:
def check_dictValueInteger (d):
    for v in d.values():
        if not isinstance(v, int):
           print('Error: Numbers of initial spots must be integers!')
           sys.exit()
check_dictValueInteger (InitialSpots)
   
# stechiometries must have length numSpecies
# and must be composed by integers and not (0,0)
def checkStoich(var,numSpecies):
    for k in var:
        if len(k) != numSpecies:
            print ('Input Error: Stoichiometries must have length equal to number of species!')
            sys.exit()
        if all( [isinstance(k[i], int) for i in range(numSpecies)] )==False:
            print('Input Error: Stoichiometries must contain integers!')
            sys.exit()
        if all ( [k[i]==0 for i in range(numSpecies)] )==True:
            print('Input Error: Stoichiometries with all zeros do not make sense!')
            sys.exit()

checkStoich(InitialSpots,numSpecies )    

