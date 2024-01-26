# -*- coding: utf-8 -*-
"""
Last Modified 25 Jan 2024

@author: Chiara
"""

# -*- coding: utf-8 -*-
"""
Simulation of diffusing molecules.
It uses the settings specified in the module inputs.py
and the class defined in SpotsClass.py 

"""
def Sim_Diffusion():
    #%% 
    import numpy as np
    import time
    start1 = time.time()
    import inputs
    import SpotsClass as Sp
    
    #%%Inputs' parameters
    InitialSpots={str(k): v for k,v in inputs.InitialSpots.items()}
        
    parameters = {'dt':inputs.dt, 't':inputs.t, 'size':inputs.size,
                  'numSpecies':inputs.numSpecies, 'D':inputs.D_species,
                  'blinkbleach':inputs.blinkbleach,'Tbleach_species':inputs.Tbleach_species,
                  'Ton_species':inputs.Ton_species, 'Toff_species':inputs.Toff_species,
                  }
    
    nstep=inputs.t/inputs.dt
    #%%create initial population of spots
    pop=Sp.Spots.CreatePopulation  (inputs.numSpecies,inputs.InitialSpots,inputs.size,inputs.size,0, 0)
    Particles=pop[0]
    monomersCount=pop[1] #total monomer count
    spotsCount=pop[2] #total spot count 
    Area=inputs.size*inputs.size
    Density=monomersCount/Area
    
    #%% inzialize blink and bleach 
    if inputs.blinkbleach:
        for sp in Particles:
            if inputs.blinkbleach:
                sp.BlinkBleach(inputs.Tbleach_species, inputs.Ton_species, inputs.Toff_species, inputs.dt)
            
    #%% inizialize variables
    DeadSpots=[] #Spots gone out surface 
    BornSpotsAll=[]
    #%% cycle
    print('start cycle')
    step=1
    while step<=nstep:       
        
        #%%motion  
        BornSpots=[]#spots entering the surface if periodic boundaries (this time step)
        DeadIdx=[]#indexes of spots that came out of the surface if periodic boundaries (all time steps)
        
        for i in range(len(Particles)):
            b=[]; d=[]
            b,d=Particles[i].stepMotion (step, (inputs.dt*10**-3),monomersCount,spotsCount, inputs.size, inputs.size, inputs.boundary)
            
            if d: #if d non empty also b non empty
                DeadIdx.append(i)
                DeadSpots.append(d)
                BornSpots.append(b)   
                BornSpotsAll.append(b)
        Particles = [Particles[i] for i in range(len(Particles)) if i not in DeadIdx]
        Particles.extend(BornSpots)
                
        #%% blink and bleach
        if inputs.blinkbleach:
            for sp in Particles:
                 if inputs.blinkbleach:
                     sp.BlinkBleach(inputs.Tbleach_species, inputs.Ton_species, inputs.Toff_species, inputs.dt)
                
        step+=1
    #%%Pool all spots existed
    All=Particles+DeadSpots
    All=[ sp for sp in All if len(sp.tstep) > 0 ]
    end1 = time.time()
    print(end1 - start1)
    #%% show all trajectories on request
#    while True:
#        show=input ('Do you want to show trajectories? (y/n) ')
#        if show in ('y' ,'yes'):
#            import matplotlib.pyplot as plt
#            for i in All:
##                if isinstance (i.pos,np.ndarray ):
##                    if i.pos.ndim==1:
##                        coordx=i.pos[0]
##                        coordy=i.pos[1]
##                    else:
##                        coordx=i.pos[:,0]
##                        coordy=i.pos[:,1]
##                elif isinstance (i.pos,list ): 
#                coordx=i.pos[0]
#                coordy=i.pos[1]
#                plt.style.use('dark_background')
#                plt.plot(coordx,coordy, 'g-')
#            break
#        elif  show in ('n' ,'no'):
#             break
#        else:
#            print('Write: ''y'' or ''n''')   
 
    return All, Density, parameters, InitialSpots, DeadSpots, BornSpots, BornSpotsAll