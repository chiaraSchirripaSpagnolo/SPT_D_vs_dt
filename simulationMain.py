# -*- coding: utf-8 -*-
"""
Last Modified 25 Jan 2024 

@author: Chiara
"""

import os
import scipy.io as sio
import math
import simulation as sim
import toMatLabNew as ToMat
import pickle
import TracksClass #it is used
import ExtractTracks_New as tr
import TracksToMatlab

n=input('Number of simulations?') #num of sim with the same parameters

for simIdx in range(int(n)):        
    res=sim.Sim_Diffusion()    
    if simIdx==0:
        Density=res[1]
        parameters=res[2]
        parameters['density']=Density
        InitialSpots=res[3]

        FolderName=''.join(( 'Dens',"{:.2f}".format(Density[0]), '_',''.join(('dt', str(parameters['dt']))), 
                             '_',''.join(('t', str(parameters['t']))),'_', 'D',str(parameters['D'][0]),
                            ))
        
        os.mkdir(FolderName)
        
        #change curr dir to FolderName directory
        os.chdir(FolderName)
        
        sio.savemat('InitialSpots.mat', {'InitialSpots':InitialSpots})
        
    All=res[0]
    
    #create a folder for each simulation
    subfoldName='_'.join(('sim', str(simIdx)))
    os.mkdir(subfoldName)
    filename = 'res'
    path=os.path.join(subfoldName,filename)
    outfile = open(path,'wb')
    pickle.dump(res,outfile)
    outfile.close()
    
    spotsXframe, spotsXframeON, idx_ON_all = ToMat.SaveSpotsXFrame (All, subfoldName, parameters['dt'], parameters['t'],parameters['blinkbleach'], nchannels=parameters['numSpecies'])
    
    #extract and save tracks
    tracks=tr.AllTracks(All,parameters['numSpecies'])
    TracksToMatlab.saveTracks(parameters['numSpecies'],tracks,subfoldName)
    
sio.savemat('parameters.mat', {'parameters':parameters})