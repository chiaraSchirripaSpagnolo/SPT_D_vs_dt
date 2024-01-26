# -*- coding: utf-8 -*-
"""
Last Modified 25 Jan 2024

@author: Chiara
"""

#from numpy import zeros, array
#from math import sqrt
from random import gauss
from math import sqrt,exp
import math
from numpy import array, zeros
from random import uniform, random
import sys

class Spots():
    #Class attributes:
    from inputs import (numSpecies, D_species,
                        Tbleach_species,Ton_species,Toff_species)
    NUM_SPECIES=numSpecies #species number
    D_SPECIES=D_species #diffusion coefficient for each species
    TBLEACH_SPECIES=Tbleach_species
    TON_SPECIES=Ton_species
    TOFF_SPECIES=Toff_species

    def __init__(self,stoich,idNum, pos, D, time, tstep, light,spotIdx, motion):
        self.stoich=stoich #stechiometry, number of molecules kind A,B etc.. Array
        self.idNum=idNum # identification number of the monomer. List of list, a list for each species
        self.pos=pos #position. List of 2 lists: x and y
        self.D=D #diffusion coefficient
        self.time=time #list (time in seconds)
        self.tstep=tstep #list
        self.light=light #List of list and the element are boolean.
        self.spotIdx=spotIdx
        self.motion=motion #kind of motion: brownian implemented
        
    def stepMotion (self, step, dt, monomersCount, spotsCount, sizeX, sizeY, boundary):
      '''step of motion of the spot self
      boundary may be 'periodic', 'reflective' or None.
      With periodic conditions if a spot goes out surface, it dies and another borns.
      dt in second, diff coeff in um^2/s'''
     
      if boundary not in ['periodic', 'reflective' , None]:
          print('boundary must be ''periodic'', ''reflective'' or None!')
          sys.exit()   
    
      if self.motion=='brownian':
          
          t=dt
          disp=[gauss(0, sqrt(2*self.D*t)),gauss(0, sqrt(2*self.D*t))]    
          currentPosX=self.pos[0][-1]
          currentPosY=self.pos[1][-1]
          newPos=[currentPosX+disp[0], currentPosY+disp[1]]                        
                                  
      #periodic boundary condition.
      born=[];dead=[];
      if boundary=='periodic':
          posOk=True
          if newPos[0] > sizeX:
              newPos[0]-=sizeX
              posOk=False
          elif newPos[0] < 0:
              newPos[0]+=sizeX
              posOk=False
          if newPos[1] > sizeY:
               newPos[1]-=sizeY
               posOk=False
          elif newPos[1] < 0:
              newPos[1]+=sizeY
              posOk=False
          if posOk: #spot hasn't gone out surface  
                  self.pos[0].append(newPos[0])
                  self.pos[1].append(newPos[1])
                  self.time.append(step*dt)
                  self.tstep.append(step)
          else:#self dies and another spot borns
              numSpecies=len(self.stoich)
              idN=[[]]*numSpecies
              dead=self
              monomersCount+=self.stoich 
              spotsCount+=1
              for j in range (numSpecies):
                  if self.stoich[j]!=0:
                      idN[j]=list(range(monomersCount[j]-self.stoich[j]+1, monomersCount[j]+1))

              born=Spots(self.stoich,idN, [[newPos[0]], [newPos[1]]], self.D, [step*dt], [step], self.light, spotsCount, self.motion)
                  
      elif boundary=='reflective':
          if newPos[0] > sizeX:
              newPos[0]=sizeX-(newPos[0]-sizeX)
          elif newPos[0] < 0:
              newPos[0]=-newPos[0]
          if newPos[1] > sizeY:
               newPos[1]=sizeY-(newPos[1]-sizeY)
          elif newPos[1] < 0:
              newPos[1]=-newPos[1]
          self.pos[0].append(newPos[0])
          self.pos[1].append(newPos[1])
          self.time.append(step*dt)
          self.tstep.append(step)
                   
      elif boundary==None:
          self.pos[0].append(newPos[0])
          self.pos[1].append(newPos[1])
          self.time.append(step*dt)
          self.tstep.append(step)
                
      return  born, dead
    
    
    def BlinkBleach (self, Tbleach, Ton, Toff,dt):
        ''' Tbleach, Ton, Toff are lists with one element for each channel.
        True=on. False = off (temporarily), None=bleached'''
        st=self.stoich
        l = [ [] for i in range (len(st)) ]       
        self.light.append(l)
        
        for ch, Nmon in enumerate(st):
            if Nmon!=0: #there are some monomers in channel ch
                for m in range(Nmon): #for each monomer in channel ch
                   
                    if self.light[-2][ch][m]==None: #monomer m in ch was already bleached
                        self.light[-1][ch].append(None)
                        
                    else: # mon m in ch was not bleached
                        #check bleaching
                        prob=math.exp(-dt/Tbleach[ch])
                        rnd=random()
                        if rnd > prob:#bleach
                            self.light[-1][ch].append(None)    
                        else: #not bleach. check blinking
                        #read if it's on/off to choose tau
    
                            lgth=self.light[-2][ch][m]

                            if lgth == True:
                                    tau=Ton[ch]
                            elif lgth == False:
                                    tau=Toff[ch]
                            #calculate probability to change the state using the chosen tauOn or tauOff 
                            prob=math.exp(-dt/tau)
                            rnd=random()
                            if rnd > prob: #change the state
                                self.light[-1][ch].append(not lgth)
                            else:
                                self.light[-1][ch].append(lgth)
                    
    @staticmethod    
    def CreatePopulation  (numSpecies, InitialSpots, sizeX, sizeY, t, tstep):
        '''create a population of spots with stechiometries specified in InitialSpots
        in a surface of dimensions sizeX and sizeY, at the time t'''
        
        monomersCount=zeros(numSpecies, dtype='int') 
        population=[]
        spotsCount=0 #total spots count
        for i in InitialSpots: #i=stoich
            idN=[[]]*numSpecies
            st=array(i,dtype='int') #steich
            n=InitialSpots[i]
            Dcoeff=Spots.Diff_Coeff_function(st)
            lght=[]
            for j in range (numSpecies):
                l=[ True for i in range(st[j]) ]
                lght.append(l)
            for m in range(n):
                monomersCount+=st
                spotsCount+=1
                for j in range (numSpecies):
                      if st[j]!=0:
                          idN[j]=list(range(monomersCount[j]-st[j]+1, monomersCount[j]+1))

                coord=[[uniform(0, sizeX)],[uniform(0, sizeY)]]
                population.append(Spots(st,idN[:],coord,Dcoeff, [t], [tstep],[lght],spotsCount, 'brownian'))
        return population, monomersCount, spotsCount
                
    @staticmethod
    def Diff_Coeff_function(stoich):
        '''For possible extensibilites. calculate the diffusion coefficient for the spots made of more than one monomer.
        It use the stoichiometry and the diffusion coefficient of each species defined in 
        the class attribute Spots.D_SPECIES'''
        if Spots.D_SPECIES.dtype != 'float64': 
             Spots.D_SPECIES.astype(float) 
        D=(sum(stoich*Spots.D_SPECIES**(-2)))**(-1/2) 
        return D            
                
  
          