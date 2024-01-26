# -*- coding: utf-8 -*-
"""
Last Modified 25 Jan 2024

@author: Chiara
"""

def AllTracks (Spots, nchannels=1):

    import TracksClass as Tracks
    
    AllchTracks=[ []for i in range (nchannels) ]
    
    for spot in Spots:
        for ch in range (nchannels):
            if spot.stoich[ch]!=0:
                
                startFrame=spot.tstep[0]
                endFrame=spot.tstep[-1]
                events=[[],[]]
                events[0]=([startFrame, 1, 0, float ('nan') ])
                events[1]=([endFrame, 2, 0, float ('nan') ])
        
                idxInMovie=[]
                
                #coord and stoichiometry info
                coord=[[],[]];
                x=spot.pos[0]; y=spot.pos[1]
                nFrames=len(spot.tstep)
                stoich=nFrames*[spot.stoich]
                coord[0]=x; coord[1]=y
                
                track=Tracks.Tracks(idxInMovie, coord, stoich, events )
                AllchTracks[ch].append(track)
                del track 
    return AllchTracks