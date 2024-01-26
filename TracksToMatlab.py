# -*- coding: utf-8 -*-
"""
Last Modified 25 Jan 2024

@author: Chiara
"""

# -*- coding: utf-8 -*-
"""
From AllchTracks which contains Trakcs objects, build
a variable saved in .mat
"""
def saveTracks (nchannels,AllchTracks, dir2save):
    
    import numpy as np
    import scipy.io as sio
    import os
    
    for ch in range(nchannels):
        
        nTr=len(AllchTracks[ch])
        tracks=np.empty((nTr,), dtype=object) #nTr tracks
        for i in range(nTr):
            #tracks[i]=np.empty((3,), dtype=object) # 3 info: movieIdx, coordStoich, events
            tracks[i]=np.empty((4,), dtype=object) # 4 info: movieIdx, coord,stoich, events
            tracks[i][0]=AllchTracks[ch][i].idxInMovie
            tracks[i][1]=AllchTracks[ch][i].coord
            tracks[i][2]=AllchTracks[ch][i].stoich
            tracks[i][3]=AllchTracks[ch][i].events
            
        DictTr=dict( tracks=tracks )
        name='TracksSim'.join(['channel_{:d}_'.format( ch),'.mat'])
        path=os.path.join(dir2save,name)
        sio.savemat(path, DictTr)