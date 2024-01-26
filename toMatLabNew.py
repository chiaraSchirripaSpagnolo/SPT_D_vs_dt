# -*- coding: utf-8 -*-
"""
Last Modified 25 Jan 2024

@author: Chiara
"""

# -*- coding: utf-8 -*-
"""
extract spots' positions for each frame;
save for matlab

"""
def SaveSpotsXFrame (All, dir2save, dt, t, blinkbleach,nchannels=1):
    
    import numpy as np
    import scipy.io as sio
    import os
    
    nf=int(t/dt)+1
    
    channels=[];
    spotsXframe=[np.empty((nf,), dtype=object) for i in range(nchannels)]
    spotsXframeON=[np.empty((nf,), dtype=object) for i in range(nchannels)]
    idx_ON_all=[np.empty((nf,), dtype=object) for i in range(nchannels)]
    
    for ch in range(1,nchannels+1):
        channels.append( 'channel_{:d}'.format( ch))
    
    for i in All: #loop on the spots
        st=i.stoich
        
        #inizialize i.idxInMovie (empty) with the correct shape, list of:
        #one element for each frame, in each element  a list for each channel
        l = [ [] for i in range (nchannels) ]
        spotLife = len(i.tstep)
        i.idxInMovie = [ l[:] for i in range(spotLife)]
        
        #loop on the spot's track
        for frIdx, frame in enumerate (i.tstep): #tstep starts from 0
            #extract positions in this frame
            x=i.pos[0][frIdx]
            y=i.pos[1][frIdx]
            for ch, n in enumerate (st):
                if blinkbleach:
                    lgth=i.light[frIdx][ch]
                    m=sum(bool(x) for x in lgth) #numbers of monomers with light On
                
                if spotsXframe[ch][frame] is None: #this frame is still empty
                    spotsXframe[ch][frame]=[n*[x],n*[y]]
                    i.idxInMovie[frIdx][ch]=list(range(n))    
                else:
                    spotsXframe[ch][frame][0].extend(n*[x]) 
                    spotsXframe[ch][frame][1].extend(n*[y])
                    i.idxInMovie[frIdx][ch]=list(range(len(spotsXframe[ch][frame][0])-n,len(spotsXframe[ch][frame][0])))
                if blinkbleach and m!= 0:
                    if spotsXframeON[ch][frame] is None: #this frame is still empty
                        spotsXframeON[ch][frame]=[m*[x],m*[y]]
                        if idx_ON_all[ch][frame] is None:
                            idx_ON_all[ch][frame]=list(range(len(spotsXframe[ch][frame][0])-m,len(spotsXframe[ch][frame][0])))
                        else:
                            idx_ON_all[ch][frame].extend(list(range(len(spotsXframe[ch][frame][0])-m,len(spotsXframe[ch][frame][0]))))
            
                    else:
                        spotsXframeON[ch][frame][0].extend(m*[x]) 
                        spotsXframeON[ch][frame][1].extend(m*[y])
                       
                        #write variable that for each on spot store the index in the full pool of spots. Index from 0.
                        if idx_ON_all[ch][frame] is None:
                            idx_ON_all[ch][frame]=list(range(len(spotsXframe[ch][frame][0])-m,len(spotsXframe[ch][frame][0])))
                        else:
                            idx_ON_all[ch][frame].extend(list(range(len(spotsXframe[ch][frame][0])-m,len(spotsXframe[ch][frame][0]))))
                   
    for i in spotsXframe:
        i[i == np.array(None)]=float('nan')        
        #None not readable by Matlab, NaN readable by Matlab
    for i in spotsXframeON:
        i[i == np.array(None)]=float('nan')
    for i in idx_ON_all:
        i[i == np.array(None)]=float('nan')
        
    DictSpXfr=dict( zip(channels,spotsXframe ))
    path=os.path.join(dir2save,'spotsXframe.mat')
    sio.savemat(path, DictSpXfr)
        
    if blinkbleach:
        DictSpXfrON=dict( zip(channels,spotsXframeON ))
        DictIdx=dict( zip(channels,idx_ON_all ))
        path=os.path.join(dir2save,'spotsXframeON.mat')
        sio.savemat(path, DictSpXfrON)
        path=os.path.join(dir2save,'idx_onEall.mat')
        sio.savemat(path, DictIdx)
          
    return spotsXframe, spotsXframeON, idx_ON_all