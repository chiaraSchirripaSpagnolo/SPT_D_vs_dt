function [tracksFinal] = TracksFinal_fromSim_dt(tracks, dt, nFrames)
%From the simulated tracks it extracts tracks with utrack formatting with
% temporal resolution dt

%INPUTS:
%tracks: simulated tracks
%dt: time resolution
%nFrames: number of desired frames in tracksFinal

%OUTPUTS
%tracksFinal tracks in utrack format

pix_size=0.16;
nTr=length(tracks);
tracksFinal=struct('tracksFeatIndxCG',[],'tracksCoordAmpCG',[],'seqOfEvents',[]);
%fields={'tracksFeatIndxCG','tracksCoordAmpCG','seqOfEvents'};

idxTrFinal=0;
for i=1:nTr % loop on tracks
    tr=tracks{i}; % 1 x 3 cell. actual track
    
    %seqOfEvents matrix
    events=reshape(tr{4},[],4);
    %add 1 to frame and track indexes (python indexed start from 0)
    events(:,[1,3,4])=events(:,[1,3,4])+1;
    events=sortrows(events,1);
    frames_dt=1:dt:dt*nFrames+1-dt;
    FirstFrameInTr=events(events(:,2)==1);
    LastFrameInTr=events(events(:,2)==2);
    FramesInTr=FirstFrameInTr:1:LastFrameInTr;
    FramesToConsider=intersect(frames_dt,FramesInTr);
    if isempty (FramesToConsider)
        % e.g. the track starts after the frames I want to consider
        continue
    end
    firstFrame_dt=FramesToConsider(1);
    lastFrame_dt=FramesToConsider(end);
    if isempty (firstFrame_dt)
        continue
    end

    idxTrFinal=idxTrFinal+1;
    events(events(:,2)==1)=firstFrame_dt;
    events(events(:,2)==2)=lastFrame_dt;
    tracksFinal(idxTrFinal).seqOfEvents=events;
 
    %ok also if simulated tracks do not start from the first frame
    nFramesTr=size(FramesToConsider,2);
    coordAmp=nan(1,8*nFramesTr);
    last=1+8*(nFramesTr-1);
    firstF_idx=find(FramesInTr==FramesToConsider(1));
    lastF_idx=find(FramesInTr==FramesToConsider(end));
    coordAmp(1,1:8:last)=tr{2}(1,firstF_idx:dt:lastF_idx)./pix_size; %x
    coordAmp(1,2:8:last+1)= tr{2}(2,firstF_idx:dt:lastF_idx)./pix_size;%y
    
    coordAmp(1,3:8:last+2)= 0;%z
    coordAmp(1,4:8:last+3)= 1; %A
    coordAmp(1,5:8:last+4)= 0;%dx
    coordAmp(1,6:8:last+5)= 0;%dy
    coordAmp(1,7:8:last+6)= 0;%dz
    coordAmp(1,8:8:last+7)= 0;%dA

    tracksFinal(idxTrFinal).tracksCoordAmpCG=coordAmp;
    tracksFinal(idxTrFinal).tracksFeatIndxCG=repmat(idxTrFinal,1,nFramesTr);
end
name=['TracksFinal_fromSim_dt', num2str(dt),'_',num2str(nFrames) 'frames.mat'];
save (name,'tracksFinal' )
end