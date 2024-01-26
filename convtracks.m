function currStr=convtracks(tracks,join)
% ntracks=length(tracks);
% trn=1:ntracks;
if nargin<2 || isempty(join)
    join=0;
end
[~, trn]=sort(arrayfun(@(t) t.frames(1),tracks));

first=min(cat(1,tracks.frames));
last=max(cat(1,tracks.frames));
lngth=last-first+1;
FICG=zeros(1,lngth);
seqOE=double.empty(0,4);
    function trntojoin=jointr(cl)
        jil=find(join(:,1)==cl);
        if length(jil)>1 %should never happen; consider to comment out
            warning('Luos:dividetracks',['subtrack to be joined with ' int2str(length(jil)) ' other one(s)']);
            jil=jil(1);
        end
        if isempty(jil)
            trntojoin=[];
        else
            trntojoin=[join(jil,2),jointr(join(jil,2))];
        end
    end
linecc=0;
end2t=double.empty(2,0);
while(~isempty(trn))
    list=[trn(1),jointr(trn(1))];
    trn(ismember(trn,list))=[];
    trc=tracks(list(1));
    if ~isempty(tracks(list(end)).end2t)
        end2t=[end2t,[list(end);tracks(list(end)).end2t]];
    end
    if isempty(trc.start2t) && ~any(end2t(2,:)==list(1))
        linecc=linecc+1;
        linec=linecc;
        convtrnum{linec}=list; 
        CACG(linec,:)=NaN(1,8*lngth);
        seqOE(2*linec-1,:)=[trc.frames(1),1,linec,trc.startt];
        seqOE(2*linec,:)=[tracks(list(end)).frames(end),2,linec,tracks(list(end)).endt];
    else
        if isempty(trc.start2t)
            temp=find(end2t(2,:)==list(1));
            tofind=end2t(1,temp);
            end2t(:,temp)=[];
            linec=find(cellfun(@(c) any(c==tofind),convtrnum));
        else
            linec=find(cellfun(@(c) any(c==trc.start2t),convtrnum));
        end
        seqOE(2*linec,4)=tracks(list(end)).endt;
        seqOE(2*linec,1)=tracks(list(end)).frames(end);
        convtrnum{linec}=[convtrnum{linec},list]; 
    end
    for trnl=list
        trc=tracks(trnl);
        frames=trc.frames-first;
        FICG(linec,frames+1)=trnl;
        
        CACG(linec,1+8*frames)=trc.pos(:,1)'; %#ok<*AGROW>
        CACG(linec,2+8*frames)=trc.pos(:,2)';
        CACG(linec,5+8*frames)=trc.dpos(:,1)';
        CACG(linec,6+8*frames)=trc.dpos(:,2)';
        CACG(linec,4+8*frames)=trc.amp';
        CACG(linec,8+8*frames)=trc.damp';
        CACG(linec,3+8*frames)=0;
        CACG(linec,7+8*frames)=0;
    end
end
linecl= find(~isnan(seqOE(:,4)));
for llc=linecl'
    seqOE(llc,4)=find(cellfun(@(c) any(c==seqOE(llc,4)),convtrnum));
end

seqOE((seqOE(:,2)==2&~isnan(seqOE(:,4))),1)=seqOE((seqOE(:,2)==2&~isnan(seqOE(:,4))),1)+1;

groups=0;
searchin=seqOE(linecl,3:4);
searchn=true(1,linecc);
while any(searchn)
    strt=find(searchn,1);
    groups=groups+1;
    [groupgroups,searchin]=findgroup(strt,searchin);
    if ~all(searchn(groupgroups)) %for debug; to be commented out
        error('Luos:convtracks','Something wrong in convtracks')
    end
    searchn(groupgroups)=false;
    group{groups}=groupgroups;
end
if groups==1
    currStr.tracksFeatIndxCG=FICG;
    currStr.tracksCoordAmpCG=CACG;
    currStr.seqOfEvents=sortrows(seqOE);
else
    for ii=1:groups
        currStr(ii).seqOfEvents=seqOE(union(2*group{ii},2*group{ii}-1),:);
        minloc=min(currStr(ii).seqOfEvents(:,1))-first+1;
        maxloc=max(currStr(ii).seqOfEvents(:,1))-first+1;
        currStr(ii).tracksFeatIndxCG=FICG(group{ii},minloc:maxloc);
        currStr(ii).tracksCoordAmpCG=CACG(group{ii},(8*minloc-7):8*maxloc);
        lc=[]; %for debugging; could be commented out.
        lc(group{ii})=1:length(group{ii});
        currStr(ii).seqOfEvents(isfinite(currStr(ii).seqOfEvents(:,4)),4)=...
            lc(currStr(ii).seqOfEvents(isfinite(currStr(ii).seqOfEvents(:,4)),4));
        currStr(ii).seqOfEvents(:,3)=lc(currStr(ii).seqOfEvents(:,3));
        currStr(ii).seqOfEvents=sortrows(currStr(ii).seqOfEvents);
    end
    currStr=reshape(currStr,[],1);
end
end

function [group,searchout]=findgroup(st1,searchin)
    group=[];
    searchout=[];
    if nargin==2 && ~isempty(searchin)
        [r,~]=find(searchin==st1);
        group=unique(searchin(r,:));
        searchin(r,:)=[];
        searchout=searchin;
        group=setdiff(group,st1);
        group=reshape(group,1,[]);
        if ~isempty(group)
            for st2=group
                [group2,searchout]=findgroup(st2,searchout);
                group=union(group,group2);
            end
        end
    end
    group=unique([st1,group]);
end
