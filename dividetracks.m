function [subtr,tracks,currStr]=dividetracks(currStr,root,count,fpath,nsubmin,div,nmin,tmin,kadd,tsigma,onlyending,ampRatioLimit)
% DIVIDETRACKS
% [subtr,tracks,currStr]=
% dividetracks(currStr,root,count,fpath,nsubmin,div,nmin,tmin,kadd,tsigma,onlyending,ampRatioLimit)
% default: 
%     root='prova';
%     count=0;
%     fpath=pwd;
%     nsubmin=+Inf;
%     div=true;
%     nmin=0;
%     tmin=0;
%     kadd=0;
%     tsigma=2.58;
%     onlyending=true;
%     ampRatioLimit=1;

if nargin<12 || isempty(ampRatioLimit)
    ampRatioLimit=1;
end
if nargin<11 || isempty(onlyending)
    onlyending=true;
end
if nargin<10 || isempty(tsigma)
    tsigma=2.58;
end
if nargin<9 || isempty(kadd)
    kadd=0;
end
if nargin<8 || isempty(tmin)
    tmin=0;
end
if nargin<6 || isempty(div)
    div=true;
end
if nargin<7 || isempty(nmin)
    nmin=0;
end
if nargin<2 || isempty(root)
    root='prova';
end
if nargin<3 || isempty(count)
    count=0;
end
if nargin<4 || isempty(fpath)
    fpath=pwd;
end
if nargin<5 || isempty(nsubmin)
    nsubmin=Inf;
end
nsubtr = size(currStr.tracksCoordAmpCG,1);
ncols=size(currStr.tracksCoordAmpCG,2)/8;
endi=find(currStr.seqOfEvents(:,2)==2);
starti=find(currStr.seqOfEvents(:,2)==1);
fmin=min(currStr.seqOfEvents(:,1));
subtr=struct('npart',num2cell(ones(nsubtr,1)),'nsubs',1,'fstart',0,'fend',0,'start',{[0;0]},'end',{[0;0]},'ampav',[]);
for c2=1:nsubtr
    currEnd=endi(find(currStr.seqOfEvents(endi,3)==c2)); %#ok<FNDSB>
    currStart=starti(find(currStr.seqOfEvents(starti,3)==c2)); %#ok<FNDSB>
    cuts=find(currStr.seqOfEvents(:,4)==c2);
    [cutframes,cutfI]=sort(currStr.seqOfEvents(cuts,1));
    framecuts=[currStr.seqOfEvents(currStart,1),cutframes',currStr.seqOfEvents(currEnd,1)];
    fcutsi=[currStart,cuts(cutfI)',currEnd];
    subtr(c2).nsubs=length(framecuts)-1;
    
    subtr(c2).joinwith=NaN(2,subtr(c2).nsubs);
    subtr(c2).splitwith=NaN(2,subtr(c2).nsubs);

    changes=[0, 2.*currStr.seqOfEvents(cuts(cutfI),2)'-3];
    
    subtr(c2).joinwith(1,end)=currStr.seqOfEvents(currEnd,4);
    subtr(c2).splitwith(1,1)=currStr.seqOfEvents(currStart,4);
    
    indJoin=find(currStr.seqOfEvents(cuts(cutfI),2)==2);
    indSplit=find(currStr.seqOfEvents(cuts(cutfI),2)==1);
    
    subtr(c2).joinwith(1,indJoin)=currStr.seqOfEvents(cuts(indJoin),3);
    subtr(c2).splitwith(1,indSplit+1)=currStr.seqOfEvents(cuts(indSplit),3);
    
    isless=currStr.seqOfEvents(fcutsi,3)'<c2;
    isless(end)=[];
    islessmn=isless & changes==-1;
    islesspiu=isless & changes==1;
    chn=ones(1,subtr(c2).nsubs);
    if any(islessmn(:))
        chn(islessmn)=arrayfun(@(s) s.npart(1),subtr(currStr.seqOfEvents(fcutsi(islessmn),3)));
    end;
    if any(islesspiu(:))
        chn(islesspiu)=arrayfun(@(s) s.npart(end),subtr(currStr.seqOfEvents(fcutsi(islesspiu),3)));
    end;
    changes=chn.*changes;
    subtr(c2).npart=cumsum(changes);
    subtr(c2).npart=1+subtr(c2).npart-min(subtr(c2).npart);

    subtr(c2).fstart=framecuts(1:end-1);
    subtr(c2).fend=framecuts(2:end)-1;
    subtr(c2).fend(end)=subtr(c2).fend(end)+1;

    subtr(c2).start=zeros(2,subtr(c2).nsubs);
    subtr(c2).end=zeros(2,subtr(c2).nsubs);

    subtr(c2).start(1,currStr.seqOfEvents(fcutsi(1:end-1),2)==2)=currStr.seqOfEvents(fcutsi(currStr.seqOfEvents(fcutsi(1:end-1),2)==2),3);
    subtr(c2).end(1,currStr.seqOfEvents(fcutsi(2:end),2)==1)=currStr.seqOfEvents(fcutsi(1+find(currStr.seqOfEvents(fcutsi(2:end),2)==1)),3);
    subtr(c2).start(1,1)=currStr.seqOfEvents(currStart,4);
    subtr(c2).end(1,end)=currStr.seqOfEvents(currEnd,4);

    samestart=find(subtr(c2).start(1,:)==0);
    if ~isempty(samestart)
        subtr(c2).start(1,samestart)=c2;
        subtr(c2).start(2,samestart)=samestart-1;
    end
    sameend=find(subtr(c2).end(1,:)==0);

    if ~isempty(sameend)
        subtr(c2).end(1,sameend)=c2;
        subtr(c2).end(2,sameend)=sameend+1;
    end
    cf=subtr(c2).start(1,1);
    if isnan(cf)
        subtr(c2).start(2,1)=NaN;
    elseif cf<c2
        subtr(c2).start(2,1)=find(subtr(cf).end(1,:)==c2,1,'first');
    else
        subtr(c2).start(2,1)=-1;
    end
    gt=subtr(c2).end(1,end);
    if isnan(gt)
        subtr(c2).end(2,end)=NaN;
    elseif gt<c2
        subtr(c2).end(2,end)=find(subtr(gt).start(1,:)==c2,1,'last');
        subtr(gt).start(2,subtr(c2).end(2,end))=size(subtr(c2).end,2);
    else
        subtr(c2).end(2,end)=-1;
    end
    startdecide=find(subtr(c2).start(2,:)==0 & subtr(c2).start(1,:)<c2);
    for cc=startdecide
        cf=subtr(c2).start(1,cc);
        subtr(c2).start(2,cc)=size(subtr(cf).end,2);
        subtr(cf).end(2,end)=cc;
    end
    enddecide=(subtr(c2).end(2,:)==0);
    subtr(c2).end(2,enddecide)=1;
    enddecide= find(enddecide & subtr(c2).end(1,:)<c2);
    for cc=enddecide
        gt=subtr(c2).end(1,cc);
        subtr(c2).end(2,cc)=1;
        subtr(gt).start(2,end)=cc;
    end

    frames=(currStr.tracksFeatIndxCG(c2,:)~=0);
    framesi=find(frames);

    subtr(c2).fname = [root 't' int2str(count) 'r' int2str(c2)];
    subtr(c2).pos = [currStr.tracksCoordAmpCG(c2,1:8:end)' currStr.tracksCoordAmpCG(c2,2:8:end)'];
    subtr(c2).pos = subtr(c2).pos(frames,:);
    subtr(c2).dpos = [currStr.tracksCoordAmpCG(c2,5:8:end)' currStr.tracksCoordAmpCG(c2,6:8:end)'];
    subtr(c2).dpos = subtr(c2).dpos(frames,:);
    subtr(c2).frames = fmin-1+framesi'; 
    subtr(c2).fpath = fpath;
    subtr(c2).amp = currStr.tracksCoordAmpCG(c2,4:8:end)';
    subtr(c2).amp = subtr(c2).amp(frames);
    subtr(c2).damp = currStr.tracksCoordAmpCG(c2,8:8:end)';
    subtr(c2).damp = subtr(c2).damp(frames);
    for c3=1:subtr(c2).nsubs
        subtr(c2).ampav = [subtr(c2).ampav mean(subtr(c2).amp(subtr(c2).frames>=subtr(c2).fstart(c3)&subtr(c2).frames<=subtr(c2).fend(c3)))];
    end
    
    for fnm=fieldnames(currStr)'
        fm=fnm{1};
        switch(fm)
            case  {'tracksCoordAmpCG','seqOfEvents'}%'tracksFeatIndxCG',
            otherwise
                sz=size(currStr.(fm),2);
                if sz==ncols  && ~ischar(currStr.(fm))
                    subtr(c2).(fm) = currStr.(fm)(c2,:)';
                    subtr(c2).(fm) = subtr(c2).(fm)(frames);
                end
        end
    end
    
    next=subtr(c2).start(:,1);
    if ~isnan(next(1)) && next(1)<c2
        cs=subtr(c2).npart(1)-(subtr(next(1)).npart(next(2))-subtr(next(1)).npart(next(2)+1));
        if cs>0
            for dummy=1:cs
                while(~isnan(next(1)) && next(1)<c2)
                    subtr(next(1)).npart(next(2))=subtr(next(1)).npart(next(2))+1;
                    next1=subtr(next(1)).start(:,next(2));
                    if next1(1)>c2 || next(2)==1 || next(1)==subtr(next(1)).start(1,next(2))
                        next=next1;
                    elseif subtr(next(1)).ampav(next(2)-1)/subtr(next(1)).npart(next(2)-1)>subtr(next1(1)).ampav(next1(2))/subtr(next1(1)).npart(next1(2))
                        next=[next(1);next(2)-1];
                    else
                        next=next1;
                    end
                end
            end
        elseif cs<0
            next2=next;
            for dummy=(1:-cs)
                next=subtr(next2(1)).end(:,next2(2)); %next=[c2;1];
                while(~isnan(next(1)) && next(1)<c2)
                    subtr(next(1)).npart(next(2))=subtr(next(1)).npart(next(2))+1;
                    next1=subtr(next(1)).end(:,next(2));
                    if next1(1)>c2 || next(2)==length(subtr(next(1)).npart) || next(1)==subtr(next(1)).end(1,next(2))
                        next=next1;
                    elseif subtr(next(1)).ampav(next(2)+1)/subtr(next(1)).npart(next(2)+1)>subtr(next1(1)).ampav(next1(2))/subtr(next1(1)).npart(next1(2))
                        next=[next(1);next(2)+1];
                    else
                        next=next1;
                    end
                end
            end
        end
    end

    next=subtr(c2).end(:,end);
    if ~isnan(next(1)) && next(1)<c2
        ce=subtr(c2).npart(end)-(subtr(next(1)).npart(next(2))-subtr(next(1)).npart(next(2)-1));
        if ce>0
            for dummy=1:ce
                while(~isnan(next(1)) && next(1)<c2)
                    subtr(next(1)).npart(next(2))=subtr(next(1)).npart(next(2))+1;
                    next1=subtr(next(1)).end(:,next(2));
                    if next1(1)>c2 || next(2)==length(subtr(next(1)).npart) || next(1)==subtr(next(1)).end(1,next(2))
                        next=next1;
                    elseif subtr(next(1)).ampav(next(2)+1)/subtr(next(1)).npart(next(2)+1)>subtr(next1(1)).ampav(next1(2))/subtr(next1(1)).npart(next1(2))
                        next=[next(1);next(2)+1];
                    else
                        next=next1;
                    end
                end
            end
        elseif ce<0
            next2=next;
            for dummy=(1:-cs)
                next=subtr(next2(1)).start(:,next2(2)); %next=[c2;size(...)];
                while(~isnan(next(1)) && next(1)<c2)
                    subtr(next(1)).npart(next(2))=subtr(next(1)).npart(next(2))+1;
                    next1=subtr(next(1)).start(:,next(2));
                    if next1(1)>c2 || next(2)==1 || next(1)==subtr(next(1)).start(1,next(2))
                        next=next1;
                    elseif subtr(next(1)).ampav(next(2)-1)/subtr(next(1)).npart(next(2)-1)>subtr(next1(1)).ampav(next1(2))/subtr(next1(1)).npart(next1(2))
                        next=[next(1);next(2)-1];
                    else
                        next=next1;
                    end
                end
            end
        end
    end
end
npartmax=0;
for c2=1:nsubtr
    npartmax=max(npartmax,max(subtr(c2).npart));
    nonansplit=~isnan(subtr(c2).splitwith(1,:));
    nonansplit(1)=false;
    nonanjoin=~isnan(subtr(c2).joinwith(1,:));
    nonanjoin(end)=false;
    subtr(c2).joinwith(2,end)=subtr(c2).end(2,end)-1;
    subtr(c2).splitwith(2,1)=subtr(c2).start(2,1)+1;
    subtr(c2).splitwith(2,nonansplit)=1;
    if any(nonanjoin)
        subtr(c2).joinwith(2,nonanjoin)=[subtr(subtr(c2).joinwith(1,nonanjoin)).nsubs];
    end
end
%%
tracks(sum([subtr.nsubs]))=struct();
k=0;

for c2=1:nsubtr
    subtr(c2).trackns=[];
    for c3=1:subtr(c2).nsubs
        frames=subtr(c2).frames>=subtr(c2).fstart(c3)&subtr(c2).frames<=subtr(c2).fend(c3);
        m=sum(frames);
        framestemp = subtr(c2).frames(frames);
        dt=framestemp(end)-framestemp(1); %***
        if (m>nmin && dt>=tmin) || div
            k=k+1;
            subtr(c2).trackns=[subtr(c2).trackns,k];
            tracks(k).frames = framestemp;
            tracks(k).fname = [subtr(c2).fname 's' int2str(c3) 'n' int2str(subtr(c2).npart(c3))];
            tracks(k).pos = subtr(c2).pos(frames,:);
            tracks(k).dpos = subtr(c2).dpos(frames,:);
            tracks(k).fpath = fpath;
            tracks(k).amp = subtr(c2).amp(frames);
            tracks(k).damp = subtr(c2).damp(frames);
            tracks(k).ampav = mean(tracks(k).amp);
            tracks(k).dampav=std(tracks(k).amp);
            tracks(k).nsubs = subtr(c2).nsubs;
            tracks(k).npart=subtr(c2).npart(c3);
            tracks(k).self=[c2; c3];
            tracks(k).start=subtr(c2).start(:,c3);
            tracks(k).end=subtr(c2).end(:,c3);
            tracks(k).npartmax=npartmax;
            tracks(k).dt=dt;
            tracks(k).n=m;
            tracks(k).start2t=[];
            tracks(k).end2t=[];
            if c3<subtr(c2).nsubs
                tracks(k).dt=dt+subtr(c2).fstart(c3+1)-subtr(c2).fend(c3);
            end
            
            for fnm=fieldnames(currStr)'
                fm=fnm{1};
                switch(fm)
                    case  {'tracksCoordAmpCG','seqOfEvents'}%'tracksFeatIndxCG',
                    otherwise
                        sz=size(currStr.(fm),2);
                        if sz==ncols && isfield(subtr(c2),fm) && ~ischar(subtr(c2).(fm))
                            tracks(k).(fm) = subtr(c2).(fm)(frames);
                        end
                end
            end
            
        else
            subtr(c2).trackns=[subtr(c2).trackns,0];
        end
    end
end;
ntracks=k;
tracks(ntracks+1:end)=[];
for k=1:ntracks
    c2=tracks(k).self(1);
    c3=tracks(k).self(2);
    if isnan(tracks(k).start(1))
        tracks(k).startt=NaN;
    else
        tracks(k).startt=subtr(tracks(k).start(1)).trackns(tracks(k).start(2));
    end
    if isnan(tracks(k).end(1))
        tracks(k).endt=NaN;
    else
        tracks(k).endt=subtr(tracks(k).end(1)).trackns(tracks(k).end(2));
    end
    if isnan(subtr(c2).joinwith(1,c3))
        tracks(k).joinwt=NaN;
    else
        tracks(k).joinwt=subtr(subtr(c2).joinwith(1,c3)).trackns(subtr(c2).joinwith(2,c3));
    end
    if isnan(subtr(c2).splitwith(1,c3))
        tracks(k).splitwt=NaN;
    else
        tracks(k).splitwt=subtr(subtr(c2).splitwith(1,c3)).trackns(subtr(c2).splitwith(2,c3));
    end
    if tracks(k).self(2)~=1 && isnan(tracks(k).splitwt)
        tracks(k).start2t=subtr(tracks(k).self(1)).trackns(tracks(k).self(2)-1);
    end
    if tracks(k).self(2)~=subtr(tracks(k).self(1)).nsubs && isnan(tracks(k).joinwt)
        tracks(k).end2t=subtr(tracks(k).self(1)).trackns(tracks(k).self(2)+1);
    end
end
if div
    shorti=find([tracks.n]<= nsubmin);
    join=zeros(1,2);
    rmtracks=double.empty(2,0);
    for k=shorti
        if ~isnan(tracks(k).splitwt) && ~isnan(tracks(k).start(1)) && (~onlyending || isnan(tracks(k).end(1)))
            with=tracks(k).splitwt; 
            res=tracks(k).startt;
            subtrleft=0;
            if isfinite(tracks(with).dampav)
                subtrleft=subtrleft+tsigma*tracks(with).dampav; %/sqrt(track(swith).n)
            end
            if isfinite(tracks(k).dampav)
                subtrleft=subtrleft+tsigma*tracks(k).dampav; %/sqrt(track(k).n); not done and not summed "in quadrature" for (over?)compensating bleaching impact
            end
            if isfinite(tracks(res).dampav)
                addright=tsigma*tracks(res).dampav; %/sqrt(track(res).n)
            else
                addright=0;
            end
            if (ampRatioLimit*(tracks(k).ampav+tracks(with).ampav-subtrleft) > ...
                    tracks(res).ampav+addright)
                corr=find(join(:,1)==res); 
                if length(corr)>1 %should never happen; consider to comment out
                    warning('Luos:dividetracks',['subtrack already joined with ' int2str(length(corr)) ' other one(s)?!?']);
                end
                if any(join(:,1)==res)
                    f=((tracks(k).ampav-tracks(res).ampav)^2+tracks(k).dampav^2)/ ...
                        ((tracks(with).ampav-tracks(res).ampav)^2+tracks(with).dampav^2);
                    if (fcdf(f,tracks(k).n,tracks(with).n)<0.5)
                        join(corr,2)=k;
                        rmtracks(1,corr)=with;
                    end
                else
                    rmtracks=[rmtracks,[k;1]]; %#ok<AGROW>
                    join=[join;res,with];  %#ok<AGROW>
                end
            end
        end
        if  ~isnan(tracks(k).joinwt) && ~isnan(tracks(k).end(1)) && (~onlyending || isnan(tracks(k).start(1)))
            with=tracks(k).joinwt;
            res=tracks(k).endt;
            subtrleft=0;
            if isfinite(tracks(with).dampav)
                subtrleft=subtrleft+tsigma*tracks(with).dampav; %/sqrt(track(swith).n)
            end
            if isfinite(tracks(k).dampav)
                subtrleft=subtrleft+tsigma*tracks(k).dampav; %/sqrt(track(k).n); not done and not summed "in quadrature" for (over?)compensating bleaching impact
            end
            if isfinite(tracks(res).dampav)
                addright=tsigma*tracks(res).dampav; %/sqrt(track(res).n)
            else
                addright=0;
            end
            if (ampRatioLimit*(tracks(k).ampav+tracks(with).ampav-subtrleft) > ...
                    tracks(res).ampav+addright)
                corr=find(join(:,2)==res); 
                if length(corr)>1 %should never happen; consider to comment out
                    warning('Luos:dividetracks',['subtrack already joined with ' int2str(length(corr)) ' other one(s)?!?']);
                end
                if any(join(:,2)==res)
                    f=((tracks(k).ampav-tracks(res).ampav)^2+tracks(k).dampav^2)/ ...
                        ((tracks(with).ampav-tracks(res).ampav)^2+tracks(with).dampav^2);
                    if (fcdf(f,tracks(k).n,tracks(with).n)<0.5)
                        join(corr,1)=k;
                        rmtracks(1,corr)=with;
                    end
                else
                    rmtracks=[rmtracks,[k;2]]; %#ok<AGROW>
                    join=[join;with,res]; %#ok<AGROW>
                end
            end
        end
    end
    rmnd=rmtracks(1,(rmtracks(2,:)==2));
    rmbgin=rmtracks(1,(rmtracks(2,:)==1));
    for rmend=rmnd
        for fldnm={'end','endt','joinwt'}
            tracks(rmend).(fldnm{1})=NaN;
        end
    end
    for rmbegin=rmbgin
        for fldnm={'start','splitwt','startt'}
            tracks(rmbegin).(fldnm{1})=NaN;
        end
    end
    join(1,:)=[];
    currStr=convtracks(tracks,join);
    subtr=[];
%     tracksold=tracks % debug; to be removed
    tracks=[];
    kmloc=0;
    for ngroup=1:length(currStr)
        [subtr1,tracks1]=dividetracks(currStr(ngroup),root,count,fpath,nsubmin,0,nmin,tmin,kmloc);
        for temp=1:length(tracks1)
            tracks1(temp).group=ngroup;
            tracks1(temp).fname=[tracks1(temp).fname,'_',int2str(ngroup)];
        end
        for temp=1:length(subtr1)
             subtr1(temp).fname=[subtr1(temp).fname,'_',int2str(ngroup)];
        end
           
        kmloc=kmloc+length(tracks1);
        subtr=catstruct(subtr,subtr1,1);
        tracks=catstruct(tracks,tracks1,2);
    end
end
if ntracks
    for trc=1:length(tracks)
        for fldnm={'start2t','end2t','startt','endt','joinwt','splitwt'}
            if ~isempty (tracks(trc).(fldnm{1})) && isfinite(tracks(trc).(fldnm{1})) && tracks(trc).(fldnm{1})
                tracks(trc).(fldnm{1})=tracks(trc).(fldnm{1})+kadd;
            end
        end
    end
else
    tracks=struct([]);
end
%'self','start','end',