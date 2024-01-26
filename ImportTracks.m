function tracks = ImportTracks(FileName,PathName,FilterIndex,nmin,tmin,fnm)
% tracks = ImportTracks(FileName,PathName,FilterIndex,nmin,tmin,fnm)
% import tracks saved in files with various formats for subsequent
% analyses.
% (c) Stefano Luin 2024 s.luin<at>sns.it
% simplified version to avoid useless file dependencies
%
% NOTE: nmin and tmin are the max number of frame and DIFFERENCE within
% last and first frame that are removed. I'm not fixing this bug for
% consistence with old analyses. Use min_n - 1 and min_t - 2 for
% results with minimum number of frames min_n and minimum duration min_t
% (last-first+1)
% optdlg.Resize='on';
% optdlg.WindowStyle='normal';
% optdlg.Interpreter='none';
%global msd_opts;
drawnow; pause(0.005);
if nargin<6
    fnm='no file';
end
if nargin<4 %||isempty(nmin)... to be finished
    %     prompt = {'nmin:','tmin:'};
    %     answer = inputdlg(prompt);
    %     nmin=str2num(answer{1}); tmin=str2num(answer{2});
    nmin=5;
    tmin=7;
end
tracks = struct([]);

if nargin<3 || isempty(FileName)
    f=fileext('traj');
    [FileName,PathName,FilterIndex] = uigetfile(f(:,1:2),'MultiSelect','on');
    FilterIndex=f{FilterIndex,3};
end

k = 0;

if isempty(FileName), return; end

if iscell(FileName)
    n=length(FileName);
else
    n = 1;
end
for i=1:n
    if iscell(FileName)
        fname = FileName{i};
    elseif isstruct(FileName)
        fname = fnm;
    else
        fname = FileName;
    end
    fpath = PathName;
    findex = FilterIndex;
    
    disp([fname ' loading... ' datestr(now)])
    
    if findex==1 % File Excel
        [workbook, ~, ~]=importdataL(fullfile(fpath,fname));
        
        if isfield(workbook,'data') && not(isfield(workbook,'Position')) && not(isfield(workbook,'Track'))
            workbook=workbook.data;
        end
        
        % Excel file output from Imaris
        if isfield(workbook,'Position')
            worksheet = workbook.Position;
            pos= worksheet(:,1:2);
            frames = worksheet(:,7);
            m = size(pos,1);
            
            amp=[];    
            if isfield(workbook, 'IntensityMeanCh0x3D1')
                worksheet = workbook.IntensityMeanCh0x3D1;
                if m==size(worksheet,1)
                    amp=worksheet(:,1);
                end
            end
            
            damp=[];
            if isfield(workbook, 'IntensityStdDevCh0x3D1')
                worksheet = workbook.IntensityStdDevCh0x3D1;
                if m==size(worksheet,1)
                    damp=worksheet(:,1);
                end
            end         
            
            if (m>nmin && (frames(end)-frames(1))>=tmin && all(diff(frames)>0)) %(frames(end)-frames(1))>20 luin28/3/2011
                k=k+1;
                tracks(k).fname = fname;
                tracks(k).pos = pos;
                tracks(k).frames = frames;
                tracks(k).fpath = fpath;
                tracks(k).amp = amp; tracks(k).damp = damp;
            end
        end
        
        % if the file has only one track
        if isfield(workbook,'Track')
            workbook=workbook.data;
            worksheet = workbook.Track;
            
            pos=worksheet(:,2:3);
            frames = worksheet(:,1);
            m = size(pos,1);
            
            amp=[];    
            if size(worksheet,2)>3
                amp=worksheet(:,4);
            end        
            damp=[];
            
            if (m>nmin && (frames(end)-frames(1))>=tmin && all(diff(frames)>0)) %(frames(end)-frames(1))>20 luin28/3/2011
                k=k+1;
                tracks(k).fname = fname;
                tracks(k).pos = pos;
                tracks(k).frames = frames;
                tracks(k).fpath = fpath;
                tracks(k).amp=amp; tracks(k).damp=damp;
            end
        end
    end
    
    if findex==2 % File CSV
        % Considering one trajectory only
        [worksheet, ~, ~]=importdataL(fullfile(fpath,fname),',',1);
        pos=worksheet.data(:,3:4);
        %       time=worksheet.data(:,2);
        frames = worksheet.data(:,1);
        m = size(pos,1);
        
        amp=[];     
        if size(worksheet.data,2)>4
            amp=worksheet.data(:,5);
        end
        damp=zeros(size(amp));
        
        dtl=mean(diff(worksheet.data(:,2)));
        global msd_opts; %#ok<*GVMIS,*TLEV>
        msd_opts.tau=dtl;
        msd_opts.dt=dtl;
        
        if (m>nmin && (frames(end)-frames(1))>=tmin && all(diff(frames)>0)) %(frames(end)-frames(1))>20 luin28/3/2011
            k=k+1;
            tracks(k).fname = fname;
            tracks(k).pos = pos;
            tracks(k).dpos=zeros(size(pos));
            tracks(k).frames = frames;
            tracks(k).fpath = fpath;
            tracks(k).amp=amp;
            tracks(k).damp=damp;
            tracks(k).dtl=dtl; %unused; just for checking
        end
    end
    
    if findex==3 || findex==5 % File DAT or TXT
        % Luin: DAT implemented as output of converter_setup.
        xdata=scanfirsts(fullfile(fpath,fname),4,1,'\t');
        pos=xdata(:,2:3);
        frames = xdata(:,1);
        m = size(pos,1);
        amp=xdata(:,4);
        damp=[];
        if (m>nmin && (frames(end)-frames(1))>=tmin && all(diff(frames)>0)) %(frames(end)-frames(1))>20 luin28/3/2011
            k=k+1;
            tracks(k).fname = fname;
            tracks(k).pos = pos;
            tracks(k).frames = frames;
            tracks(k).fpath = fpath;
            tracks(k).amp=amp;
            tracks(k).damp=damp;
        end
    end
    
    if findex==4 % File MAT from loadImaris
        fileStr=load(fullfile(fpath,fname));
        if isfield(fileStr,'trcks')
            nStr=length(fileStr.trcks);
            for count=1:nStr
                m = size(fileStr.trcks(count).pos,1);
                if (m>nmin && (fileStr.trcks(count).frames(end)-fileStr.trcks(count).frames(1))>tmin && all(diff(fileStr.trcks(count).frames)>0)) %(frames(end)-frames(1))>20 luin28/3/2011
                    k=k+1;
                    fileStr.trcks(count).fpath=fpath;
                    %                    tracks(k)=fileStr.trcks(count);
                    tracks(k).fname = [fileStr.trcks(count).fname '.mat'];
                    tracks(k).pos = double(fileStr.trcks(count).pos);
                    tracks(k).frames = double(fileStr.trcks(count).frames);
                    tracks(k).fpath = fileStr.trcks(count).fpath;
                    tracks(k).amp=double(fileStr.trcks(count).amp);
                    tracks(k).damp=double(fileStr.trcks(count).damp);
                    
                end
            end
        else
            disp([fname ' is not a valid file']);
        end
        clear fileStr; 
    end
    %% From here on: output from utrack or use utrack to determine the trajectories from movies
    exUTrack=false;
    strdt=datestr(now,'yymmddHHMMSS'); %#ok<*TNOW1,*DATST>
    if findex==10 || findex==11 % movie to be analyzed with utrack
        error('Luin:ImportTracks','Option not implemented in this simplified version of ImportTracks')
    end
    
    if findex==8 || findex==9 % File MAT output of u-track.
        if isequal(PathName,0)
            fileStr.tracksFinal=FileName;
            exUTrack=1;
        end
        if ~exUTrack, fileStr=load(fullfile(fpath,fname)); end
        if isfield(fileStr,'tracks'),trcks='tracks';
        elseif isfield(fileStr,'tracksFinal'),trcks='tracksFinal';
        else trcks=''; %#ok<*SEPEX>
        end
        if ~isempty(trcks)
            nStr=length(fileStr.(trcks));
            tracksFinal=struct([]);
            [partdir,uudir,~]=fileparts(fileparts(fileparts(fullfile(fpath,fname))));
            if strcmp(uudir,'TrackingPackage')
                [~, root, ~]=fileparts(partdir);
            else
                [~, root, ~]=fileparts(fname);
            end
            root=regexprep(root,'\W','');
            for count=1:nStr %trajectory number
                currStr=fileStr.(trcks)(count);%struct corresponding to trajecotory count
                nsubtr = size(currStr.tracksCoordAmpCG,1);%number of subtrajectories
                ncols=size(currStr.tracksCoordAmpCG,2)/8;
                if findex==8
                    endi=find(currStr.seqOfEvents(:,2)==2);
                    starti=find(currStr.seqOfEvents(:,2)==1);
                    fmin=min(currStr.seqOfEvents(:,1));%starting frame
                    for c2=1:nsubtr %subtrajectory identifier within the trajectory
                        currR=currStr.seqOfEvents(starti(c2),3);
                        currEnd=endi(find(currStr.seqOfEvents(endi,3)==currR)); %#ok<FNDSB>
                        frames=(currStr.tracksFeatIndxCG(currR,:)~=0);
                        framesi=find(frames);
                        m=length(framesi);
                        dt=currStr.seqOfEvents(currEnd,1)-currStr.seqOfEvents(starti(c2),1);
                        if (m>nmin && dt>tmin)
                            k=k+1;
                            currStr.fpath=fpath;
                            %       rootn=[rootn extn];
                            rootn=[root 't' int2str(count) 'r' int2str(currR)];
                            if isfinite(currStr.seqOfEvents(starti(c2),4))
                                rootn=[rootn 's' int2str(currStr.seqOfEvents(starti(c2),4))]; %#ok<*AGROW>
                            end
                            if isfinite(currStr.seqOfEvents(currEnd,4))
                                rootn=[rootn 'm' int2str(currStr.seqOfEvents(currEnd,4))];
                            end
                            %                    tracks(k)=currStr;
                            tracks(k).fname = [rootn '.utm'];
                            tracks(k).pos = [currStr.tracksCoordAmpCG(currR,1:8:end)' currStr.tracksCoordAmpCG(currR,2:8:end)'];
                            tracks(k).pos = tracks(k).pos(frames,:);
                            tracks(k).dpos = [currStr.tracksCoordAmpCG(currR,5:8:end)' currStr.tracksCoordAmpCG(currR,6:8:end)'];
                            tracks(k).dpos = tracks(k).dpos(frames,:);
                            tracks(k).frames = fmin-1+framesi';
                            tracks(k).fpath = currStr.fpath;
                            tracks(k).amp = currStr.tracksCoordAmpCG(currR,4:8:end)';
                            tracks(k).amp = tracks(k).amp(frames);
                            tracks(k).damp = currStr.tracksCoordAmpCG(currR,8:8:end)';
                            tracks(k).damp = tracks(k).damp(frames);
                            tracks(k).utrackn=count;
                            
                            tracks(k).ParticleIndex = currStr.tracksFeatIndxCG(currR,:)'; %is probably a copy of 'tracksFeatIndxCG';
                            tracks(k).ParticleIndex = tracks(k).ParticleIndex(frames);    
                            
                            for fnm=fieldnames(currStr)'
                                fm=fnm{1};
                                switch(fm)
                                    case  {'tracksCoordAmpCG','seqOfEvents'}%'tracksFeatIndxCG',
                                    otherwise
                                        sz=size(currStr.(fm),2);
                                        if sz==ncols && ~ischar(currStr.(fm))
                                            tracks(k).(fm) = currStr.(fm)(currR,:)';
                                            tracks(k).(fm) = tracks(k).(fm)(frames);
                                        end
                                end
                            end
                            
                        end
                    end
                else
                    ask_divide;
                    [~,tracksout,currStr]=dividetracks(currStr,root,count,fpath,nsubmin,div,nmin,tmin,k,tsigma,onlyending,ampRatioLimit);
                    if ~isempty(tracksout)
                        [tracksout.utrackn]=deal(count);
                    end
                    tracks=catstruct(tracks,tracksout,2);
                    k=k+length(tracksout);
                    if div, tracksFinal=catstruct(tracksFinal,currStr,1); end
                end
            end
            if ~isempty(tracksFinal) && (~exist('sButtonSave','var') || (exist('sButtonSave','var') && ~strcmpi(sButtonSave,'no')))
                [fp,fn,fx]=fileparts(fname);
                global msd_opts;
                if exist('oldpre','var')
                    oldname=[regexprep(fn,'\d{12}','') oldpre fx];
                else
                    oldname=[fn '_old' strdt fx];
                end
                
                if isfield(msd_opts,'newpre')
                    newname=['Channel_1_tracking_result' msd_opts.newpre datestr(now,'yymmddHHMMSS') '.mat'];
                else
                    newname=[fn '.mat'];
                end
                
                if ~strcmpi(fullfile(fp,oldname),newname) && exist(fullfile(fpath,fname),'file') && (strcmpi(pref('get','LuosTrack','SaveOldDir','ask'),'yes') || ~isfield(msd_opts,'newpre'))
                    copyfile(fullfile(fpath,fname),fullfile(fpath,fp,oldname));
                end
                save(fullfile(fpath,fp,newname),'tracksFinal');
            end
        else
            disp([fname ' is not a valid file']);
        end
        clear fileStr tracksFinal;
    end
    %%
    if findex==7 % File MAT containing the variable tracks
        try
            fileStr=load(fullfile(fpath,fname));
            tracks=catstruct(tracks,fileStr.tracks,2);
        catch MloadtracksfromMAT
            display(['problems in ' fname]);
            disp(getReport(MloadtracksfromMAT,'basic'));
        end
        clear fileStr
    end

    if findex==6 % File MAT containing the results of STALL_Analysis; divide the subtrajectories of TAD.
        error('Luin:ImportTracks','Option not implemented in this simplified version of ImportTracks')
    end
%%
end
return