function [ outfile, centers, im_nonoise, img ] = RandomMovie_dt_sumFrame(dt, FrameTime, nf, outfile, varargin )
%[ out, centers, im_nonoise ] = RandomMovie( nf, outfile, szs,n,w,amp,gaussnoise,nwh,nbl,bg,rel_da,rel_dw,poiss,whl,bll,maxn,minn)

%Create movie with time resolution dt and frame time FrameTime

%INPUTS:
% dt: time resolution
% FrameTime : frame time
% nf : number of frames
% outfile : path for saving
% varargin : variables defined in randspots.m function

%OUTPUTS ( TODO defined in randspots.m function ? )
% outfile : 
% centers :
% im_nonoise :
% img :

if nargin<3 || isempty(nf)
    nf=1;
end
if nargin>5 && iscell(varargin{2})
    nf=length(varargin{2});
end

nFr_Dt=floor(nf/(dt/FrameTime));
movie=zeros(varargin{1}, varargin{1},nf);
movieDt=zeros(varargin{1}, varargin{1},nFr_Dt);

if nargin<4 || isempty(outfile)
    
    [outfile,pp]=uiputfile('*.tif','image or movie output file');
    if ~isequal(outfile,0)
        outfile=fullfile(pp,outfile);
    end
end

if nargout>1
    centers=cell(1,nf);
    if nargout>2
        im_nonoise=centers;
        if nargout>3
            img=im_nonoise;
        end
    end
    for ii=1:nf
        if nargin>5 && iscell(varargin{2})
            vararg=[varargin(1),{varargin{2}(ii)},varargin(3:end)];
        else
            vararg=varargin;
        end
        
        if nargout>2
            [frame,centers{ii},im_nonoise{ii}]=randspots(vararg{:});
        else
            [frame,centers{ii}]=randspots(vararg{:});
        end
        frame=uint16(frame);
        if ~isequal(outfile,0)
            imwrite(frame,outfile,'tif','WriteMode','append');
        end
        if nargout>3
            img{ii}=frame;
        end
    end
elseif ~isequal(outfile,0)
    for ii=1:nf
        if nargin>5 && iscell(varargin{2})
            vararg=[varargin(1),{varargin{2}(ii)},varargin(3:end)];
        else
            vararg=varargin;
        end
        frame=uint16(randspots(vararg{:}));
        movie(:,:,ii)=frame;
    end
    
    startFr=1; endFr=dt/FrameTime;
for fr=1:nFr_Dt
    movieDt(:,:,fr)=sum(movie(:,:,startFr:endFr),3);
    startFr=startFr+dt/FrameTime;
    endFr=endFr+dt/FrameTime;
    imwrite(uint16(movieDt(:,:,fr)),...
        outfile,'tif','WriteMode','append');
end
end

end

