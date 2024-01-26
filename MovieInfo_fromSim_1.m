function [movieInfo]= ...
    MovieInfo_fromSim_1(pix_size, A, dA, dr, spotsXframe, dir2save)

%INPUTS
% pix_size: pixel size. Mandatory.
% A: if provided, the amplitude is assigned with a gaussian distribution
%   with mean A(1) and standard deviation A(2). If not provided, amplitude
%   is set to 1.
% dA: if provided, the amplitude uncertainty is assigned with a gaussian distribution
%   with mean dA(1) and standard deviation dA(2). If not provided, amplitude
%   uncertainty is set to 0
% dr: coordinates uncertainty. dx=dy=dr. If not provided, it is set to 0.
% spotsXframe: coordinates of spots in each frame, outuput of simulations;
%   if not provided, ask for opening.
% dir2save: directory for saving

%OUTPUTS
% MovieInfo: variable describing spots coordinates and amplitudes for each
% frame, analogous to the variable produced by utrack at the end of the
% detection process.

%% inputs
if nargin < 2
    A=[]; dA=[]; dr=[];
elseif nargin < 3
    dA=[]; dr=[];
elseif nargin < 4
    dr=[];
end
if nargin<4 || isempty(dr)
    dr=0;
end
dx=dr; dy=dr;

%% spotsXframe
if nargin < 5 || isempty (spotsXframe)
    [file,path] = uigetfile ('spotsXframe.mat', 'open spotsXframe');
    filepath=fullfile(path, file);
    cd(path)
    dir2save=path;
    allChCoord=load(filepath);
else
    allChCoord=struct('channel_1',spotsXframe(1),...
        'channel_2',spotsXframe(2));
end
chNames=fieldnames(allChCoord);

%% construct movieInfo for each channel
%extract number of frames
name=chNames(1);
nfr=length(getfield(allChCoord,name{1}));
%inizialize variables
movieInfo(nfr,1) = struct('xCoord', [], 'yCoord', [], 'amp', []);

name=chNames(1);
oneChCoord=getfield(allChCoord,name{1});

for fr=1:nfr %for each frame

    nspot=length(oneChCoord{fr});
    %coordinates
    movieInfo(fr).xCoord=oneChCoord{fr}(1,:)'/pix_size + dx*randn(nspot,1);
    movieInfo(fr).yCoord=oneChCoord{fr}(2,:)'/pix_size + dy*randn(nspot,1);
    if dr %coordinates' uncertainties
        movieInfo(fr).xCoord(:,2)=dx;
        movieInfo(fr).yCoord(:,2)=dy;
    else
        movieInfo(fr).xCoord(:,2)=zeros(nspot,1);
        movieInfo(fr).yCoord(:,2)=zeros(nspot,1);
    end

    if A %amplitude
        movieInfo(fr).amp=randn(nspot,1)*A(2)+A(1);
    else
        movieInfo(fr).amp=ones(nspot,1);
    end
    if dA %amplitude uncertainty
        movieInfo(fr).amp(:,2)=randn(nspot,1)*dA(2)+dA(1);
    else
        movieInfo(fr).amp(:,2)=zeros(nspot,1);
    end
end

%% save movieInfo variables
name=strcat('movieInfo', '_dr',num2str(dr),'.mat');
folder=fullfile(dir2save, name);
save (folder,'movieInfo')

end