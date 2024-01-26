function [movieInfo_dt]= movieInfo_sim_dt (dt,Nframe)

% load movieInfo variable and create an analogous variable
% with time resolution dt and total number of frames Nframe

%INPUT
% dt: time resolution
% Nframe : desired number of frames

%OUTPUT
% movieInfo_dt : variable describing spots coordinates and amplitudes for each
%   frame considering the chosed time resolution
%   (the variable format is analogous to that of the variable produced by utrack
%   at the end of the detection process). 

load('movieInfo_dr0.mat','movieInfo' );

movieInfo_dt=struct([]);
fr_read=1;
for fr_write=1:Nframe
    movieInfo_dt=[movieInfo_dt;movieInfo(fr_read)];
    fr_read=fr_read+dt;
end
name=strcat('movieInfo','_dt',num2str(dt),'.mat');
movieInfo=movieInfo_dt;
save (name,'movieInfo')

end