function diffCoeff_IstoMultiSim (pix_size,dt)

%the function calculates the MSD and then the diffusion coefficient
% from the first two points.
% It saves an histograhm of diffusion coefficients of all loaded tracks,
% taking into account the uncertainties on D and weighting each trajectory
% with the number of frames when the spots have been detected. 

%INPUTS
% pix_size: pixel size
% dt: time resolution

tracks=ImportTracks();
pn =cd;
fn=['Tracks_dt', num2str(dt), '_5sim'];
%%
All_D_12=zeros(length(tracks),1);
All_sigmaD=zeros(length(tracks),1);
All_N=All_D_12;
All_dt=All_D_12;
for trIdx=1:length(tracks)
    [ msd, sigma, N] = MSD_RobustMSD( 1:2, tracks(trIdx).frames, tracks(trIdx).pos*pix_size, 0);%pos=nframe*2; t= indice del frame;
    msd_delta=sigma./sqrt(N);%err standard
    D_12=(msd(2)-msd(1))/(4.0*(dt/1000));
    All_D_12(trIdx)=D_12;
    sigmaD=sqrt(msd_delta(1)^2+msd_delta(2)^2)/(4.0*(dt/1000));
    All_sigmaD(trIdx)=sigmaD;
    %positions in micron. dt: frame time
    All_N(trIdx)=tracks(trIdx).n;
    All_dt(trIdx)=tracks(trIdx).dt;
end
makehist([All_D_12,All_sigmaD,All_N],1,fullfile(pn,fn),'D12',2);