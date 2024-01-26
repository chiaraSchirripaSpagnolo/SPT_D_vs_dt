%% this file reports the functions used in the different steps of the analysis
% (besides original utrack functions)

%% analysis of exact trajectories
%from simulated tracks, write tracks with the chosen 
% time resolution in utrack format
[tracksFinal] = TracksFinal_fromSim_dt(tracks, dt, nFrames);

%% analysis on exact positions
% from the simulated positions, save a variable in a .mat file describing
% spots coordinates for each frame with utrack formatting.
% note: only dr=0 (or not given) is used in this case, the variable is
% saved in the 'movieInfo_dr0.mat' file
[movieInfo]= ...
    MovieInfo_fromSim_1(pix_size, A, dA, dr, spotsXframe, dir2save);
%from the constructed variable, read from the file 'movieInfo_dr0.mat',
% write an analogous variable extracting the positions with the desiderd temporal sampling
[movieInfo_dt]= movieInfo_sim_dt (dt,Nframe);
% The produced .mat file is analogue to those produced by 
% uTrack at the end of the detection process;
% this file could be used as input for the tracking step.

%% analysis on simulated movies
% creation of a movie from simulated positions with background and noise in
% the images
[ outfile, centers, im_nonoise, img ] = RandomMovie_dt_sumFrame(dt, FrameTime, nf, outfile, varargin );
% on the simulated movies utrack detection and tracking
% processing can be applied

%% calculation of diffusion coefficient distribution (on the tracks resulted for all kinds of analysis above)
diffCoeff_IstoMultiSim (pix_size,dt)
