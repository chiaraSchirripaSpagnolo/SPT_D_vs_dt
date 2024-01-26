function [ msd, sigma, N, C, NC, rho] = MSD_RobustMSD( dt, t, pos, fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [msd, sigma] = MSD_RobustMSD( track, opt)
% 
% Calculates the theoretical MSD for a trace, this version consider the
% possibility of missing points along the trajectory
% 
% INPUT
% 
%  dt : time lag for msd calculation (in frames)
%  t : times over which the trajectory is defined (in frames)
%  pos: positions in the trajectory at times t
%  fig: if greater than zero, it plots the trajectory and the results
% 
% OUTPUT
% 
%  msd: Mean Square Displacement for the specified parameters
%  sigma: standard deviation for each point of the msd.
%  N: number of copules of points used in calculating each point of the msd
%  C: Estimated covariance
%  NC: number of elements used for each covariance point
% 
% ATTENTION: this function consider as unit of time the interframe time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:     Alessandro Duci 
% Correction: Stefano Luin
% Copyright:  (c) 2009 Alessandro Duci; 2024 Stefano Luin s.luin<at>sns.it 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Simplified version, without plots or covariance.

msd     = [];
sigma   = [];
N       = [];
C       = [];
NC      = [];

for i = 1:length(dt) 
    dt_ = dt(i);
    % shifted times
    t_ = t + dt_;
    % checking which t_ are in t
    [tf,loc] = ismember(t_, t);
    idx  = find( tf );
    idx_ = loc( idx );

    if length(idx)>=1

        N(i) = length(idx); %#ok<*AGROW>
        dpos = pos( idx_, : ) - pos(idx, : );
        dpos2 = sum( dpos .^2, 2);
        msd(i) = mean(dpos2) ;
        sigma(i) = sqrt( sum( (dpos2 - msd(i)).^2 ) / ( N(i) - 1 )  );
        
    else
        N(i) = 0;
        msd(i) = NaN;
        sigma(i) = NaN;
    end
    
end

% covariance
if nargout>3
    C = [];
    NC = [];
    rho=NaN;
    warning('Luin:MSD_RobustMSD','Covariance not implemented in this simplified version of MSD_RobustMSD')
end

if (fig>0)
    warning('Luin:MSD_RobustMSD','Plots not implemented in this simplified version of MSD_RobustMSD');
end

return