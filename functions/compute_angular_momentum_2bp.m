function h_vec_hist = compute_angular_momentum_2bp(x_vec_hist)
%==========================================================================
%
% Computes the specific angular momentum vector history for a trajectory
% in the dimensionless Two-Body Problem (2BP).
%
% MODEL DESCRIPTION:
% The specific angular momentum vector is defined as:
%
%       h_vec = r_vec x v_vec
%
% where:
%   - r_vec is the inertial position vector
%   - v_vec is the inertial velocity vector
%
% In the Two-Body Problem, the angular momentum vector is conserved due
% to the central nature of the gravitational force.
%
% STATE DEFINITION:
% The Cartesian state vector is defined as:
%
%   x_vec = [x; y; z; vx; vy; vz]
%
% where:
%   - [x,y,z]    are the inertial position coordinates
%   - [vx,vy,vz] are the inertial velocity coordinates
%
% Author: G. Montseny
% Date: May 18, 2026
%
% INPUT:               Description                                   Units
%
%  x_vec_hist - trajectory state history                              -
%
% OUTPUT:              Description                                   Units
%
%  h_vec_hist - specific angular momentum vector history              -
%
%==========================================================================

    % Preallocation
    N = size(x_vec_hist,1);
    h_vec_hist = zeros(N, 3);

    % Loop
    for i = 1 : N
        
        % Extract positions and velocities
        r_vec = x_vec_hist(i, 1 : 3);
        v_vec = x_vec_hist(i, 4 : 6); 
        
        % Compute and store angular momentum vector 
        h_vec_hist(i, :) = cross(r_vec, v_vec);
    
    end

end