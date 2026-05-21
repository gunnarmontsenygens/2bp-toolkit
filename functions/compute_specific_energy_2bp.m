function E_hist = compute_specific_energy_2bp(x_vec_hist, mu)
%==========================================================================
%
% Computes the specific mechanical energy history for a trajectory in the
% dimensionless Two-Body Problem (2BP).
%
% MODEL DESCRIPTION:
% The specific mechanical energy is defined as:
%
%       E = (1/2)v^2 - mu/r
%
% where:
%   - r is the distance from the central body
%   - v is the inertial speed
%   - mu is the gravitational parameter
%
% In the Two-Body Problem, the specific mechanical energy is conserved.
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
%  mu         - gravitational parameter                               -
%
% OUTPUT:              Description                                   Units
%
%  E_hist     - specific mechanical energy history                    -
%
%==========================================================================

    % Preallocation
    N = size(x_vec_hist,1);
    E_hist = zeros(N, 1);

    % Loop
    for i = 1 : N

        % Extract positions and velocities
        r_vec = x_vec_hist(i, 1 : 3); r = norm(r_vec);
        v_vec = x_vec_hist(i, 4 : 6); v = norm(v_vec);
        
        % Compute and store energy 
        E_hist(i, :) = 0.5 * v^2  -  mu / r;
    end


end