function x0_vec = generate_IC_2bp(e, w, Omega, i, f0)
%==========================================================================
%
% Generates Cartesian initial conditions for the dimensionless Two-Body
% Problem (2BP) from the classical orbital elements.
%
% MODEL DESCRIPTION:
% The Two-Body Problem describes the motion of a particle under the
% gravitational influence of a single central body. The trajectory is a
% conic section fully determined by the classical orbital elements.
%
% NORMALIZATION:
% The system is nondimensionalized such that:
%   - Semi-major axis: a = 1
%   - Gravitational parameter: mu = 1
%
% Under this normalization:
%   - The orbital period is T = 2*pi
%   - The orbital radius is:
%
%         r = (1 - e^2)/(1 + e*cos(f))
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
%  e       - eccentricity                                             -
%  w       - argument of periapsis                                    deg
%  Omega   - longitude of ascending node                              deg
%  i       - inclination                                              deg
%  f0      - true anomaly at initial time                             deg
%
% OUTPUT:              Description                                   Units
%
%  x0_vec  - Cartesian initial state vector                           -
%
%==========================================================================

    % Compute distance
    r = (1 - e^2) / (1+e*cosd(f0));

    % Compute eccentricity vectors
    e_hat = compute_e_hat(i,w,Omega);
    e_hat_perp = compute_e_hat_perp(i,w,Omega);

    % Compute position vector
    r_vec = r * (cosd(f0) * e_hat + sind(f0) * e_hat_perp);

    % Compute velocity vector
    v_vec = (- sind(f0) * e_hat + (e + cosd(f0)) * e_hat_perp) / sqrt(1 - e^2);


    % Return output
    x0_vec = [r_vec(:); v_vec(:)];


    function e_hat = compute_e_hat(i,w,Omega)
    %==========================================================================
    %
    % Calculates normalized eccentricity vector.
    %
    % Author: G. Montseny
    % Date: October 4, 2025
    %
    %
    % INPUT:               Description                                   Units
    %
    %  i    -           inclination                                        deg
    %  w    -           argument of periapsis                              deg
    %  Omega -          longitude of the ascending node                    deg
    %
    % OUTPUT:       
    %    
    %  e_hat    -       normalized eccentricity vector                     -
    %==========================================================================
    e_hat = [cosd(w)*cosd(Omega)-cosd(i)*sind(w)*sind(Omega),
            cosd(w)*sind(Omega)+cosd(i)*sind(w)*cosd(Omega),
            sind(w)*sind(i)];
    end

    function e_hat_perp = compute_e_hat_perp(i,w,Omega)
    %==========================================================================
    %
    % Calculates normalized vector perpendicular to the eccentricity vector.
    %
    % Author: G. Montseny
    % Date: October 4, 2025
    %
    %
    % INPUT:               Description                                   Units
    %
    %  i    -           inclination                                        deg
    %  w    -           argument of periapsis                              deg
    %  Omega -          longitude of the ascending node                    deg
    %
    % OUTPUT:       
    %    
    %  e_hat_perp  - norm. perpendicular vector to the eccentricity vector   -
    %==========================================================================
    e_hat_perp = [-sind(w)*cosd(Omega)-cosd(i)*cosd(w)*sind(Omega),
            -sind(w)*sind(Omega)+cosd(i)*cosd(w)*cosd(Omega),
            cosd(w)*sind(i)];
    end

end