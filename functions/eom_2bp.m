function dx_dt_vec =  eom_2bp(t, x_vec, params)
%==========================================================================
%
% Computes the dimensionless Two-Body Problem (2BP) equations of motion
% in an inertial frame, using the Lagrangian (position–velocity)
% formulation.
%
% MODEL DESCRIPTION:
% The 2BP describes the motion of a particle under the gravitational
% influence of a central body. The motion follows Keplerian dynamics and
% evolves in a fixed orbital plane due to conservation of angular momentum.
%
% NORMALIZATION:
% The system is nondimensionalized such that:
%   - Gravitational parameter: mu = GM = 1
%   - Semi-major axis: a = 1
%
% STATE DEFINITION (LAGRANGIAN FORM):
%   x_vec = [x; y; z; v_x; v_y; v_z]
%
% where (x,y,z) is the inertial position vector and
% (v_x,v_y,v_z) is the inertial velocity vector.
%
%
%
% EQUATIONS OF MOTION:
%
%   x_ddot = -x / r^3
%   y_ddot = -y / r^3
%
% where:
%
%   r = sqrt(x^2 + y^2)
%
% NOTES:
% - This is the LAGRANGIAN (position–velocity) formulation of the 2BP.
% - The system is autonomous and conservative.
% - The motion is inherently planar.
% - Initial conditions should be selected such that:
%
%       a = 1
%
%   where a is the semi-major axis of the Keplerian orbit.
%
% - The equivalent Hamiltonian formulation is:
%
%       H = 1/2*(p_x^2 + p_y^2) - 1/r
%
% MODEL ID: 2BP_LAG_INERTIAL_ND
%
% Author: G. Montseny
% Date: May 14, 2026
%
% INPUT:               Description                                   Units
%
%  t         -   time (unused, included for ODE solver compatibility) [-]
%  x_vec     -   state vector (6x1)                                   [-]
%  params    -   parameter struct (unused)                            [-]
%
% OUTPUT:              Description                                   Units
%
%  dx_dt_vec -   time derivative of state (6x1)                       [-]
%
%==========================================================================

    % Initialization
    x_vec = x_vec(:);
    dx_dt_vec = zeros(size(x_vec));

    % Extract values from x_vec
    x = x_vec(1);
    y = x_vec(2);
    z = x_vec(3);
    v_x = x_vec(4);
    v_y = x_vec(5);
    v_z = x_vec(6);

    % Calculate important quantities
    r_vec = [x; y; z];
    r = norm(r_vec);

    % Lagrangian EoM
    dx_dt_vec(1) = v_x;
    dx_dt_vec(2) = v_y;
    dx_dt_vec(3) = v_z;
    dx_dt_vec(4) = - x / r^3;
    dx_dt_vec(5) = - y / r^3;
    dx_dt_vec(6) = - z / r^3;

end