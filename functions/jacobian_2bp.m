function A = jacobian_2bp(t, x_vec, params)
%==========================================================================
%
% Computes the Jacobian matrix of the dimensionless Two-Body Problem (2BP)
% equations of motion in an inertial frame, using the Lagrangian
% (position–velocity) formulation.
%
% MODEL DESCRIPTION:
% The 2BP describes the motion of a particle under the gravitational
% influence of a central body. The motion follows Keplerian dynamics and
% evolves in a fixed orbital plane due to conservation of angular momentum.
%
% NORMALIZATION:
% The system is nondimensionalized such that:
%   - Gravitational parameter: mu = GM = 1
%
% STATE DEFINITION (LAGRANGIAN FORM):
%   x_vec = [x; y; z; v_x; v_y; v_z]
%
% where (x,y,z) is the inertial position vector and
% (v_x,v_y,v_z) is the inertial velocity vector.
%
% EQUATIONS OF MOTION:
%
%   x_dot = v_x
%   y_dot = v_y
%   z_dot = v_z
%
%   v_x_dot = -x / r^3
%   v_y_dot = -y / r^3
%   v_z_dot = -z / r^3
%
% where:
%
%   r = sqrt(x^2 + y^2 + z^2)
%
% JACOBIAN MATRIX:
% The Jacobian matrix is defined as:
%
%   A = df/dx
%
% where f(x) is the 2BP vector field. This matrix governs the variational
% dynamics and is used to propagate the State Transition Matrix (STM):
%
%   Phi_dot = A * Phi
%
% NOTES:
% - Although all Keplerian motion is planar, the equations are written in
%   full 3D Cartesian coordinates for compatibility with the toolkit
%   architecture and STM propagation routines.
% - Planar trajectories are recovered by selecting:
%
%       z_0   = 0
%       v_z_0 = 0
%
% - The lower-left block of A corresponds to the Hessian of the
%   gravitational potential:
%
%       U = -1/r
%
% MODEL ID: 2BP_LAG_INERTIAL_ND
%
% Author: G. Montseny
% Date: May 14, 2026
%
% INPUT:               Description                                   Units
%
%  t         -   time (unused, included for compatibility)            [-]
%  x_vec     -   state vector (6x1)                                  [-]
%  params    -   parameter struct (unused)                            [-]
%
% OUTPUT:              Description                                   Units
%
%  A         -   Jacobian matrix of the vector field (6x6)           [-]
%
%==========================================================================

    % Initialization
    x_vec = x_vec(:);

    % Extract values from x_vec
    x = x_vec(1);
    y = x_vec(2);
    z = x_vec(3);
    
    % Calculate important quantities
    r_vec = [x; y; z];
    r = norm(r_vec);
    

    % Calculate lower left components
    df4dx = -1/r^3 + 3*x^2/r^5;
    df4dy = 3*x*y/r^5;
    df4dz = 3*x*z/r^5;

    df5dx = 3*x*y/r^5;
    df5dy = -1/r^3 + 3*y^2/r^5;
    df5dz = 3*y*z/r^5;

    df6dx = 3*x*z/r^5;
    df6dy = 3*y*z/r^5;
    df6dz = -1/r^3 + 3*z^2/r^5;


    % Dynamics matrix A
    A = [0,    0,    0,    1,  0, 0;
         0,    0,    0,    0,  1, 0;
         0,    0,    0,    0,  0, 1;
         df4dx, df4dy, df4dz, 0,  0, 0;
         df5dx, df5dy, df5dz, 0, 0, 0;
         df6dx, df6dy, df6dz, 0,  0, 0];

end