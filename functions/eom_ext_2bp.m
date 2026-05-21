function dX_dt_vec = eom_ext_2bp(t, X_vec, params)
%==========================================================================
%
% Computes the dimensionless Two-Body Problem (2BP) equations of motion
% in an inertial frame together with the State Transition Matrix (STM)
% dynamics, using the Lagrangian (position–velocity) formulation.
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
% STATE DEFINITION (AUGMENTED SYSTEM):
% The augmented state combines the physical state and the STM:
%
%   X_vec = [x; y; z; v_x; v_y; v_z; vec(Phi)]
%
% where:
%   - x_vec = [x; y; z; v_x; v_y; v_z] is the 6x1 Lagrangian state
%   - Phi is the 6x6 State Transition Matrix
%   - vec(Phi) is the column-wise vectorization of Phi (36x1)
%
% Total state dimension: 6 + 36 = 42
%
% DYNAMICS:
% The physical state evolves according to:
%
%   x_dot = f(x)
%
% and the STM evolves according to:
%
%   Phi_dot = A(x) * Phi
%
% where:
%   - f(x) is the 2BP vector field
%   - A(x) = df/dx is the Jacobian matrix
%
% NOTES:
% - This is the LAGRANGIAN formulation (not Hamiltonian).
% - The system is autonomous and conservative.
% - Although Keplerian motion is planar, the equations and STM are written
%   in full 3D Cartesian coordinates for compatibility with the toolkit.
% - Planar trajectories are recovered by selecting:
%
%       z_0   = 0
%       v_z_0 = 0
%
% - The STM captures the linearized flow and is used for:
%     • stability analysis
%     • monodromy matrix computation
%     • sensitivity analysis
%
% MODEL ID: 2BP_LAG_INERTIAL_ND
%
% Author: G. Montseny
% Date: May 14, 2026
%
% INPUT:               Description                                   Units
%
%  t         -   time (unused, included for ODE solver compatibility) [-]
%  X_vec     -   augmented state (42x1)                              [-]
%  params    -   parameter struct (unused)                            [-]
%
% OUTPUT:              Description                                   Units
%
%  dX_dt_vec -   time derivative of augmented state (42x1)           [-]
%
%==========================================================================

    % Initialization
    X_vec = X_vec(:);
    
    % Extract state and STM vector
    x_vec = X_vec(1:6);
    Phi_vec = X_vec(7:42);

    % EoM
    dx_dt_vec = eom_2bp(t, x_vec, params);

    % STM
    Phi_mtx = reshape(Phi_vec, 6, 6);
    A_t = jacobian_2bp(t, x_vec, params);
    dPhi_dt_mtx = A_t*Phi_mtx;

    % Put vectors back into Y
    dx_dt_vec = dx_dt_vec(:);
    dPhi_dt_vec = dPhi_dt_mtx(:);
    dX_dt_vec = [dx_dt_vec; dPhi_dt_vec];
end