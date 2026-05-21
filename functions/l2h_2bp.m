function [t_hist_h, x_vec_hist_h, Phi_mtx_hist_h, i_e_h, t_e_h, x_e_vec_h, Phi_mtx_e_h, params_h] ...
    = l2h_2bp(t_hist_l, x_vec_hist_l, Phi_mtx_hist_l, i_e_l, t_e_l, x_e_vec_l, Phi_mtx_e_l, params_l)
%==========================================================================
%
% Converts trajectory data for the Two-Body Problem (2BP) from
% Lagrangian coordinates to Hamiltonian coordinates.
%
% MODEL DESCRIPTION:
% The 2BP may be formulated using either Lagrangian position-velocity
% coordinates or Hamiltonian canonical coordinates. In the inertial frame
% and normalized units, the canonical momenta are equal to the inertial
% velocities.
%
% This function maps trajectory data from the Lagrangian state definition
%
%   x_l_vec = [x; y; z; v_x; v_y; v_z]
%
% to the Hamiltonian state definition
%
%   x_h_vec = [x; y; z; p_x; p_y; p_z]
%
% using the inertial-frame momentum relations
%
%   p_x = v_x
%   p_y = v_y
%   p_z = v_z
%
% Therefore, the state transformation is the identity map:
%
%   x_h_vec = T * x_l_vec
%
% where:
%
%   T = I
%
% The transformation is applied to:
%   - Time histories
%   - State histories
%   - State Transition Matrix (STM) histories
%   - Event data
%   - Parameter structure formulation flag
%
% STM TRANSFORMATION:
% Since the Lagrangian-to-Hamiltonian map is linear and time-independent,
% the State Transition Matrix transforms by similarity:
%
%   Phi_h = T * Phi_l * T^{-1}
%
% For the inertial 2BP, T = I, so:
%
%   Phi_h = Phi_l
%
% NOTES:
% - Time histories and event times are unchanged by this transformation.
% - Event indices are unchanged.
% - State histories and STM histories are unchanged for the inertial 2BP.
% - The only practical change is the interpretation of the last three
%   state components as canonical momenta instead of velocities.
% - Eigenvalues of the STM are invariant under this similarity
%   transformation.
%
% MODEL ID: 2BP_LAG_INERTIAL_ND -> 2BP_HAM_INERTIAL_ND
%
% Author: G. Montseny
% Date: May 14, 2026
%
% INPUT:                 Description                                   Units
%
%  t_hist_l       -  time history in Lagrangian formulation (Nx1)     [-]
%  x_vec_hist_l   -  Lagrangian state history (Nx6)                   [-]
%  Phi_mtx_hist_l -  Lagrangian STM history (Nx6x6)                   [-]
%  i_e_l          -  event indices                                    [-]
%  t_e_l          -  event times                                      [-]
%  x_e_vec_l      -  Lagrangian event states (Ne x 6)                 [-]
%  Phi_mtx_e_l    -  Lagrangian event STMs (Ne x 6 x 6)               [-]
%  params_l       -  parameter struct in Lagrangian formulation       [-]
%
% OUTPUT:                Description                                   Units
%
%  t_hist_h       -  time history in Hamiltonian formulation (Nx1)    [-]
%  x_vec_hist_h   -  Hamiltonian state history (Nx6)                  [-]
%  Phi_mtx_hist_h -  Hamiltonian STM history (Nx6x6)                  [-]
%  i_e_h          -  event indices                                    [-]
%  t_e_h          -  event times                                      [-]
%  x_e_vec_h      -  Hamiltonian event states (Ne x 6)                [-]
%  Phi_mtx_e_h    -  Hamiltonian event STMs (Ne x 6 x 6)              [-]
%  params_h       -  parameter struct with updated formulation flag   [-]
%                    params_h.model.formulation = 'hamiltonian'
%
%==========================================================================

    % Initialize
    t_hist_h = zeros(size(t_hist_l));
    x_vec_hist_h = zeros(size(x_vec_hist_l));
    Phi_mtx_hist_h = zeros(size(Phi_mtx_hist_l));
    i_e_h = zeros(size(i_e_l));
    t_e_h = zeros(size(t_e_l));
    x_e_vec_h = zeros(size(x_e_vec_l));
    Phi_mtx_e_h = zeros(size(Phi_mtx_e_l));
    N = length(t_hist_l);
    Ne = length(t_e_l);

    % Variables that do not change 
    i_e_h = i_e_l;
    t_hist_h = t_hist_l;
    t_e_h = t_e_l;

    % Calculate Transformation matrix T
    T_mtx = eye(6);
    T_inv_mtx = eye(6);


    % Transform coordinates
    for i = 1 : N
        x_vec_h = T_mtx*x_vec_hist_l(i, :)';
        x_vec_hist_h(i, :) = x_vec_h';
    end

    for i = 1 : Ne
        x_e_vec = T_mtx*x_e_vec_l(i, :)';
        x_e_vec_h(i, :) = x_e_vec';
    end

    % Transform the STM
    for i = 1 : N
        Phi_mtx_l = squeeze(Phi_mtx_hist_l(i, :,:));
        Phi_mtx_hist_h(i, :,:) = T_mtx*Phi_mtx_l*T_inv_mtx;
    end

    for i = 1 : Ne
        Phi_mtx_l = squeeze(Phi_mtx_e_l(i, :,:));
        Phi_mtx_e_h(i, :,:) = T_mtx*Phi_mtx_l*T_inv_mtx;
    end

    % Changing the model of the parameters
    params_h = params_l;
    params_h.model.formulation = 'hamiltonian';

end