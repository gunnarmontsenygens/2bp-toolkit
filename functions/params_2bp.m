function params = params_2bp()
%==========================================================================
%
% Builds parameter struct for the dimensionless Two-Body Problem (2BP)
% in an inertial frame, using nondimensional units.
%
% MODEL ID: 2BP_LAG_INERTIAL_ND
%
% Author: G. Montseny
% Date: May 14, 2026
%
% OUTPUT:              Description                                   Units
%
%  params     -   parameter struct                                  [-]
%
%==========================================================================


    %-------------------------------
    % Model definition
    %-------------------------------
    params.model.name = '2BP';
    params.model.formulation = 'lagrangian';
    params.model.frame = 'inertial';
    params.model.units = 'nd';

    %-------------------------------
    % Functions
    %-------------------------------
    params.fun.eom = @eom_2bp;

    %-------------------------------
    % ODE options
    %-------------------------------
    params.ode.options = odeset( ...
        'RelTol', 1e-12, ...
        'AbsTol', 1e-12);

end