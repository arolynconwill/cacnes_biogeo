function [ prob_in_pore, diff_eqn_solution, t_vals, x_vals ] = ...
    compute_diffusion_case( ...
    diffusion_coefficient, drift_speed, delta_t, delta_x, num_disks, pop_initial_pdf )


%% Summary

% This function solves the diffusion equation for the provided set of
% parameters.


%% Step size for solving the diffusion quation

% Length and timescales for solving the diffusion equation (can be more
% fine than delta_x and delta_t)
dx_diffeq = .1;%1; % um % spatial step for diffusion equation solution
dt_diffeq = 0.001;%0.001; % days % time step diffusion equation solution

% Calculate how much these differ from delta_x and delta_t
dx_factor = ceil(delta_x/dx_diffeq);
dt_factor = max(ceil(delta_t/dt_diffeq),3); % a minimum of 3 time steps are required for pdepe in compute_diffusion


%% Compute

% Solve diffusion(-drift) equation
[ prob_in_pore, ~, ~, diff_eqn_solution, t_vals, x_vals ] = compute_diffusion( pop_initial_pdf, ...
    diffusion_coefficient, drift_speed, ...
    delta_t, dt_factor, ...
    delta_x, num_disks, dx_factor );

% Note: These originate on absorbing boundary conditions (since no cells
% can end up there, attempt at normalization results in division by zero).

% pop_final_pdf( isnan( pop_final_pdf ) ) = 0;
% if prob_in_pore<0
%     prob_in_pore = 0;
% end


end