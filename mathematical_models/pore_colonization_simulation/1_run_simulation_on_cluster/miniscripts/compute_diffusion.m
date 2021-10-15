function [ prob_in_pore, prob_in_pore_final, pop_final_pdf, diff_eqn_solution, t_vals, x_vals ] = compute_diffusion( pop_initial_pdf, ...
    diffusion_coefficient, drift_speed, ...
    delta_t, dt_factor, ...
    delta_x, num_disks, dx_factor )

%% Summary

% Solves drift-diffusion equation
% "Absorbing" boundary condition at the top of the pore
% "Reflecting" boundary condition at the bottom of the pore

% Returns probability density function and fraction of cells left inside pore


%% Check for no diffusion or drift

if diffusion_coefficient == 0 && drift_speed == 0
    pop_diffused = pop_initial_pdf;
    return
end


%% Simulate diffusion and drift

% Discretize time and space
pop_initial_fine = reshape( ...
    (pop_initial_pdf.*(ones(dx_factor,num_disks)/dx_factor)), ...
    [1 numel(pop_initial_pdf)*dx_factor]); % convert to finer spatial scale
x_vals = linspace( 0, num_disks*delta_x, num_disks*dx_factor );
t_vals = linspace( 0, delta_t, dt_factor );

% Numerically solve drift-diffusion equation
m = 0; % dimension
options = odeset('RelTol',1e-10,'AbsTol',1e-16); % defaults: 1e-3, 1e-6
diff_eqn_solution = pdepe(m,@pdefun,@icfun,@bcfun,x_vals,t_vals,options);

% Compute probability density function over disks 
pop_final_fine = diff_eqn_solution(end,:);
pop_final = sum( reshape(pop_final_fine,[dx_factor,num_disks]) ); % convert back to coarser scale
pop_final_pdf = pop_final/sum(pop_final);
% Compute probability in pore
prob_in_pore = sum( diff_eqn_solution,2 );
prob_in_pore_final = prob_in_pore(end); % since some may have diffused out of pore


% % % % % % % % % % % % % % % % %
% Define differential equation! %
% % % % % % % % % % % % % % % % %

% Define diffusion equation
function [c,f,s] = pdefun(x,t,u,dudx)
    c = 1; % diffusion of cells
    f = diffusion_coefficient*dudx;
    s = drift_speed*dudx; % sebum flow
%    s = 0; % no drift
end

% Define initial conditions
function u0 = icfun(x)
    u0 = pop_initial_fine( x==x_vals ); % use initial condition from function input
end

% Define boundary conditions
function [pL,qL,pR,qR] = bcfun(xL,uL,xR,uR,t)
    % top of pore ("left boundary"): absorbing u(x=top)=0
    pL = uL;
    qL = 0;
    % bottom of pore ("right boundary"): reflecting du/dx(x=bottom)=0
    pR = drift_speed*uR;
    qR = 1;
%     pR = 0;
%     qR = 1;
end

end