%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMPLE PORE DIFFUSION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup

% Add miniscripts folder
path('miniscripts',path)

% Add directory for saving figures
dir_figs = 'figs';
if ~exist( dir_figs, 'dir' )
    mkdir( dir_figs )
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER ESTIMATIONS AND INITIAL CONDITIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter estimations

% % diffusion_coefficient: describes diffusion of a cell in sebum (um2/day)
% % drift_speed: describes the speed of sebum moving up the pore (um/day)
% % gamma_anaerobic: growth rate in an anaerobic environment (/day)

% Diffusion parameters
diffusion_coefficient = estimate_diffusion_coefficient_of_cell_in_pore(); % um^2/day
drift_speed = estimate_sebum_flow_speed(); % um/day

% Generation time in aerobic conditions
t_generation = estimate_aerobic_generation_time();

% Pore height
pore_height = 1000;
delta_x = 2;
num_disks = pore_height/delta_x;

% Initial population
pop_initial_depth = 10; % um


%% Define initial location of cells

pop_initial_pdf = zeros( 1,num_disks );
pop_initial_pdf(1:1:floor(pop_initial_depth/delta_x)) = 1/(floor(pop_initial_depth/delta_x)); % all in top disk


%% Run and plot base case 

% % Solve diffusion equation
% [ prob_in_pore, diff_eqn_solution, t_vals, x_vals ] = ...
%     compute_diffusion_case( diffusion_coefficient, drift_speed, t_generation, delta_x, num_disks, pop_initial_pdf );
% % [ prob_in_pore, diff_eqn_solution, t_vals, x_vals ] = ...
% %     compute_diffusion_case( diffusion_coefficient, (2^4)*drift_speed, t_generation, delta_x, num_disks, pop_initial_pdf );
% 
% figure(11); scatter( t_vals, prob_in_pore ); set(gca, 'YScale', 'log');


%% Ranges of generation times, diffusion coefficients, and sebum speeds to test

% Parameters to test (tgen)
parameter_test_multiplier_tgen = 2;
parameter_test_factors_tgen = 2.^(-parameter_test_multiplier_tgen:1:parameter_test_multiplier_tgen);
t_generations_to_test = t_generation*parameter_test_factors_tgen;

% Parameters to test (v and D)
parameter_test_multiplier_vD = 4;
parameter_test_factors_vD = 2.^(-parameter_test_multiplier_vD:1:parameter_test_multiplier_vD);
diffusion_coefficients_to_test = diffusion_coefficient*parameter_test_factors_vD;
drift_speeds_to_test = drift_speed*parameter_test_factors_vD;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLVE DIFFUSION EQUATION FOR CASES TO TEST %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve diffusion equation for all cases

% Save data
prob_in_pore_mat = cell( numel(diffusion_coefficients_to_test), numel(drift_speeds_to_test) ); % initialize

fprintf(1,['Total time: ' num2str(max(t_generations_to_test)) ' days \n'])
for d=1:numel(diffusion_coefficients_to_test)
    fprintf(1,['Diffusion coefficient: ' num2str(diffusion_coefficients_to_test(d)) ' um^2/day \n'])
    for v=1:numel(drift_speeds_to_test)
        fprintf(1,['Sebum speed: ' num2str(drift_speeds_to_test(v)) ' um/day \n'])
        [ prob_in_pore_mat{d,v}, ~, t_vals, ~ ] = ...
            compute_diffusion_case( diffusion_coefficients_to_test(d), drift_speeds_to_test(v), max(t_generations_to_test), delta_x, num_disks, pop_initial_pdf );
    end
end


%% Save results

save('diffusion_equation_results.mat', ...
    't_generations_to_test', 'diffusion_coefficients_to_test', 'drift_speeds_to_test', ...
    'prob_in_pore_mat', 't_vals' ); % doesn't include 'diff_eqn_solutions' because it's massive


