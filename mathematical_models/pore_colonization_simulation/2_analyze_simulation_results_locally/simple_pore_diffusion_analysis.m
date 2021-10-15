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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD DIFFUSION EQUATION SOLUTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% solved on the cluster

load('data/diffusion_equation_results.mat')


%%%%%%%%%%%%%%%%%%
%% PLOT RESULTS %%
%%%%%%%%%%%%%%%%%%

%% Plotting info

% Minimum probability to show on scales
min_prob_to_show = -5; % log10; 5 is around the number of cells that could fit in the cross-sectional area of a pore opening
min_prob_to_show_extended = min_prob_to_show - 3;

% Color scheme for heatmap
% https://colorbrewer2.org/#type=sequential&scheme=PuBu&n=9
my_colors_heatmap_0 = [ 255,247,251;
    236,231,242;
    208,209,230;
    166,189,219;
    116,169,207;
    54,144,192;
    5,112,176;
    4,90,141;
    2,56,88 ]/255;
my_colors_heatmap = zeros( 2*size(my_colors_heatmap_0,1)-1,3 );
for i=1:size(my_colors_heatmap_0,1)
    my_colors_heatmap( 2*i-1,:) = my_colors_heatmap_0(i,:);
    if i<size(my_colors_heatmap_0,1)
        my_colors_heatmap(2*i,:) = mean( my_colors_heatmap_0(i:i+1,:) );
    end
end

% Color scheme for lines
% https://colorbrewer2.org/#type=diverging&scheme=PRGn&n=9
my_colors_lines = [...
    217,240,211;
    166,219,160;
    90,174,97;
    27,120,55;
    100,100,100;
    118,42,131;
    153,112,171;
    194,165,207;
    231,212,232]/255;


%% Heatmap of prob in pore for base t_gen

% Find time index closest to t_generation
[ ~, tgen_index ] = min( abs(t_vals-t_generation) );
prob_in_pore_at_tgen = zeros( numel(diffusion_coefficients_to_test), numel(drift_speeds_to_test) );
for d=1:numel(diffusion_coefficients_to_test)
    for v=1:numel(drift_speeds_to_test)
        this_prob_in_pore = prob_in_pore_mat{ d,v };
        prob_in_pore_at_tgen(d,v) = this_prob_in_pore( tgen_index );
    end
end

% Make plot
make_heatmap_fig( prob_in_pore_at_tgen, min_prob_to_show, ...
    parameter_test_multiplier_vD, drift_speeds_to_test, diffusion_coefficients_to_test, ...
    my_colors_heatmap, dir_figs )


%% Probability over time plots

% Make probability plots
make_prob_plots( prob_in_pore_mat, ...
    t_vals, t_generation, tgen_index, ...
    drift_speed, drift_speeds_to_test, diffusion_coefficient, ...
    diffusion_coefficients_to_test, parameter_test_multiplier_vD, ...
    min_prob_to_show_extended, ...
    my_colors_lines, dir_figs )

