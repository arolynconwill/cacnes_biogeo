function [ out_surface_frac_timecourse ] = ...
    competition_main_fcn_nopores( ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type )


%% COMPETITION SIMULATIONS

% Outline:
% % Two subpopulations with a fitness difference
% % Inhabit skin surface ONLY (no pores)


%% PARAMETERS

% Numbers
num_pores = 10000; % on face
num_cells_in_pore = 50000; % empirical estimate
ratio_surface_pop_to_pore_pop = 1/10; % ??? 1/10 as many cells on surface as in pores ???
num_cells_on_surface = ratio_surface_pop_to_pore_pop*num_pores*num_cells_in_pore; 

% % Migration rates
% not relevant because no pores in this case

% Time
time_step = 1; % 1 day
time_final = 10*365; % duration of simulation = 2 years
num_time_steps = time_final/time_step;

% % Repopulation time
% not relevant because no pores in this case

% Initial conditions
frac_initial = 0.5; % start populations at 50/50 on surface and inside pores


%% GROWTH RATES

fitness_ratio = 1.5; % fitness differencne between two populations

% Pore growth: whether or not nutrients determine max growth rate inside pores
if pores_are_nutrient_limited
    t_gen_pores = pore_gentime_multiplier*30; % days; time for sebum to "refill" the pore
    t_gen_pores_fit = t_gen_pores; % no relevant fitness difference btwn strains
else
    t_gen_pores = log(2)/4.8; % generation time (days) inside pores (anaerobic; fast)
    t_gen_pores_fit = t_gen_pores/fitness_ratio; % fitness difference btwn strains
end

% Surface growth: none, slower than pores, same as pores
if isequal( surface_growth_type, 'none' ) % A. No surface growth
    t_gen_surface = 100*time_final;
elseif isequal( surface_growth_type, 'slow' ) % % B. Slower surface growth
    t_gen_ratio = 5;
    t_gen_surface = t_gen_pores*t_gen_ratio; % growth rate slowdown on surface (aerobic; slow)
elseif isequal( surface_growth_type, 'same' ) % % C. Same surface growth as in pores
    t_gen_surface = t_gen_pores; 
else
    fprintf( 1, [ 'Unrecognized surface_growth_type: ' surface_growth_type '\n' ] )
end
% Get t_gen for fitter pop
t_gen_surface_fit = t_gen_surface/fitness_ratio;


%% SIMULATION

% Data structures to save
out_surface_frac_timecourse = zeros( num_time_steps, 1 );
% out_pore_frac_timecourse = zeros( num_time_steps, 1 );

% 1. Initialize populations
% 1a. Surface
pop_frac_surface = frac_initial;
% % 1b. Pores
% not relevant because no pores in this case

% % LOOP TIME STEPS % %
% one day time steps
for t=1:num_time_steps
    
    % 2. Growth
    % 2a. Surface
    pop_frac_surface_grown = pop_frac_surface*2^(time_step/t_gen_surface_fit) ...
        / ( pop_frac_surface*2^(time_step/t_gen_surface_fit) + (1-pop_frac_surface)*2^(time_step/t_gen_surface) );
    % 2b. Pores
    % not relevant because no pores in this case
    
    %%
    
    % 3. Migration
    % not relevant because no pores in this case

    %%
    
    % 4. Stochastic pore repopulation
    % not relevant because no pores in this case

    %%
    
    % 5. Update population
    pop_frac_surface = pop_frac_surface_grown;
    % Save
    out_surface_frac_timecourse(t) = pop_frac_surface;
    

end
% % END TIME STEPS % %


end