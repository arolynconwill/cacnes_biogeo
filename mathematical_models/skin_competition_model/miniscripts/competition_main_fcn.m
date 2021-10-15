function [ out_surface_frac_timecourse, out_pore_frac_timecourse, out_pore_fracs_final ] = ...
    competition_main_fcn( num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type )


%% COMPETITION SIMULATIONS

% Outline:
% % Two subpopulations with a fitness difference
% % Inhabit skin surface and skin pores


%% PARAMETERS

% Numbers
num_pores = 10000; % on face
num_cells_in_pore = 50000; % empirical estimate
ratio_surface_pop_to_pore_pop = 1/10; % ??? 1/10 as many cells on surface as in pores ???
num_cells_on_surface = ratio_surface_pop_to_pore_pop*num_pores*num_cells_in_pore; 

% Migration rates
% % Pores to surface:
migration_pores_to_surface = 1/30; % per day; assume pore turnover is one month
migration_pores_to_surface_adjusted = migration_pores_to_surface*num_pores*num_cells_in_pore/num_cells_on_surface; % 1/3 of cells on surface replaced
% Surface to pores: 
migration_surface_to_pores = 0.01; % per day; assume 1/100 of cells are from surface??

% Time
time_step = 1; % 1 day
time_final = 10*365; % duration of simulation = 2 years
num_time_steps = time_final/time_step;

% Repopulation time
t_pore_repop = 365; % on average a pore is destabilized once per year

% Initial conditions
frac_initial = 0.5; % start populations at 50/50 on surface and inside pores


%% GROWTH RATES

fitness_ratio = 1.5; % fitness differencne between two populations

% Pore growth: whether or not nutrients determine max growth rate inside pores
if pores_are_nutrient_limited
    t_gen_pores = pore_gentime_multiplier*30; % days; time for sebum to "refill" the pore
    t_gen_pores_fit = t_gen_pores; % no relevant fitness difference btwn strains
else
%     t_gen_pores = log(2)/4.8; % generation time (days) inside pores (anaerobic; fast)
    t_gen_pores = pore_gentime_multiplier*30; % days; time for sebum to "refill" the pore
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
out_pore_frac_timecourse = zeros( num_time_steps, 1 );

% 1. Initialize populations
% 1a. Surface
pop_frac_surface = frac_initial;
% 1b. Pores
pop_frac_pores = zeros(num_pores,1); % initialize
for n=1:num_pores
    pop_frac_pores(n) = binornd(num_cells_colonizing_a_pore,frac_initial)/num_cells_colonizing_a_pore;
end

% % LOOP TIME STEPS % %
% one day time steps
for t=1:num_time_steps
    
    % 2. Growth
    % 2a. Surface
    pop_frac_surface_grown = pop_frac_surface*2^(time_step/t_gen_surface_fit) ...
        / ( pop_frac_surface*2^(time_step/t_gen_surface_fit) + (1-pop_frac_surface)*2^(time_step/t_gen_surface) );
    % 2b. Pores
    pop_frac_pores_grown = pop_frac_pores*2^(time_step/t_gen_pores_fit) ...
        ./ ( pop_frac_pores*2^(time_step/t_gen_pores_fit) + (1-pop_frac_pores)*2^(time_step/t_gen_pores) );

    %%
    
    % 3. Migration
    % 3a. Migration from pores to surface
    pop_frac_surface_migrated = (1-migration_pores_to_surface_adjusted)*pop_frac_surface_grown + migration_pores_to_surface_adjusted*mean(pop_frac_pores_grown);
    % 3b. Migration from surface to pores
    if allow_migration_into_pores
        pop_frac_pores_migrated = (1-migration_surface_to_pores)*pop_frac_pores_grown + migration_surface_to_pores*pop_frac_surface_grown;
    else
        pop_frac_pores_migrated = pop_frac_pores_grown;
    end

    %%
    
    % 4. Stochastic pore repopulation
    % 4a. Number of pores that get repopulated at this time step
    % Expected number of pores
    num_pores_for_repop_expected = num_pores*(time_step/t_pore_repop);
    % Draw from poisson
    num_pores_for_repop_rand = poissrnd( num_pores_for_repop_expected );
    % Choose specific pores that get repopulated 
    % NOTE: not considering history of pore
    pores_for_repop = randperm( num_pores, num_pores_for_repop_rand );
    % 4b. Repopulate depending on whether or not we're allowing mixed
    % colonization
    pop_frac_pores_repop = pop_frac_pores_migrated; % initialize from prev step
    for n=1:num_pores_for_repop_rand
        % relative abundances in recolonized pore reflects population on surface
        pop_frac_pores_repop(pores_for_repop(n)) = binornd(num_cells_colonizing_a_pore,pop_frac_surface_migrated)/num_cells_colonizing_a_pore;
    end

    %%
    
    % 5. Update population
    pop_frac_surface = pop_frac_surface_migrated;
    pop_frac_pores = pop_frac_pores_repop;
    % Save
    out_surface_frac_timecourse(t) = pop_frac_surface;
    out_pore_frac_timecourse(t) = mean( pop_frac_pores );
    

end
% % END TIME STEPS % %

% Save
out_pore_fracs_final = pop_frac_pores;



end