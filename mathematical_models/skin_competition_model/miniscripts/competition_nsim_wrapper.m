function [ out_surface_frac_endpoints, out_surface_frac_timecourse ] = competition_nsim_wrapper( num_sims, ...
    num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type )

% Initialize

out_surface_frac_endpoints = zeros( num_sims,1 );

% Do num_sims simulations and save the final timepoint

for n=1:num_sims
    
    fprintf(1,[ num2str(n) '\n' ] ) % print index of simulation

    [ out_surface_frac_timecourse, ~, ~ ] =  ...
        competition_main_fcn( num_cells_colonizing_a_pore, allow_migration_into_pores, ...
        pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );
    
    out_surface_frac_endpoints(n) = out_surface_frac_timecourse(end);

end

% Also returns the last simulation

end