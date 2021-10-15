%% OVERVIEW

% This script runs a skin competition model where we compete two strains,
% one of which has a big fitness advantage, under varying assumptions.


%% SETUP

% Add functions to path
path(path,'miniscripts')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPETITION SIMULATIONS: "OUR MODEL" %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Runs one instance of our base model

% Pore initialization: how many cells colonize a pore
num_cells_colonizing_a_pore = 1;
% Migration: whether or not cells can get into filled pores
allow_migration_into_pores = 0;
% Nutrient limitation: whether or not growth in pores is limited by sebum
pores_are_nutrient_limited = 1;
pore_gentime_multiplier = 1; % % 1, 1/10, 1/30 multiplier for base doubling time of 30 days (only used in nutrient limited case)
% Surface growth type
surface_growth_options = { 'none', 'slow', 'same' };
surface_growth_type = 'slow';

[ out_surface_frac_timecourse, out_pore_frac_timecourse, out_pore_fracs_final ] =  ...
    competition_main_fcn( num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );


fig=figure(1);
clf(1)
% Fractions in time
subplot(1,2,1)
hold on
box on
plot( out_pore_frac_timecourse, 'LineWidth', 1 )
plot( out_surface_frac_timecourse, 'LineWidth', 1 )
legend( 'in pores (mean)', 'on surface', 'Location', 'SouthEast' )
ylim([0 1])
xlabel('days')
ylabel('fraction fitter strain')
set(gca,'FontSize',14)
hold off
% Final pore populations
subplot(1,2,2)
box on
histogram( out_pore_fracs_final, 0:0.05:1 )
xlabel('fraction fitter strain')
ylabel('number of pores')
set(gca,'FontSize',14)


fig_string = [ 'fig_porecol-' num2str(num_cells_colonizing_a_pore) ...
    '_poremig-' num2str(allow_migration_into_pores) ...
    '_porelim-' num2str(pores_are_nutrient_limited) ...
    '_poregrow-' num2str(pore_gentime_multiplier) ...
    '_surfgrow-' surface_growth_type '.png' ];
print(fig, fig_string, '-r400','-dpng')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPETITION SIMULATIONS: VARY ASSUMPTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Runs num_sims simulations for each set of model assumptions we're
% testing; then makes some summary plots

% Number of simulations per case
num_sims = 25;


%% BASE MODEL

% Pore initialization: how many cells colonize a pore
num_cells_colonizing_a_pore = 1;
% Migration: whether or not cells can get into filled pores
allow_migration_into_pores = 0;
% Nutrient limitation: whether or not growth in pores is limited by sebum
pores_are_nutrient_limited = 1;
pore_gentime_multiplier = 1; % % 1, 1/10, 1/30 multiplier for base doubling time of 30 days (only used in nutrient limited case)
% Surface growth type
surface_growth_options = { 'none', 'slow', 'same' };
surface_growth_type = 'slow';

[ out_surface_frac_endpoints_base, out_surface_frac_timecourse_base ] = competition_nsim_wrapper( num_sims, ...
    num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );

%% BASE + NUMBER OF CELLS COLONIZING

% Pore initialization: how many cells colonize a pore
num_cells_colonizing_a_pore = 2;
% Migration: whether or not cells can get into filled pores
allow_migration_into_pores = 0;
% Nutrient limitation: whether or not growth in pores is limited by sebum
pores_are_nutrient_limited = 1;
pore_gentime_multiplier = 1; % % 1, 1/10, 1/30 multiplier for base doubling time of 30 days (only used in nutrient limited case)
% Surface growth type
surface_growth_options = { 'none', 'slow', 'same' };
surface_growth_type = 'slow';

[ out_surface_frac_endpoints_numcells2, out_surface_frac_timecourse_numcells2 ] = competition_nsim_wrapper( num_sims, ...
    num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );

% Pore initialization: how many cells colonize a pore
num_cells_colonizing_a_pore = 10;
% Migration: whether or not cells can get into filled pores
allow_migration_into_pores = 0;
% Nutrient limitation: whether or not growth in pores is limited by sebum
pores_are_nutrient_limited = 1;
pore_gentime_multiplier = 1; % % 1, 1/10, 1/30 multiplier for base doubling time of 30 days (only used in nutrient limited case)
% Surface growth type
surface_growth_options = { 'none', 'slow', 'same' };
surface_growth_type = 'slow';

[ out_surface_frac_endpoints_numcells10, out_surface_frac_timecourse_numcells10 ] = competition_nsim_wrapper( num_sims, ...
    num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );

%% BASE + HIGHER GROWTH RATES

% Pore initialization: how many cells colonize a pore
num_cells_colonizing_a_pore = 1;
% Migration: whether or not cells can get into filled pores
allow_migration_into_pores = 0;
% Nutrient limitation: whether or not growth in pores is limited by sebum
pores_are_nutrient_limited = 1;
pore_gentime_multiplier = 1/10; % % 1, 1/10, 1/30 multiplier for base doubling time of 30 days (only used in nutrient limited case)
% Surface growth type
surface_growth_options = { 'none', 'slow', 'same' };
surface_growth_type = 'slow';

[ out_surface_frac_endpoints_fast, out_surface_frac_timecourse_fast ] = competition_nsim_wrapper( num_sims, ...
    num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );

% Pore initialization: how many cells colonize a pore
num_cells_colonizing_a_pore = 2;
% Migration: whether or not cells can get into filled pores
allow_migration_into_pores = 0;
% Nutrient limitation: whether or not growth in pores is limited by sebum
pores_are_nutrient_limited = 1;
pore_gentime_multiplier = 1/10; % % 1, 1/10, 1/30 multiplier for base doubling time of 30 days (only used in nutrient limited case)
% Surface growth type
surface_growth_options = { 'none', 'slow', 'same' };
surface_growth_type = 'slow';

[ out_surface_frac_endpoints_fast_numcells2, out_surface_frac_timecourse_fast_numcells2 ] = competition_nsim_wrapper( num_sims, ...
    num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );

% Pore initialization: how many cells colonize a pore
num_cells_colonizing_a_pore = 10;
% Migration: whether or not cells can get into filled pores
allow_migration_into_pores = 0;
% Nutrient limitation: whether or not growth in pores is limited by sebum
pores_are_nutrient_limited = 1;
pore_gentime_multiplier = 1/10; % % 1, 1/10, 1/30 multiplier for base doubling time of 30 days (only used in nutrient limited case)
% Surface growth type
surface_growth_options = { 'none', 'slow', 'same' };
surface_growth_type = 'slow';

[ out_surface_frac_endpoints_fast_numcells10, out_surface_frac_timecourse_fast_numcells10 ] = competition_nsim_wrapper( num_sims, ...
    num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );

% Pore initialization: how many cells colonize a pore
num_cells_colonizing_a_pore = 1;
% Migration: whether or not cells can get into filled pores
allow_migration_into_pores = 0;
% Nutrient limitation: whether or not growth in pores is limited by sebum
pores_are_nutrient_limited = 1;
pore_gentime_multiplier = 1/30; % % 1, 1/10, 1/30 multiplier for base doubling time of 30 days (only used in nutrient limited case)
% Surface growth type
surface_growth_options = { 'none', 'slow', 'same' };
surface_growth_type = 'slow';

[ out_surface_frac_endpoints_ultra, out_surface_frac_timecourse_ultra ] = competition_nsim_wrapper( num_sims, ...
    num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );

%% BASE + EXTRAS

% Pore initialization: how many cells colonize a pore
num_cells_colonizing_a_pore = 1;
% Migration: whether or not cells can get into filled pores
allow_migration_into_pores = 1;
% Nutrient limitation: whether or not growth in pores is limited by sebum
pores_are_nutrient_limited = 1;
pore_gentime_multiplier = 1; % 1, 1/10, 1/30 multiplier for base doubling time of 30 days (only used in nutrient limited case)
% Surface growth type
surface_growth_options = { 'none', 'slow', 'same' };
surface_growth_type = 'slow';

[ out_surface_frac_endpoints_mig, out_surface_frac_timecourse_mig ] = competition_nsim_wrapper( num_sims, ...
    num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );

% Pore initialization: how many cells colonize a pore
num_cells_colonizing_a_pore = 1;
% Migration: whether or not cells can get into filled pores
allow_migration_into_pores = 1;
% Nutrient limitation: whether or not growth in pores is limited by sebum
pores_are_nutrient_limited = 1;
pore_gentime_multiplier = 1/10; % % 1, 1/10, 1/30 multiplier for base doubling time of 30 days (only used in nutrient limited case)
% Surface growth type
surface_growth_options = { 'none', 'slow', 'same' };
surface_growth_type = 'slow';

[ out_surface_frac_endpoints_fast_mig, out_surface_frac_timecourse_fast_mig ] = competition_nsim_wrapper( num_sims, ...
    num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );

% Pore initialization: how many cells colonize a pore
num_cells_colonizing_a_pore = 1;
% Migration: whether or not cells can get into filled pores
allow_migration_into_pores = 0;
% Nutrient limitation: whether or not growth in pores is limited by sebum
pores_are_nutrient_limited = 1;
pore_gentime_multiplier = 1; % % 1, 1/10, 1/30 multiplier for base doubling time of 30 days (only used in nutrient limited case)
% Surface growth type
surface_growth_options = { 'none', 'slow', 'same' };
surface_growth_type = 'same';

[ out_surface_frac_endpoints_surfsame, out_surface_frac_timecourse_surfsame ] = competition_nsim_wrapper( num_sims, ...
    num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );

% Pore initialization: how many cells colonize a pore
num_cells_colonizing_a_pore = 2;
% Migration: whether or not cells can get into filled pores
allow_migration_into_pores = 0;
% Nutrient limitation: whether or not growth in pores is limited by sebum
pores_are_nutrient_limited = 0;
pore_gentime_multiplier = 1; % % 1, 1/10, 1/30 multiplier for base doubling time of 30 days (only used in nutrient limited case)
% Surface growth type
surface_growth_options = { 'none', 'slow', 'same' };
surface_growth_type = 'slow';

[ out_surface_frac_endpoints_numcells2_porecomp, out_surface_frac_timecourse_numcells2_porecomp ] = competition_nsim_wrapper( num_sims, ...
    num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );


%% CONTROLS

% "Island model"
% Pore initialization: how many cells colonize a pore
num_cells_colonizing_a_pore = 1;
% Migration: whether or not cells can get into filled pores
allow_migration_into_pores = 1;
% Nutrient limitation: whether or not growth in pores is limited by sebum
pores_are_nutrient_limited = 0;
pore_gentime_multiplier = 1; % % 1, 1/10, 1/30 multiplier for base doubling time of 30 days (only used in nutrient limited case)
% Surface growth type
surface_growth_options = { 'none', 'slow', 'same' };
surface_growth_type = 'same';

[ out_surface_frac_endpoints_porecomp_mig_surfsame, out_surface_frac_timecourse_porecomp_mig_surfsame ] = competition_nsim_wrapper( num_sims, ...
    num_cells_colonizing_a_pore, allow_migration_into_pores, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );

%%

% "No pores"
% Nutrient limitation: whether or not growth in pores is limited by sebum
pores_are_nutrient_limited = 1;
pore_gentime_multiplier = 1; % % 1, 1/10, 1/30 multiplier for base doubling time of 30 days (only used in nutrient limited case)
% Surface growth type
surface_growth_options = { 'none', 'slow', 'same' };
surface_growth_type = 'same';

[ out_surface_frac_endpoints_nopores, out_surface_frac_timecourse_nopores ] = competition_nsim_wrapper_nopores( num_sims, ...
    pores_are_nutrient_limited, pore_gentime_multiplier, surface_growth_type );


%%%%%%%%%%%%%%%%%%%%%
%% SUMMARY FIGURES %%
%%%%%%%%%%%%%%%%%%%%%

% Colors
my_colors = lines(8);
col_1 = my_colors(1,:);
col_2 = my_colors(2,:);
col_3 = my_colors(3,:);
col_4 = my_colors(4,:);
col_5 = my_colors(5,:);
col_6 = my_colors(6,:);
col_7 = my_colors(7,:);

fig=figure(2);
clf(2)
% Surface fractions in time
hold on
box on
lw = 2;
plot( out_surface_frac_timecourse_base, '-', 'LineWidth', lw*2, 'Color', 'k' )
plot( out_surface_frac_timecourse_numcells2, '-', 'LineWidth', lw, 'Color', col_2 )
plot( out_surface_frac_timecourse_numcells10, '-', 'LineWidth', lw, 'Color', col_3 )
plot( out_surface_frac_timecourse_mig, '-', 'LineWidth', lw, 'Color', col_4 )
plot( out_surface_frac_timecourse_fast, '--', 'LineWidth', lw, 'Color', col_1 )
plot( out_surface_frac_timecourse_fast_numcells2, '--', 'LineWidth', lw, 'Color', col_2 )
plot( out_surface_frac_timecourse_fast_numcells10, '--', 'LineWidth', lw, 'Color', col_3 )
plot( out_surface_frac_timecourse_fast_mig, '--', 'LineWidth', lw, 'Color', col_4 )
plot( out_surface_frac_timecourse_ultra, ':', 'LineWidth', lw, 'Color', col_1 )
plot( out_surface_frac_timecourse_surfsame, ':', 'LineWidth', lw, 'Color', col_5 )
plot( out_surface_frac_timecourse_numcells2_porecomp, '-', 'LineWidth', lw, 'Color', col_6 )
plot( out_surface_frac_timecourse_porecomp_mig_surfsame, '--', 'LineWidth', lw, 'Color', col_6 )
plot( out_surface_frac_timecourse_nopores, '-', 'LineWidth', lw, 'Color', col_7 )
series_names = { 'base_model', ...
    'num_cells=2', ...
    'num_cells=10', ...
    'mig_into_pores=true', ...
    'growth=fast', ...
    'growth=fast__num_cells=2', ...
    'growth=fast__num_cells=10', ...
    'growth=fast__mig_into_pores=true', ...
    'growth=ultrafast', ...
    'surface-growth=same', ...
    'num_cells=2__pore_comp_true', ...
    'mig_into_pores=true__surface-growth=same__pore_comp=true', ...
    'no_pores' };
legend( series_names, 'Location', 'NorthEastOutside', 'FontSize', 12, 'Interpreter', 'none' )
ylim([0 1])
xlim([0 10*365])
xticks(0:365:365*10); xticklabels(0:1:10);
xlabel('years')
ylabel('fraction fitter strain on surface')
set(gca,'FontSize',16)
hold off

fig_string = [ 'fig_surfcoex-allmodels.png' ];
print(fig, fig_string, '-r400','-dpng')

%%

col_1 = 0.5*[ 1 1 1 ];
col_2 = col_1;
col_3 = col_1;
col_4 = col_1;
col_5 = col_1;
col_6 = col_1;
col_7 = col_1;

fig=figure(3);
clf(3)
% Surface fractions in time
hold on
box on
sz=100;
fa=0.2;
scatter( 1*ones(num_sims,1)+rand(num_sims,1)/2-0.25, out_surface_frac_endpoints_base, sz, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', fa )
scatter( 2*ones(num_sims,1)+rand(num_sims,1)/2-0.25, out_surface_frac_endpoints_numcells2, sz, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_2, 'MarkerFaceAlpha', fa )
scatter( 3*ones(num_sims,1)+rand(num_sims,1)/2-0.25, out_surface_frac_endpoints_numcells10, sz, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_3, 'MarkerFaceAlpha', fa )
scatter( 4*ones(num_sims,1)+rand(num_sims,1)/2-0.25, out_surface_frac_endpoints_mig, sz, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_4, 'MarkerFaceAlpha', fa )
scatter( 5*ones(num_sims,1)+rand(num_sims,1)/2-0.25, out_surface_frac_endpoints_fast, sz, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_1, 'MarkerFaceAlpha', fa )
scatter( 6*ones(num_sims,1)+rand(num_sims,1)/2-0.25, out_surface_frac_endpoints_fast_numcells2, sz, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_2, 'MarkerFaceAlpha', fa )
scatter( 7*ones(num_sims,1)+rand(num_sims,1)/2-0.25, out_surface_frac_endpoints_fast_numcells10, sz, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_3, 'MarkerFaceAlpha', fa )
scatter( 8*ones(num_sims,1)+rand(num_sims,1)/2-0.25, out_surface_frac_endpoints_fast_mig, sz, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_4, 'MarkerFaceAlpha', fa )
scatter( 9*ones(num_sims,1)+rand(num_sims,1)/2-0.25, out_surface_frac_endpoints_ultra, sz, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_1, 'MarkerFaceAlpha', fa )
scatter( 10*ones(num_sims,1)+rand(num_sims,1)/2-0.25, out_surface_frac_endpoints_surfsame, sz, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_5, 'MarkerFaceAlpha', fa )
scatter( 11*ones(num_sims,1)+rand(num_sims,1)/2-0.25, out_surface_frac_endpoints_numcells2_porecomp, sz, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_6, 'MarkerFaceAlpha', fa )
scatter( 12*ones(num_sims,1)+rand(num_sims,1)/2-0.25, out_surface_frac_endpoints_porecomp_mig_surfsame, sz, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_6, 'MarkerFaceAlpha', fa )
scatter( 13*ones(num_sims,1)+rand(num_sims,1)/2-0.25, out_surface_frac_endpoints_nopores, sz, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_7, 'MarkerFaceAlpha', fa )
%legend( series_names, 'Location', 'NorthEastOutside', 'FontSize', 12, 'Interpreter', 'none' )
ylim([0 1])
xlim([0 14])
xlabel('model')
ylabel({'fraction of the more fit strain', 'on the surface after 10 years'})
xticks( 1:1:numel(series_names) )
xticklabels( series_names ); xtickangle( 90 );
set(gca,'TickLabelInterpreter','none')
set(gca,'FontSize',16)
ax = get(gca); 
ax.XAxis.FontSize = 12;
hold off

fig_string = [ 'fig_surfcoex-allmodels-endpoints.png' ];
print(fig, fig_string, '-r400','-dpng')
