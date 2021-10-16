%% FIGURE 5 SCRIPT

% Spatial confinement analysis for facial zones
% uses strips AND extracts


%% ENVIRONMENT

% Add my scripts
dir_my_scripts = [ pwd '/' 'scripts/myscripts_denovomuts'];
path(dir_my_scripts,path);

% Dir to save stuff
dir_save = '6_misc_spatial_confinement';
if ~exist( dir_save, 'dir' )
    mkdir( '6_misc_spatial_confinement' )
end


%% DATA

% Cluster info
load('2_snvs/cluster_names')
clusters_all_slst_super = cellfun(@(x) x(1), clusters_all_slst);
% Sample info
SampleNamesLong_all = load('data/cluster_step_variables','SampleNamesLong_final'); SampleNamesLong_all = SampleNamesLong_all.SampleNamesLong_final;
sample_clusters = cellfun(@(x) str2num(x(end-1:end)), SampleNamesLong_all )';
SampleNames_all = load('data/cluster_step_variables','SampleNames_final'); SampleNames_all = SampleNames_all.SampleNames_final;

% Metadata
types_all = load('data/cluster_step_variables','types_final'); types_all = types_all.types_final;
zones_all = load('data/cluster_step_variables','zones_final'); zones_all = zones_all.zones_final;
specimen_number_all = load('data/cluster_step_variables','specimen_number_final'); specimen_number_all = specimen_number_all.specimen_number_final;

% Load mutliplicity of pore specimens (strips and extracts)
load('data/spec_mult.mat')


%% PROCESS DATA

my_clade_list = [ 1 2 5 8 ]; % 1, 2, 5, 8

for c=1:numel(my_clade_list)
    
    my_clade = my_clade_list(c); 
    this_cluster_name = cluster_names{ my_clade };
    % Note: 
    % % Lineage A-2 basically doesn't have strip samples? (only 2)
    % % Lineage A-3 is on the back so no pore strip samples

    % Which pore strip samples are in this lineage
    this_clade_name = cluster_names{ my_clade }; 
    my_specimens = unique( specimen_number_all( ... % from pore strip and in this lineage
        ( ( types_all == 2 | types_all == 1 ) & ( zones_all <= 5 ) & ( sample_clusters' == my_clade ) ) ... % strips or extracts on FACE
        ) );
    % Which pores are also monophyletic (found manually by looking at trees but should write a function to do this)
    if my_clade == 1 % extracts and strips for subject A
        is_monophyletic = [ 196, 197, 201, 209, 214, 215, 216, 217, 218, 220, 221, 226, 229, 235, 236, 239, 240, 241, 242, 245, 248, 249, 253, 254, 264, ... % strips
            1, 2, 3, 4, 5, 9, 11, 15, 29, 33, 43, 266, 272, 274 ]; % extracts % includes samples with single colony only
        not_monophyletic = [ 199, 210, 211, ... % strips
            7, 19, 30, 34, 44, 268, ... % extracts
            31, 41 ]; % exclude hypermutators
    elseif my_clade == 2 % extracts and strips for subject A
        is_monophyletic = [ 222, 232, ... % strips
            7, 13, 14, 16, 19, 20, 23, 25, 29, 32, 33, 34, 38, 39, 40, 44, 265, 267, 272, 275, 279 ]; % extracts % includes samples with single colony only
        not_monophyletic = [ ... % strips
            24, 27, 37, 276, 280 ]; % extracts
    elseif my_clade == 5 % pore strip only for Subject F
        is_monophyletic = [ 101, 103, 105, 107, 113, 123, 140, 148, 160, 161, 168, 172 ];
        not_monophyletic = [ 166, 170 ];
    elseif my_clade == 8 % pore strip only for subject F
        is_monophyletic = [ 102, 103, 106, 109, 114, 115, 120, 128, 129, 132, 141, 142, 143, 145, 147 ];
        not_monophyletic = [ 127 ];
    else
        fprintf(1,'Warning! Pore strip samples not available for this lineage!')
    end
    my_specimens = intersect( intersect( my_specimens, spec_mult_specnums( spec_mult_porenums == 1 ) ), is_monophyletic ); % only single pore specimens
    numel(my_specimens) 

    %%

    % Genetic distance between pore pairs
    % Load clade data
    load( ['2_snvs/' this_cluster_name '/data_' this_cluster_name '.mat'], 'Calls_for_analysis')
    load( ['2_snvs/' this_cluster_name '/data_' this_cluster_name '.mat'], 'specimen_numbers')
    load( ['2_snvs/' this_cluster_name '/data_' this_cluster_name '.mat'], 'anc_nti_goodpos')
    % Find MRCA for each pore
    pore_mrcas = zeros( numel(my_specimens), numel(anc_nti_goodpos) );
    for p=1:numel(my_specimens)
        spec = my_specimens(p);
        spec_indices = find(specimen_numbers==spec);
        pore_calls = Calls_for_analysis(:,spec_indices);
        pore_calls_mode = mode( pore_calls, 2 );
        pore_calls_nonvariable = ( ...
            sum( pore_calls == 1, 2 ) == numel(spec_indices) | ...
            sum( pore_calls == 2, 2 ) == numel(spec_indices) | ...
            sum( pore_calls == 3, 2 ) == numel(spec_indices) | ...
            sum( pore_calls == 4, 2 ) == numel(spec_indices) ...
            );
        % Infer pore mrca
        pore_mrcas( p, pore_calls_nonvariable ) = pore_calls_mode( pore_calls_nonvariable ); % nonvariable positions from mode
        pore_mrcas( p, ~pore_calls_nonvariable ) = anc_nti_goodpos( ~pore_calls_nonvariable ); % variable positions from ancestor of clade
    end
    % Genetic distance matrix
    dist_genetic = zeros( numel(my_specimens) ); % initialize
    for i=1:numel(my_specimens)
        for j=1:numel(my_specimens)
            pair_notN = ( pore_mrcas(i,:) > 0 ) & ( pore_mrcas(j,:) > 0 );
            pair_dmrca = sum( pore_mrcas(i,pair_notN) ~= pore_mrcas(j,pair_notN) );
            dist_genetic(i,j) = pair_dmrca;
            dist_genetic(j,i) = pair_dmrca;
        end
    end

    % Physical distance between pore pairs
    % Find coordinates for each pore
    my_specimens_zone = zeros( numel(my_specimens),1 ); % initialize
    for p=1:numel(my_specimens)
        spec = my_specimens(p);
        my_specimens_zone(p) = mode( zones_all( specimen_number_all == spec ) );
    end
    % Physical distance matrix
    % binary matrix: 0 = same zone; 1 = diff zones
    dist_physical = zeros( numel(my_specimens) ); % initialize
    for i=1:numel(my_specimens)
        for j=1:numel(my_specimens)
            next_dist = ( my_specimens_zone(i)~=my_specimens_zone(j) );
            dist_physical(i,j) = next_dist;
            dist_physical(j,i) = next_dist;
        end
    end


    %% SPATIAL CONFINEMENT ANALYSIS

    % Define error bar width 
    alpha = 0.05; % for 95% confidence interval

    % Range of genetic distances in data
    list_genetic_distances = 1:1:max(max(dist_genetic));

    % Compute spatial confinement
    list_spatial_confinement_obs = zeros( size(list_genetic_distances) );
    list_spatial_confinement_obs_lower = zeros( size(list_genetic_distances) );
    list_spatial_confinement_obs_upper = zeros( size(list_genetic_distances) );
    for d=1:numel(list_genetic_distances)
        % Pairs to evaluate
        pairs_to_keep = ~eye(numel(my_specimens)) & ... % not i,i
            triu(numel(my_specimens)) & ... % i,j but not j,i
            my_specimens_zone<=5 == my_specimens_zone'<=5 & ... % on face
            dist_genetic <= d;
        % Physical distance between those pairs
        dists_to_keep = dist_physical(pairs_to_keep);
        % Number of pairs within or outside same site radius
        num_same_site = sum( dists_to_keep == 0 );
        num_diff_site = sum( dists_to_keep == 1 );
        % Get probability of related pairs being close together with 95% CI
        x = num_same_site;
        ntrial = num_same_site + num_diff_site;
        [phat,pci] = binofit(x,ntrial,alpha);
        list_spatial_confinement_obs(d) = phat;
        list_spatial_confinement_obs_lower(d) = pci(1);
        list_spatial_confinement_obs_upper(d) = pci(2);
    end

    % Null model for spatial confinement
    num_trials = 1000; % number of reshuffles
    list_spatial_confinement_null_all = zeros( num_trials, numel(list_genetic_distances) );
    for ntrial=1:num_trials
        % Reshuffle specimen numbers
        reshuffle = randperm( numel(my_specimens) );
        dist_genetic_shuffled = dist_genetic( reshuffle, reshuffle );
        for d=1:numel(list_genetic_distances)
            % Pairs to evaluate
            pairs_to_keep = ~eye(numel(my_specimens)) & ... % not i,i
                triu(numel(my_specimens)) & ... % i,j but not j,i
                my_specimens_zone<=5 == my_specimens_zone'<=5 & ... % on face
                dist_genetic_shuffled <= d;
            % Physical distance between those pairs
            dists_to_keep = dist_physical(pairs_to_keep);
            % Number of pairs within or outside same site radius
            num_same_site = sum( dists_to_keep == 0 );
            num_diff_site = sum( dists_to_keep == 1 );
            % Get probability of related pairs being close together with 95% CI
            x = num_same_site;
            n = num_same_site + num_diff_site;
            [phat,~] = binofit(x,n,alpha);
            list_spatial_confinement_null_all(ntrial,d) = phat;
        end
    end

    % Get median vals and 95% CI on null model
    list_spatial_confinement_null = median(list_spatial_confinement_null_all,'omitnan'); 
    list_spatial_confinement_null_mean = mean(list_spatial_confinement_null_all,'omitnan'); % 2020.12.11 
    list_spatial_confinement_null_lower = zeros( size( list_spatial_confinement_null ) ); % initialize
    list_spatial_confinement_null_upper = zeros( size( list_spatial_confinement_null ) ); % initialize
    for d=1:numel(list_genetic_distances)
        next_vals = sort( list_spatial_confinement_null_all(:,d), 'ascend' );
        next_vals = next_vals( ~isnan( next_vals ) ); % remove nans
        if isempty(next_vals)
            fprintf(1,num2str(d))
            list_spatial_confinement_null_lower(d) = 0; % 95% CI
            list_spatial_confinement_null_upper(d) = 1; % 95% CI
        else
            list_spatial_confinement_null_lower(d) = next_vals( round(0.025*numel(next_vals)) ); % 95% CI
            list_spatial_confinement_null_upper(d) = next_vals( round(0.975*numel(next_vals)) ); % 95% CI
        end
    end


    %%

    % Make figure
    fig=figure(2);
    clf(2)
    hold on
    box on
    color_obs = rgb('LightBlue');
    color_null = rgb('LightGray');
    lw = 2;
    alpha_val = 0.2;
    % Null model 
    plot(list_genetic_distances,list_spatial_confinement_null, 'Color', color_null, 'LineWidth', lw)
    % Observation
    plot(list_genetic_distances,list_spatial_confinement_obs, 'Color', color_obs, 'LineWidth', lw )
    % Null model uncertainty
    % plot(list_genetic_distances,list_spatial_confinement_null_lower./list_spatial_confinement_null, 'Color', color_null, 'LineWidth', lw )
    % plot(list_genetic_distances,list_spatial_confinement_null_upper./list_spatial_confinement_null, 'Color', color_null, 'LineWidth', lw)
    patch([list_genetic_distances fliplr(list_genetic_distances)], [list_spatial_confinement_null_lower fliplr(list_spatial_confinement_null_upper)], ...
        color_null, 'EdgeColor', color_null, 'FaceAlpha', alpha_val, 'EdgeAlpha', alpha_val, 'LineWidth', lw );
    % Observed uncertainty
    % plot(list_genetic_distances,list_spatial_confinement_obs_lower./list_spatial_confinement_null, 'Color', color_obs, 'LineWidth', lw )
    % plot(list_genetic_distances,list_spatial_confinement_obs_upper./list_spatial_confinement_null, 'Color', color_obs, 'LineWidth', lw )
    patch([list_genetic_distances fliplr(list_genetic_distances)], [list_spatial_confinement_obs_lower fliplr(list_spatial_confinement_obs_upper)], ...
        color_obs, 'EdgeColor', color_obs, 'FaceAlpha', alpha_val, 'EdgeAlpha', alpha_val, 'LineWidth', lw );
    % X
    xlabel('Genetic distance btwn pores \leq{\it d} SNVs') 
    xlim( [ 0.5 max(max(dist_genetic))+0.5 ] ) 
    % Y
    ylabel(['Probability on same facial zone'])
    ylim([0 1])
    % Title
    title(cluster_names_new(my_clade),'Interpreter','none')
    % Legend
    legend( {'null model','observed'} ) 
    % Appearance
    set(gca,'FontSize',20)
    hold off 

    % Save figure
    print([ dir_save '/' 'PatchPlot_SpatialConfinement_Lineage-' num2str(my_clade) '_zones-face' '_raw.png'],'-dpng')

end