%% MOLECULAR CLOCK ANALYSIS

%% Summary

% This script attempts to fit a molecular clock for lineages with samples
% at multiple time points (#SNVs to lineage ancestor vs sampling time).


%% Directory setup

% Directory for Lieberman Lab scripts:
dir_scripts_liebermanlab = '../lab_scripts';
path(dir_scripts_liebermanlab,path);

% Add my scripts
dir_my_scripts = [ pwd '/' 'scripts/myscripts_misc'];
path(dir_my_scripts,path);

% Where to find cluster data
dir_lineages = '2_snvs';

% Names of lineages
load( [ dir_lineages '/' 'cluster_names' ] )
num_lineages = numel(cluster_names);


%% Load coverage data

% For filtering individual calls
min_cov_each_strand_for_call = 3; % cov = coverage % 3 in identify_clusters
% For filtering positions
min_median_coverage_position = 12; %15; % across samples per position 

% Coverage data across all colonies
load('data/coverage_matrix.mat','all_coverage_per_bp');
load('data/coverage_matrix.mat','SampleNames');
is_control = cellfun(@(x) contains(x,'Control'), SampleNames );
all_coverage_per_bp = all_coverage_per_bp( ~is_control,: );
% Sample names across all colonies
SampleNames_parsed = load('1_clustering/data_clusterfiltering.mat','SampleNames_unfiltered'); SampleNames_parsed = SampleNames_parsed.SampleNames_unfiltered;

% For filtering plasmid positions
load('data/plasmid_pos_to_mask.mat','plasmid_pos_on_genome'); % positions on reference genome with homology to plasmid gene content



%% MOLECULAR CLOCK ANALYSIS BY LINEAGE

% Directory to store info
dir_lineage_clocks = '6_misc_molec_clock/molec_clock_by_lineage';
if ~exist([dir_lineage_clocks],'dir')
    mkdir([dir_lineage_clocks])
end
dir_lineage_dmrcas = '6_misc_molec_clock/dmrcas_by_lineage';
if ~exist([dir_lineage_dmrcas],'dir')
    mkdir([dir_lineage_dmrcas])
end

%% Load data for next lineage
% uses both polyfit and fit (originally polyfit; then switched to fit but
% didn't delete the old stuff)

% For saving slopes for fits (samples; normalized by coverage; per year)
slopes_list = zeros( numel(cluster_names),1 ); % initialize
slopes_interval_list = zeros( numel(cluster_names),1 ); % initialize

for c=1:numel(cluster_names)
    
    % Lineage info
    next_lineage = cluster_names{c};
    fprintf(1, [ 'Analyzing molecular clock for lineage ' next_lineage '...\n' ] )

    % Load data
    load([ dir_lineages '/' next_lineage '/' 'data_' next_lineage '.mat' ])

    % Remove outgroup colonies and hypermutator colonies
    samples_to_keep = ~outgroup_isolates & ~( specimen_numbers==31 | specimen_numbers==41 );
    Calls_lineage = Calls_for_analysis( :,samples_to_keep ); 
    specimen_numbers_lineage = specimen_numbers( samples_to_keep ); 
    zones_lineage = zones( samples_to_keep );
    times_lineage = times( samples_to_keep );
    times_lineage = times_lineage - min( times_lineage ); % start at zero
    names_lineage = SampleNames( samples_to_keep );
    names_lineage_simple = SampleNamesSimple( samples_to_keep );

    
    %% Compute dMRCA for each colony
    
    if numel(goodpos)>1

        diff_mrca=((Calls_lineage~=repmat(anc_nti_goodpos,1,sum(samples_to_keep))) & Calls_lineage>0);
        dmrca_by_colony = sum( diff_mrca )';

        % Write csv (for putting branch lengths on trees)
        fid = fopen( [ dir_lineage_dmrcas '/' next_lineage '_dmrcas.csv' ], 'w' );
        for i=1:numel(names_lineage)
            fprintf(fid, [ names_lineage{i} ',' names_lineage_simple{i} ',' num2str(dmrca_by_colony(i)) ',' '\n'] );
        end
        fclose(fid);
        
    end
    
    fprintf( 1, [ next_lineage '\n' ] )
    fprintf( 1, ['Max dMRCA = ' num2str(max(dmrca_by_colony)) '\n' ] )
    

    %% Skip molecular clock fitting if insufficient data for the lineage
    
    if numel(unique(times_lineage))<2 || numel(unique(specimen_numbers_lineage))<10
        fprintf(1, [ 'Insufficient data for lineage ' next_lineage '...\n' ] )        
        continue
    end
    
    
    %% Compute SNVs per mb for each colony
    
    % Coverage info
    [ ~, indices_in_all ] = ismember( names_lineage, SampleNames_parsed );
    coverage_lineage = all_coverage_per_bp( indices_in_all,: );
    coverage_lineage_median = median(coverage_lineage);
    
    % Positions to remove
    indices_in_plasmid = zeros( size( coverage_lineage_median ), 'logical' );
    indices_in_plasmid(plasmid_pos_on_genome) = 1;
    indices_low_cov = ( coverage_lineage_median < min_median_coverage_position );
    indices_pos_to_remove = indices_in_plasmid | indices_low_cov;
    
    % Remove those positions from coverage
    coverage_lineage_filtered = coverage_lineage( :,~indices_pos_to_remove );
    
    % Portion of these positions sufficiently covered per sample
    num_pos_covered_by_colony = sum( coverage_lineage_filtered >= 2*min_cov_each_strand_for_call, 2 );
    num_pos_covered_by_colony_mb = num_pos_covered_by_colony/10^6;
    
    dmrca_by_colony_normalized = dmrca_by_colony./num_pos_covered_by_colony_mb; % SNVs per megabases covered
    
    %% Molecular clock: linear fit to dMRCA vs sampling time for all colonies
    % all colonies (though they are far from "independent" if from the same sample)
    
    % Fit molecular clock
    dataset_name = '_all-cols';
    fit_clock( times_lineage', dmrca_by_colony, dir_lineage_clocks, next_lineage, dataset_name, 1 )
    dataset_name = '_all-cols-norm';
    fit_clock( times_lineage', dmrca_by_colony_normalized, dir_lineage_clocks, next_lineage, dataset_name, 2 )
    

    %% Molecular clock: linear fit to dMRCA vs sampling time for all samples
    % average dMRCA over all colonies in sample
    
    % Get average dMRCA for each specimen (since colonies from the same
    % specimen aren't independent)
    specimen_list = unique(specimen_numbers_lineage);
    num_specimens = numel(specimen_list);
    times_lineage_spec = arrayfun(@(x) mode(times_lineage(x==specimen_numbers_lineage)), specimen_list );
    dmrca_by_spec = arrayfun(@(x) mean(dmrca_by_colony(x==specimen_numbers_lineage)), specimen_list );
    dmrca_by_spec_normalized = arrayfun(@(x) mean(dmrca_by_colony_normalized(x==specimen_numbers_lineage)), specimen_list );

    % Fit molecular clock
    dataset_name = 'avg-samp';
    fit_clock( times_lineage_spec', dmrca_by_spec', dir_lineage_clocks, next_lineage, dataset_name, 3 )
    dataset_name = 'avg-samp-norm';
    [ slope, slope_inverval ] = fit_clock( times_lineage_spec', dmrca_by_spec_normalized', dir_lineage_clocks, next_lineage, dataset_name, 4 );
    slopes_list(c) = slope;
    slopes_interval_list(c) = slope_inverval;
    

    %% Molecular clock: face/back for Subj A
    
    if c<=3
        
        %% Molecular clock: linear fit to dMRCA vs sampling time for face samples
        % average dMRCA over all colonies in sample for face only

        % Get average dMRCA for face sample
        specimen_list = unique(specimen_numbers_lineage( zones_lineage <= 6 ));
        num_specimens = numel(specimen_list);
        times_lineage_spec = arrayfun(@(x) mode(times_lineage(x==specimen_numbers_lineage)), specimen_list );
        dmrca_by_spec = arrayfun(@(x) mean(dmrca_by_colony_normalized(x==specimen_numbers_lineage)), specimen_list );

        % Fit molecular clock
        dataset_name = 'faceonly-norm';
        fit_clock( times_lineage_spec', dmrca_by_spec', dir_lineage_clocks, next_lineage, dataset_name, 5 )


        %% Molecular clock: linear fit to dMRCA vs sampling time for back samples
        % average dMRCA over all colonies in sample for back only

        % Get average dMRCA for each back sample
        specimen_list = unique(specimen_numbers_lineage( zones_lineage == 7 ));
        num_specimens = numel(specimen_list);
        times_lineage_spec = arrayfun(@(x) mode(times_lineage(x==specimen_numbers_lineage)), specimen_list );
        dmrca_by_spec = arrayfun(@(x) mean(dmrca_by_colony_normalized(x==specimen_numbers_lineage)), specimen_list );

        % Fit molecular clock
        dataset_name = 'backonly-norm';
        fit_clock( times_lineage_spec', dmrca_by_spec', dir_lineage_clocks, next_lineage, dataset_name, 6 )

    end
    
    
    %% Molecular clock: split lineage 2 subclades
    
    if c==2
        
        % Load info about subclades
        lin2_table = readtable('data/data_extra/lineage2_subcladeinfo.csv');
        lin2_names = lin2_table.sample_name;
        lin2_subclades = lin2_table.sub_clade_num;
        num_subclades = max( lin2_subclades );

        % Initialize
        times_subclades_by_colony = []; % initialize
        dmrca_subclades_by_colony = []; % initialize
        times_subclades_by_sample = []; % initialize
        dmrca_subclades_by_sample = []; % initialize

        for i=1:num_subclades
            % find subclade samples in lineage data
            subclade_samples = ismember( names_lineage, lin2_names( i==lin2_subclades ) );
            num_subclade_colonies = sum(subclade_samples);
            % get mutation info for these samples
            diff_mrca_subclade = diff_mrca( :,subclade_samples );
            subclade_varpos = ( sum( diff_mrca_subclade,2 ) ~= 0 ) ...
                & ( sum( diff_mrca_subclade,2 ) ~= num_subclade_colonies ); % positions that vary within subclade only
            dmrca_this_subclade = sum( diff_mrca_subclade( subclade_varpos,: ) );
            num_pos_covered_by_colony_mb_this_subclade = num_pos_covered_by_colony_mb( subclade_samples );
            dmrca_this_subclade_norm = dmrca_this_subclade./num_pos_covered_by_colony_mb_this_subclade';
            % get other info for for these samples
            times_this_subclade = times_lineage( subclade_samples );
            spec_this_subclade = specimen_numbers_lineage( subclade_samples );
            spec_this_subclade_list = unique( spec_this_subclade );
            % record info by colony
            times_subclades_by_colony = [ times_subclades_by_colony, [ times_this_subclade ] ];
            dmrca_subclades_by_colony = [ dmrca_subclades_by_colony, [ dmrca_this_subclade_norm ] ];
            % record info by specimen
            % (average across colonies from each specimen in this subclade)
            for n=1:numel( spec_this_subclade_list )
                subclade_sample_bool = ( spec_this_subclade_list(n) == spec_this_subclade );
                times_subclades_by_sample( end+1 ) = mode( times_this_subclade( subclade_sample_bool ) );
                dmrca_subclades_by_sample( end+1 ) = mean( dmrca_this_subclade_norm( subclade_sample_bool ) );
            end
        end

        % Fit molecular clock: colonies case
        dataset_name = 'subclades-cols-norm';
        fit_clock( times_subclades_by_colony', dmrca_subclades_by_colony', dir_lineage_clocks, next_lineage, dataset_name, 5 )

        % Fit molecular clock: samples case
        dataset_name = 'subclades-samp-norm';
        fit_clock( times_subclades_by_sample', dmrca_subclades_by_sample', dir_lineage_clocks, next_lineage, dataset_name, 6 )

    end
    
end


%% Make a bar chart with fits for first three lineages

x_axis_label = 'lineage';
x_tick_labels = { 'A-1', 'A-2', 'A-3' };
dir_save = dir_lineage_clocks;
plot_title = '';
plot_file_name = 'MolecClocl_Bar_avg-samp-norm';
plot_bar_with_ebars( slopes_list(1:3), slopes_interval_list(1:3), ...
    plot_title, x_axis_label, x_tick_labels, dir_save, plot_file_name )


