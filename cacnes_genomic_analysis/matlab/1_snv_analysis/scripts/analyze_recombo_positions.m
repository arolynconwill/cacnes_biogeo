%% Summary

% This script examines positions masked due to suspected recombination...
% ...and makes a supplemental table.



%% Set up directories

dir_tables = '8_tables';
if ~exist( dir_tables, 'dir' )
    mkdir( dir_tables )
end

% Directory for Lieberman Lab scripts:
dir_scripts_liebermanlab = '../lab_scripts';
path(dir_scripts_liebermanlab,path);

% Where to find the reference genome on the Lieberman Lab Dropbox:
dir_ref_genome = [ pwd '/reference_genomes/Pacnes_C1'];


%% Parameters

recombination_block_size=500;



%% LOAD INFORMATION ON SUSPECTED RECOMBINATION


%% Get reference genome length

[~, GenomeLength, ~, ~] = genomestats(dir_ref_genome);


%% Scan through clusters

% Load cluster names
dir_clusters = '2_snvs/';
load( [ dir_clusters 'cluster_names.mat' ],'cluster_names') % cluster names
load( [ dir_clusters 'cluster_names.mat' ],'cluster_names_new') % cluster names
load( [ dir_clusters 'cluster_names.mat' ],'clusters_all') 
cluster_sizes = cellfun(@(x) numel(x), clusters_all);
num_clusters = numel(cluster_names);

% Initialize
p_recombo_by_lineage = {};
p_recombo_by_lineage_samples = {};
p_recombo_by_lineage_samples_simple = {};
p_recombo_all = [];
num_events_total = 0;

% Scan each cluster
for c=1:num_clusters

    % Load data from this cluster
    load( [ dir_clusters cluster_names{c} '/data_' cluster_names{c} '_recombination.mat' ], 'p_involved_in_non_snp_event' )
    load( [ dir_clusters cluster_names{c} '/data_' cluster_names{c} '_recombination.mat' ], 'p_goodpos_preliminary' )
    load( [ dir_clusters cluster_names{c} '/data_' cluster_names{c} '_recombination.mat' ], 'fixedmutation_preliminary' )
    load( [ dir_clusters cluster_names{c} '/data_' cluster_names{c} '_recombination.mat' ], 'SampleNames' )
    load( [ dir_clusters cluster_names{c} '/data_' cluster_names{c} '_recombination.mat' ], 'SampleNamesSimple' )
    load( [ dir_clusters cluster_names{c} '/data_' cluster_names{c} '_recombination.mat' ], 'p' )
    load( [ dir_clusters cluster_names{c} '/data_' cluster_names{c} '_recombination.mat' ], 'Calls' ) % uncomment if need base calls
    load( [ dir_clusters cluster_names{c} '/data_' cluster_names{c} '_recombination.mat' ], 'anc_nti' ) % uncomment if need base calls
    
    % Update variables
    p_recombo_by_lineage{c} = p_involved_in_non_snp_event;
    p_recombo_all = [ p_recombo_all; p_involved_in_non_snp_event ];
    num_events_total = num_events_total + numel(p_involved_in_non_snp_event);

    % Determine samples with fixed mutations at each position
    involved_samples = {};
    for pos=1:numel(p_involved_in_non_snp_event)
        involved_samples{pos} = SampleNames( fixedmutation_preliminary( p_involved_in_non_snp_event(pos)==p,: ) );
        involved_samples_simple{pos} = SampleNamesSimple( fixedmutation_preliminary( p_involved_in_non_snp_event(pos)==p,: ) );
    end
    p_recombo_by_lineage_samples{end+1} = involved_samples;
    p_recombo_by_lineage_samples_simple{end+1} = involved_samples_simple;

    if c==1
        load( [ dir_clusters cluster_names{c} '/data_' cluster_names{c} '_recombination.mat' ], 'p_on_plasmid' )
    end

end



%% BREAK EVENTS INTO DIFFERENT CATEGORIES

% plasmid-associated regions
% two adjacent "SNPs"
% other


%% Remove plasmid-associated positions

% Initialize
p_recombo_by_lineage_noplasmid = {};
p_recombo_noplasmid = [];
p_recombo_plasmid = [];
num_events_plasmid = 0;

% Scan each cluster
for c=1:num_clusters
    p_next = p_recombo_by_lineage{c};
    p_recombo_by_lineage_noplasmid{c} = setdiff( p_next, p_on_plasmid );
    p_recombo_noplasmid = [ p_recombo_noplasmid; setdiff( p_next, p_on_plasmid ) ];
    p_recombo_plasmid = [ p_recombo_plasmid; p_next( ismember( p_next, p_on_plasmid ) ) ];
    num_events_plasmid = num_events_plasmid + numel( intersect( p_recombo_by_lineage{c}, p_on_plasmid ) );
end


%% Investigate remaining positions

% Note: no non-plasmid positions were repeated
bool_no_repeat_pos = numel(p_recombo_noplasmid) == numel(unique(p_recombo_noplasmid));

% How many adjacent pairs?
adj_to_right = [ abs( p_recombo_noplasmid(1:end-1)-p_recombo_noplasmid(2:end) ) == 1; 0 ];
adj_to_left = [ 0; abs( p_recombo_noplasmid(1:end-1)-p_recombo_noplasmid(2:end) ) == 1 ];
num_events_adjacent = sum( adj_to_right | adj_to_left );
p_recombo_noplasmid_adjacent = p_recombo_noplasmid( adj_to_right | adj_to_left );

% Other types
num_events_other = num_events_total - num_events_plasmid - num_events_adjacent;
p_recombo_noplasmid_other = p_recombo_noplasmid( ~( adj_to_right | adj_to_left ) );


%% Histogram

figure(1)
clf(1)
n_bins = 500;
y_max = 125;
% All
subplot(4,1,1)
histogram(p_recombo_all,1:GenomeLength/(n_bins-1):GenomeLength)
title([ 'all events (n=' num2str(num_events_total) ')'])
ylim([0 y_max])
% Plasmid
subplot(4,1,2)
histogram(p_recombo_plasmid,1:GenomeLength/(n_bins-1):GenomeLength)
ylim([0 y_max])
title([ 'plasmid (n=' num2str(num_events_plasmid) ')'])
% Adjacent
subplot(4,1,3)
histogram(p_recombo_noplasmid_adjacent,1:GenomeLength/(n_bins-1):GenomeLength)
ylim([0 y_max])
title([ 'adjacent (n=' num2str(num_events_adjacent) ')'])
% Other
subplot(4,1,4)
histogram(p_recombo_noplasmid_other,1:GenomeLength/(n_bins-1):GenomeLength)
ylim([0 y_max])
xlabel('position on reference genome')
title([ 'other (n=' num2str(num_events_other) ')'])

% Save image
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 8]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 6 8]);
print([dir_tables '/' 'Hist_recombo-summary.png'],'-dpng')


%% MAKE A TABLE

% Initialize
region_counter = 0;
region_lineage = [];
region_pos = {};
region_pos_first = [];
region_samples_lists = {};
region_samples_lists_simple = {};

for c=1:num_clusters

    % Get recombo positions for this lineage
    p_recombo_this_lineage = p_recombo_by_lineage_noplasmid{c};
    p_recombo_this_lineage_samples = p_recombo_by_lineage_samples{c};
    p_recombo_this_lineage_samples_simple = p_recombo_by_lineage_samples_simple{c};

    % Compute how far apart they are along the reference genome
    p_gaps = p_recombo_this_lineage(2:end) - p_recombo_this_lineage(1:end-1);

    % Determine how may regions there are
    p_regions = zeros( size( p_recombo_this_lineage ) ); % initialize
    for pos=1:numel(p_recombo_this_lineage)
        if pos==1
            region_counter = region_counter + 1;
        elseif p_gaps(pos-1) > recombination_block_size
            region_counter = region_counter + 1;
        end
        p_regions(pos) = region_counter;        
    end

    % Save region info
    for r=min(p_regions):1:max(p_regions)
        region_lineage(end+1) = c;
        p_recombo_this_region = p_recombo_this_lineage( p_regions==r );
        region_pos{end+1} = p_recombo_this_region;
        region_pos_first(end+1) = min( p_recombo_this_region );
        region_samples = {}; 
        region_samples_simple = {};
        region_samples_cells = p_recombo_this_lineage_samples( p_regions==r );
        region_samples_cells_simple = p_recombo_this_lineage_samples_simple( p_regions==r );
        for i=1:numel(region_samples_cells)
            region_samples_cell = region_samples_cells{i};
            for j=1:numel(region_samples_cell)
                region_samples{end+1} = region_samples_cell{j};
            end
        end
        for i=1:numel(region_samples_cells_simple)
            region_samples_cell_simple = region_samples_cells_simple{i};
            for j=1:numel(region_samples_cell_simple)
                region_samples_simple{end+1} = region_samples_cell_simple{j};
            end
        end
        region_samples_lists{end+1} = region_samples;
        region_samples_lists_simple{end+1} = region_samples_simple;
    end
end


%% Write a CSV

% Reorder by first position on genome
[~,reorder]=sort(region_pos_first);

% Make csv
fid = fopen([ dir_tables '/supp_table_recombo.csv' ],'w');
fprintf(fid,'Region index, Position on reference genome, Number of positions on reference genome, Lineage name, Sample names, Sample names (simple), \n');
%fprintf(fid,'region_index, first_pos_on_ref_genome, positions_on_ref_genome, num_pos_on_ref_genome, lineage_name, sample_names, sample_names_simple, \n');
for i=1:region_counter
    r=reorder(i);
    fprintf(fid, [ num2str(i) ',' ]);
    %fprintf(fid, [ num2str(region_pos_first(r)) ',' ]);
    next_pos = region_pos{r};
    next_pos_str = [];
    for i=1:numel(next_pos)
        next_pos_str = [ next_pos_str ' ' num2str(next_pos(i)) ];
    end
    fprintf(fid, [ next_pos_str ',' ]);
    fprintf(fid, [ num2str(numel(next_pos)) ',' ]);
    %fprintf(fid, [num2str(region_lineage(r)) ','] );
    fprintf(fid, [cluster_names_new{region_lineage(r)} ','] );
    % regular sample names
    next_pos_samples = region_samples_lists{r};
    next_pos_samples_str = [];
    for i=1:numel(next_pos_samples)
        next_pos_samples_str = [ next_pos_samples_str ' ' next_pos_samples{i} ];
    end
    fprintf(fid, [ next_pos_samples_str ',' ]);
    % regular sample names
    next_pos_samples = region_samples_lists_simple{r};
    next_pos_samples_str = [];
    for i=1:numel(next_pos_samples)
        next_pos_samples_str = [ next_pos_samples_str ' ' next_pos_samples{i} ];
    end
    fprintf(fid, [ next_pos_samples_str ',' ]);
    fprintf(fid, [ '\n' ]);
end


%% Other notes

% L11 = 1793873 - 1793881