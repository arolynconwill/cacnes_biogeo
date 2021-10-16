%% SUMMARY

% This script computes the observed mutation spectrum.



%% Directory setup

% Directory for Lieberman Lab scripts:
dir_scripts_liebermanlab = '../lab_scripts';
path(dir_scripts_liebermanlab,path);

% Where to find cluster data
dir_lineages = '2_snvs';


%% Load data

% Lineage info
load( [ dir_lineages '/' 'cluster_names' ] )
clusters_all_slst_super = cellfun(@(x) x(1), clusters_all_slst);
num_clusters = numel(cluster_names);

% Lineage and sample info
clusters_all_subjects = load('data/cluster_step_variables','clusters_final_subjects'); clusters_all_subjects = clusters_all_subjects.clusters_final_subjects;
SampleNames_all = load('data/cluster_step_variables','SampleNames_final'); SampleNames_all = SampleNames_all.SampleNames_final;
SampleNamesLong_all = load('data/cluster_step_variables','SampleNamesLong_final'); SampleNamesLong_all = SampleNamesLong_all.SampleNamesLong_final;
SampleNamesSimple_all = load('data/cluster_step_variables','SampleNamesSimple_final'); SampleNamesSimple_all = SampleNamesSimple_all.SampleNamesSimple_final;
clusters_all = load('data/cluster_step_variables','clusters_final'); clusters_all = clusters_all.clusters_final;
unclustered_all = load('data/cluster_step_variables','unclustered_final'); unclustered_all = unclustered_all.unclustered_final;
subjects_all = load('data/cluster_step_variables','subjects_final'); subjects_all = subjects_all.subjects_final;
types_all = load('data/cluster_step_variables','types_final'); types_all = types_all.types_final;
specimen_number_all = load('data/cluster_step_variables','specimen_number_final'); specimen_number_all = specimen_number_all.specimen_number_final;

% Pore strip coordinates
load('data/spec_coor') % 'spec_coor_specnums','spec_coor_stripnum','spec_coor_x','spec_coor_y'

% Load specimen multiplicity
load( 'data/spec_mult' ) % spec_mult_specnums, spec_mult_porenums

% Load SLST colors
load( 'data/fig_format_slst_colors.mat' )


%% Analyze by pore strip

% Pore strips to analyze
list_subjects = { 'A', 'A', 'F', 'F' };
list_subjects_nums = [ 1 1 6 6 ];
list_strip_nums = [ 1 2 1 2 ];

for my_strip = 1:numel(list_subjects)

    % Subject and clusters
    this_subject = list_subjects{my_strip}; 
    this_subject_num = list_subjects_nums(my_strip);
    this_strip_num = list_strip_nums(my_strip); 
    clusters_list = find( cell2mat(clusters_all_subjects) == this_subject );
    clusters_list_superslst = clusters_all_slst_super( cell2mat(clusters_all_subjects) == this_subject );

    % Get specimens of interest
    % only keep this subject's pore strips
    specimens_list = unique( specimen_number_all( subjects_all == this_subject_num & types_all == 2 ) ); 
    % remove samples that correspond to more than one pore
    [ ~, indices_for_mult ] = ismember( specimens_list, spec_mult_specnums );
    specimens_mult = spec_mult_porenums( indices_for_mult );
    specimens_list =  specimens_list( specimens_mult == 1 ); 
    % only keep if on this pore strip
    [ ~, indices_for_coor ] = ismember( specimens_list, spec_coor_specnums );
    specimens_stripnum = spec_coor_stripnum( indices_for_coor );
    specimens_list = specimens_list( specimens_stripnum == this_strip_num );


    % Create a table with the number of colonies from each pore strip specimen
    % in each cluster for a given subject

    specimen_cluster_table = zeros( numel(clusters_list)+1, numel(specimens_list) ); % initialize
    % loop through clusters
    for i=1:numel(clusters_list)
        next_cluster = clusters_all{ clusters_list(i) };
        next_cluster_specimens = specimen_number_all( next_cluster );
        specimen_cluster_table(i,:) = arrayfun(@(x) sum( x==next_cluster_specimens ), specimens_list );
    end
    % add unclustered
    specimen_cluster_table(end,:) = arrayfun(@(x) sum( x==specimen_number_all(unclustered_all) ), specimens_list );

    % Get coordinates for these specimens
    [ ~, indices_for_coor ] = ismember( specimens_list, spec_coor_specnums );
    specimens_x = spec_coor_x( indices_for_coor );
    specimens_y = spec_coor_y( indices_for_coor );


    %% Make a plot

    skip_pores_with_one_colony_only = true;

    if this_subject == 'A'
        jitter_factor = 0.15;
    elseif this_subject == 'F'
        if this_strip_num == 1
            jitter_factor = 0.5;
        else
            jitter_factor = 0.15;
        end
    end

    f=figure(1);
    clf(1)
    hold on
    box on
    % Show location of each pore
    dot_size = 250;
    title( [ 'subject ' this_subject ': pore strip ' num2str(this_strip_num) ] ) 
    if this_subject == 'A'
        xlabel( 'x position (au)' )
        ylabel( 'y position (au)' )
    elseif this_subject == 'F'
        xlabel( 'x position (mm)' )
        ylabel( 'y position (mm)' )
    end
    if skip_pores_with_one_colony_only
        scatter( specimens_x(sum(specimen_cluster_table)>1), -specimens_y(sum(specimen_cluster_table)>1), dot_size, 'k' )
    else
        scatter( specimens_x, -specimens_y, dot_size, 'k' )
    end
    x_min = 5*floor(min( specimens_x )/5 );
    x_max = 5*ceil(max( specimens_x )/5 );
    y_min = -5*floor(min( specimens_y )/5 );
    y_max = -5*ceil(max( specimens_y )/5 );
    xlim( [ x_min x_max ] )
    ylim( [ y_max y_min ] )
    xticks( x_min:5:x_max )
    yticks( y_max:5:y_min )
    axis equal
    grid on
    set(gca,'FontSize',20)
    % Draw colonies from each cluster
    colony_size = 50;
    colony_alpha = 0.75;
    % Cluster colors and symbols
    cluster_colors = cell2mat( arrayfun(@(x) slst_list_colors_short(x==slst_list_short,:)', clusters_list_superslst, 'UniformOutput', false ) )'; % by super SLST
    cluster_colors = [ cluster_colors; 0.5 0.5 0.5 ]; % add gray for unclustered
    symbol_list = 'odsp'; % circle, diamond, star, square
    cluster_symbols = [];
    for k=1:numel(clusters_list_superslst)
        cluster_symbols(end+1) = symbol_list( sum( clusters_list_superslst(k) == clusters_list_superslst(1:k) ));
    end
    cluster_symbols(end+1) = symbol_list(1);
    % Plot by specimen
    for i=1:numel(specimens_list)
        this_pore_x = specimens_x(i);
        this_pore_y = -specimens_y(i);
        this_pore_cols = specimen_cluster_table( :,i );
        this_pore_num_cols = sum( this_pore_cols );
        if this_pore_num_cols == 1 && skip_pores_with_one_colony_only
            continue
        end
        if this_pore_num_cols == 1
            jitter = 0*(rand( this_pore_num_cols, 2 )-0.5);
        else
            jitter = jitter_factor*([ sin( (1:1:this_pore_num_cols)*(2*pi/this_pore_num_cols) ); cos( (1:1:this_pore_num_cols)*(2*pi/this_pore_num_cols) ) ]);
            jitter = jitter';
            %jitter = 1*(rand( this_pore_num_cols, 2 )-0.5)
        end
        j=1; % count pores plotted so far to get jitter right
        for k=1:numel(clusters_list)+1
            num = this_pore_cols(k);
            for n=1:num
                scatter( this_pore_x+jitter(j,1), this_pore_y+jitter(j,2), colony_size, char(cluster_symbols(k)), 'filled', ...
                    'MarkerFaceColor', cluster_colors(k,:), 'MarkerEdgeColor', 'k', ...
                    'MarkerFaceAlpha', colony_alpha, 'MarkerEdgeAlpha', 1 )
                j=j+1;
            end
        end
    end
    hold off

    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [12 8]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 12 8]);
    if skip_pores_with_one_colony_only
        print(['3_straintype_porestripvis' '/' 'Subj-' this_subject '_PoreStrip-' num2str(this_strip_num) '_nosingles.png'],'-r800','-dpng')
    else
        print(['3_straintype_porestripvis' '/' 'Subj-' this_subject '_PoreStrip-' num2str(this_strip_num) '.png'],'-r800','-dpng')
    end

    % Draw legend
    figure(2)
    clf(2)
    hold on
    legend_dot_size = 250;
    legend_x = 1;
    legend_x_text = 1.25;
    legend_y = numel(clusters_list)+2;
    xlim( [ 0 5 ] )
    ylim( [ 0 legend_y+1 ] )
    axis off
    for k=1:numel(clusters_list)+1
        scatter( legend_x, legend_y-k, legend_dot_size, char(cluster_symbols(k)), 'filled', ...
            'MarkerFaceColor', cluster_colors(k,:), 'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', colony_alpha, 'MarkerEdgeAlpha', 1 )
        if k<=numel(clusters_list)
            next_lineage_name = cluster_names_new{clusters_list(k)};
            text( legend_x_text, legend_y-k, [ next_lineage_name(end-1:end) ], 'FontSize', 24, 'Interpreter', 'none' );
    %        text( legend_x_text, legend_y-k, [ 'Lineage ' next_lineage_name(end-2:end) ' colony' ], 'FontSize', 20, 'Interpreter', 'none' );
            %text( legend_x_text, legend_y-k, [ 'cluster ' num2str( clusters_list(k) ) ' colony' ], 'FontSize', 20 );
        else
            text( legend_x_text, legend_y-k, [ 'unclustered' ], 'FontSize', 24 );
    %        text( legend_x_text, legend_y-k, [ 'unclustered colony' ], 'FontSize', 20 );
        end
    end
    hold off

    print(['3_straintype_porestripvis' '/' 'Subj-' this_subject '_PoreStrip-' num2str(this_strip_num) '_legend.png'],'-r400','-dpng')

end


%% Make supplemental table

% Make csv
fid = fopen(['3_straintype_porestripvis' '/' 'supp_table_porestrips.csv'],'w');
fprintf(fid,'Subject, Sample, Strip number, X coordinate, Y coordinate, \n');
for i=1:numel(spec_coor_specnums)
    next_spec = spec_coor_specnums(i);
    next_spec_subjects = subjects_all( next_spec == specimen_number_all );
    if ~isempty( next_spec_subjects) % only record specimens that have colonies that passed filtering
        % Get subject
        next_spec_subject = mode( next_spec_subjects );
        next_spec_subject_char = char( next_spec_subject+64 );
        % Add row to csv
        fprintf(fid, [ ...
            next_spec_subject_char ', ' ...
            num2str(next_spec) ', ' ...
            num2str(spec_coor_stripnum(i)) ', ' ...
            num2str(spec_coor_x(i)) ', ' ...
            num2str(spec_coor_y(i)) ', ' ...        
            '\n' ] );
    end
end
fclose(fid);