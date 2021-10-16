%% FIGURE 3 SCRIPT


%% ENVIRONMENT

dir_scripts_aro = [ 'scripts/myscripts_denovomuts' ];
path(dir_scripts_aro,path);


%% DATA

% Cluster info
load('2_snvs/cluster_names')
clusters_all_slst_super = cellfun(@(x) x(1), clusters_all_slst);
% Sample info
SampleNamesLong_all = load('data/cluster_step_variables','SampleNamesLong_final'); SampleNamesLong_all = SampleNamesLong_all.SampleNamesLong_final;
sample_clusters = cellfun(@(x) str2num(x(end-1:end)), SampleNamesLong_all )';
SampleNames_all = load('data/cluster_step_variables','SampleNames_final'); SampleNames_all = SampleNames_all.SampleNames_final;

% Metadata
subjects_all = load('data/cluster_step_variables','subjects_final'); subjects_all = subjects_all.subjects_final;
specimen_number_all = load('data/cluster_step_variables','specimen_number_final'); specimen_number_all = specimen_number_all.specimen_number_final;
types_all = load('data/cluster_step_variables','types_final'); types_all = types_all.types_final;
zones_all = load('data/cluster_step_variables','zones_final'); zones_all = zones_all.zones_final;

% Load mutliplicity of pore specimens (strips and extracts)
load('data/spec_mult.mat')


%% COLLECTOR'S CURVE FOR PORE VS SURFACE SAMPLES (LITERAL SPECIMENS)
% num of lineages vs num of colonies per specimen

% 3D: Compares diversity of pore vs coarse samples

% Implementation notes
% % Counts unclustered samples as another lineage
% % % but technically counts multiple unclustered as same lineage even though this isn't really true
% % Skips extract or strip specimens originating from >1 pore
% % Bins specimens with >= max_num_col_on_xaxis colonies

legend_bool = false;

max_num_col_on_xaxis_default = 5;

my_subjects_sets = { unique(subjects_all), [1], [2], [6], setdiff( unique(subjects_all), [1]) };
my_subjects_sets_names = { 'All', 'SubjectA', 'SubjectB', 'SubjectF', 'All-A' };

for my_set=1%:numel(my_subjects_sets)

    my_subjects = my_subjects_sets{my_set};
    my_specimens = unique( specimen_number_all( ismember(subjects_all,my_subjects) ) );
    max_num_col_per_spec = max( arrayfun(@(x) sum(specimen_number_all==x), my_specimens ) );
    
    max_num_col_on_xaxis = min( max_num_col_on_xaxis_default, max_num_col_per_spec );
    x_tick_labels = [ 2:1:max_num_col_on_xaxis ];
    x_tick_labels = arrayfun(@(x) {num2str(x)}, x_tick_labels );
    x_tick_labels{end} = [ num2str(max_num_col_on_xaxis) '+' ];

    % initialize
    collec_curve_porespec = zeros( max_num_col_per_spec,2 );
    collec_curve_surfspec = zeros( max_num_col_per_spec,2 );
    % loop through specimens
    for i=1:numel(my_specimens)
        next_spec = my_specimens(i);
        next_type = mode(types_all( specimen_number_all == next_spec ));
        next_subject = mode(subjects_all( specimen_number_all == next_spec ));
        if ismember( next_type, [ 1 2 ] )
            is_pore = 1;
            % check if pore multiplicity > 1
            if spec_mult_porenums( spec_mult_specnums == next_spec ) > 1
                continue % skip if more than one pore
            end
        else
            is_pore = 0;
        end
        % clade assignment for each colony:
        next_clades = sample_clusters( specimen_number_all == next_spec );
        % number of colonies
        next_num_colonies = numel( next_clades );
        % number of clades
        next_num_clades = numel( unique( next_clades ) ); % minor modification needed for unclustered
        is_single_clade = ( next_num_clades == 1 );
        % save data
        if is_pore % pore specimen
            if is_single_clade
                collec_curve_porespec( next_num_colonies,1 ) = collec_curve_porespec( next_num_colonies,1 ) + 1;
            else % multiple clades
                collec_curve_porespec( next_num_colonies,2 ) = collec_curve_porespec( next_num_colonies,2 ) + 1;
            end
        else % surface specimen
            if is_single_clade
                collec_curve_surfspec( next_num_colonies,1 ) = collec_curve_surfspec( next_num_colonies,1 ) + 1;
            else % multiple clades
                collec_curve_surfspec( next_num_colonies,2 ) = collec_curve_surfspec( next_num_colonies,2 ) + 1;
            end
        end
        % print outlier pores
        if my_set==1 && is_pore && ~is_single_clade
            fprintf(1,[ 'Subject ' char(64+next_subject) ', Specimen ' num2str(next_spec) ' has multiple lineages (' num2str(next_num_colonies) ' colonies, ' num2str(next_num_clades) ' lineages).\n' ])
        end
    end
    
    % Bin pore that have a lot of colonies, and skip 1 colony
    collec_curve_porespec_trunc = collec_curve_porespec(2:max_num_col_on_xaxis,:);
    collec_curve_porespec_trunc(max_num_col_on_xaxis-1,:) = sum( collec_curve_porespec(max_num_col_on_xaxis:end,: ));
    collec_curve_surfspec_trunc = collec_curve_surfspec(2:max_num_col_on_xaxis,:);
    collec_curve_surfspec_trunc(max_num_col_on_xaxis-1,:) = sum( collec_curve_surfspec(max_num_col_on_xaxis:end,: ));

    % Figure prep
    fs = 32;
    lw = 1.5;
    
    % Figure: surface ONLY
    fig=figure(6);
    clf(6)
    %
    % surface
    hold on
    box on
    b=bar( collec_curve_surfspec_trunc, 'stacked', ...
        'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', lw );
    b(1).CData = 0.25*[ 1 1 1 ]; % single
    b(2).CData = 0.85*[ 1 1 1 ]; % multiple
    % Axes
    ax = gca;
    ax.XAxis.FontSize = fs-6;
    ax.YAxis.FontSize = fs-6;
    % X
    xlabel('# colonies per sample','FontSize',fs)
    xlim( [ 0+0.25 max_num_col_on_xaxis-0.25 ] )
    xticks( 1:1:max_num_col_on_xaxis-1 )
    xticklabels( x_tick_labels )
    % Y
    h_ylabel = ylabel('#{\bf scrape} samples','FontSize',fs);
    ylim( [ 0 max(sum( collec_curve_surfspec(2:end,:),2 ))+3 ] ) 
    % Legend
    if legend_bool
        legend( { 'single lineage', 'multiple lineages'}, 'Location', 'northeastoutside', 'FontSize', fs )
    end
    hold off

    % Save fig
    if legend_bool
        print(fig,[ '4_pore_biogeo/' 'Bar_CollecCurve_PoreSurface_' my_subjects_sets_names{my_set} '_surfaceonly_wleg.png' ],'-dpng')
    else
        print(fig,[ '4_pore_biogeo/' 'Bar_CollecCurve_PoreSurface_' my_subjects_sets_names{my_set} '_surfaceonly.png' ],'-dpng')
    end
    
    % Figure
    fig=figure(7);
    clf(7)
    %
    % pore
    hold on
    box on
    b=bar( collec_curve_porespec_trunc, 'stacked', ...
        'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', lw );
    b(1).CData = 0.25*[ 1 1 1 ]; % single
    b(2).CData = 0.85*[ 1 1 1 ]; % multiple
    % Axes
    ax = gca;
    ax.XAxis.FontSize = fs-6;
    ax.YAxis.FontSize = fs-6;
    % X
    xlabel('# colonies per sample','FontSize',fs)
    xlim( [ 0+0.25 max_num_col_on_xaxis-0.25 ] )
    xticks( 1:1:max_num_col_on_xaxis-1 )
    xticklabels( x_tick_labels )
    % Y
    h_ylabel = ylabel('#{\bf pore} samples','FontSize',fs);
    ylim( [ 0 max(sum( collec_curve_porespec(2:end,:),2 ))+1 ] ) 
    % Legend
    if legend_bool
        legend( { sprintf(' single\n lineage'), sprintf(' multiple\n lineages')}, 'Location', 'northeastoutside', 'FontSize', fs )
        %legend( { 'single lineage', 'multiple lineages'}, 'Location', 'northeastoutside', 'FontSize', fs )
    end
    hold off

    % Save fig
    if legend_bool
        print(fig,[ '4_pore_biogeo/' 'Bar_CollecCurve_PoreSurface_' my_subjects_sets_names{my_set} '_poresonly_wleg.png' ],'-r400','-dpng')
    else
        print(fig,[ '4_pore_biogeo/' 'Bar_CollecCurve_PoreSurface_' my_subjects_sets_names{my_set} '_poresonly.png' ],'-r400','-dpng')
    end
end

% manually make Fig 6/7 windows a little narrower



%% PORE SAME LINEAGE PROBABILITY
% probability pair of pore colonies from same pore is in same lineage; reshuffle
% probability pair of pore colonies from same zone is in same lineage; reshuffle

% 3B: Compares pore vs zone diversity

% Implementation notes
% % Only colonies originating from pores
% % Face only

% Basics
my_subject = 6; % A or F 
my_samples = ( subjects_all == my_subject ) ...
    & ( zones_all <= 6 ); % facial zones only
my_specimen_numbers = specimen_number_all(my_samples);
my_zones = zones_all( my_samples );
my_sample_clades = sample_clusters( my_samples );
num_clades = numel( unique( sample_clusters( my_samples ) ) );

% For reshuffling
num_trials = 10000;
bw = 0.025; % bind width for histogram

% pores
my_specimen_list = unique( my_specimen_numbers );
% generate list of specimen number and clade number
pore_colony_clades = []; % initialize
c=0; % how many colonies
n_pores = 0;
for s=1:numel(my_specimen_list)
    next_spec = my_specimen_list(s);
    next_type = mode(types_all( specimen_number_all == next_spec ));
    if ~ismember( next_type, [ 1 2 ] )
        continue % skip if not a pore specimen
    end
    if spec_mult_porenums( spec_mult_specnums == next_spec ) > 1
        continue % skip if more than one pore
    end
    n_pores = n_pores+1;
    % clade assignment for each colony:
    next_clades = sample_clusters( specimen_number_all == next_spec );
    %next_clades = next_clades( next_clades ~= 0 ); % ignore unclustered samples
    % number of colonies
    next_num_colonies = numel( next_clades );
    for n=1:next_num_colonies
        c=c+1;
        pore_colony_clades(c,1) = next_spec;
        pore_colony_clades(c,2) = next_clades(n);
    end
end
% compute prob that pair of colonies from same pore is in the same clade
same_spec_matrix = pore_colony_clades(:,1) == pore_colony_clades(:,1)';
same_clade_matrix = pore_colony_clades(:,2) == pore_colony_clades(:,2)';
same_spec_same_clade = same_clade_matrix( same_spec_matrix & ~eye(c) );
prob_pore_pairs_observed = sum(same_spec_same_clade)/numel(same_spec_same_clade)
% compute prob when reshuffling
prob_pore_pairs_reshuffled = zeros( num_trials,1 );
for n=1:num_trials
    pore_colony_clades_reshuffled = [ pore_colony_clades( randperm(c),1 ), pore_colony_clades( :,2 ) ];
    same_spec_matrix = pore_colony_clades_reshuffled(:,1) == pore_colony_clades_reshuffled(:,1)';
    same_clade_matrix = pore_colony_clades_reshuffled(:,2) == pore_colony_clades_reshuffled(:,2)';
    same_spec_same_clade = same_clade_matrix( same_spec_matrix & ~eye(c) );
    prob_pore_pairs = sum(same_spec_same_clade)/numel(same_spec_same_clade);
    prob_pore_pairs_reshuffled(n) = prob_pore_pairs;
end
mean( prob_pore_pairs_reshuffled)
% pval
n_pores
p_shuffle_pores = sum(prob_pore_pairs_observed < prob_pore_pairs_reshuffled)/num_trials

% zones
my_specimen_list = unique( my_specimen_numbers );
% generate list of specimen number and clade number
zone_colony_clades = []; % initialize
c=0; % how many colonies
for s=1:numel(my_specimen_list)
    next_spec = my_specimen_list(s);
    next_zone = mode(zones_all( specimen_number_all == next_spec ));
    next_type = mode(types_all( specimen_number_all == next_spec ));
    if ~ismember( next_type, [ 1 2 ] )
        continue % skip if not a pore specimen
    end
    % clade assignment for each colony:
    next_clades = sample_clusters( specimen_number_all == next_spec );
    %next_clades = next_clades( next_clades ~= 0 ); % ignore unclustered samples
    % number of colonies
    next_num_colonies = numel( next_clades );
    for n=1:next_num_colonies
        c=c+1;
        zone_colony_clades(c,1) = next_zone;
        zone_colony_clades(c,2) = next_clades(n);
    end
end
% compute prob that pair of colonies from same zone is in the same clade
same_spec_matrix = zone_colony_clades(:,1) == zone_colony_clades(:,1)';
same_clade_matrix = zone_colony_clades(:,2) == zone_colony_clades(:,2)';
same_spec_same_clade = same_clade_matrix( same_spec_matrix & ~eye(c) );
prob_zone_pairs_observed = sum(same_spec_same_clade)/numel(same_spec_same_clade)
% compute prob when reshuffling
prob_zone_pairs_reshuffled = zeros( num_trials,1 );
for n=1:num_trials
    zone_colony_clades_reshuffled = [ zone_colony_clades( randperm(c),1 ), zone_colony_clades( :,2 ) ];
    same_spec_matrix = zone_colony_clades_reshuffled(:,1) == zone_colony_clades_reshuffled(:,1)';
    same_clade_matrix = zone_colony_clades_reshuffled(:,2) == zone_colony_clades_reshuffled(:,2)';
    same_spec_same_clade = same_clade_matrix( same_spec_matrix & ~eye(c) );
    prob_pore_pairs = sum(same_spec_same_clade)/numel(same_spec_same_clade);
    prob_zone_pairs_reshuffled(n) = prob_pore_pairs;
end
mean( prob_zone_pairs_reshuffled);
% pval
p_shuffle_zones = sum( prob_zone_pairs_observed < prob_zone_pairs_reshuffled )/num_trials


% prob_pore_pairs_observed =
% 
%     0.9144
% 
% 
% mean( prob_pore_pairs_reshuffled)
% 
% ans =
% 
%     0.4427
% 
% n_pores =
% 
%     75
% 
% 
% p_shuffle_pores =
% 
%      0
% 
% 
% prob_zone_pairs_observed =
% 
%     0.5095
% 
% 
% p_shuffle_zones =
% 
%     0.0317
    
    
%%

% fig=figure(8);
% clf(8)
% fs=14;
% lw = 4;
% subplot(2,1,1)
% % pores
% hold on
% box on
% histogram( prob_pore_pairs_reshuffled, 0:bw:1, 'EdgeColor', 'k', 'FaceColor', 0.9*[1 1 1] )
% xlim([0 1])
% xlabel('probability colony pairs from same pore are in same lineage')
% %yticks([])
% set( gca, 'FontSize', fs )
% line( [prob_pore_pairs_observed,prob_pore_pairs_observed], ylim, 'LineWidth', lw, 'Color', 'k' )
% hold off
% % zones
% subplot(2,1,2)
% hold on
% box on
% histogram( prob_zone_pairs_reshuffled, 0:bw:1, 'EdgeColor', 'k', 'FaceColor', 0.9*[1 1 1] )
% xlim([0 1])
% xlabel('probability colony pairs from same facial region are in same lineage')
% %yticks([])
% set( gca, 'FontSize', fs )
% line( [prob_zone_pairs_observed,prob_zone_pairs_observed], ylim, 'LineWidth', lw, 'Color', 'k' )
% hold off
% 
% % Save fig
% print(fig,[ 'Bar_PorePairsReshuffle_Subj-' char(64+my_subject) '.png' ],'-dpng')

legend_bool = false;
color_hist = 0.25*[1 1 1];

% Figure with just pores
fig=figure(9);
clf(9)
fs=32;
lw = 4;
% pores
hold on
box on
h=histogram( prob_pore_pairs_reshuffled, 0:bw:1, 'EdgeColor', 'k', 'FaceAlpha', 1, 'FaceColor', color_hist);
% Axes
ax = gca;
ax.XAxis.FontSize = fs-6;
ax.YAxis.FontSize = fs-6;
xlim([0 1])
xlabel({'fraction of{\bf intrapore}';'colony pairs in the same lineage'},'FontSize',fs)
%xlabel({'probability that colony pairs from the';'{\bf same pore} are in same lineage'})
ylabel('% of simulations','FontSize',fs)
ylim([0 num_trials/4])
yticks([0 num_trials/4])
yticklabels([0 num_trials/4/100])
line( [prob_pore_pairs_observed,prob_pore_pairs_observed], ylim, 'LineWidth', 4, 'Color','k' )
if legend_bool
    legend({'reshuffled','observed'},'Location','northeastoutside','FontSize',fs)
end
hold off

% Save fig
if legend_bool
    print(fig,[ '4_pore_biogeo/' 'Bar_PorePairsReshuffle_Subj-' char(64+my_subject) '_poreonly_wleg.png' ],'-r400','-dpng')
else
    print(fig,[ '4_pore_biogeo/' 'Bar_PorePairsReshuffle_Subj-' char(64+my_subject) '_poreonly.png' ],'-r400','-dpng')
end

% Figure with just pores
fig=figure(10);
clf(10)
fs=32;
lw = 4;
% pores
hold on
box on
histogram( prob_zone_pairs_reshuffled, 0:bw:1, 'EdgeColor', 'k', 'FaceAlpha', 1, 'FaceColor', color_hist )
% Axes
ax = gca;
ax.XAxis.FontSize = fs-6;
ax.YAxis.FontSize = fs-6;
xlim([0 1])
xlabel({'fraction of{\bf intraregion}';'colony pairs in the same lineage'},'FontSize',fs)
%xlabel({'probability that colony pairs from the';'{\bf same facial region} are in same lineage'})
ylabel('% of simulations','FontSize',fs)
ylim([0 num_trials])
yticks([0 num_trials])
yticklabels([0 num_trials/100])
line( [prob_zone_pairs_observed,prob_zone_pairs_observed], ylim, 'LineWidth', 4, 'Color', 'k' )
if legend_bool
    legend({' reshuffled',' observed'},'Location','northeastoutside','FontSize',fs)
end
hold off

% Save fig
if legend_bool
    print(fig,[ '4_pore_biogeo/' 'Bar_PorePairsReshuffle_Subj-' char(64+my_subject) '_zoneonly_wleg.png' ],'-r400','-dpng')
else
    print(fig,[ '4_pore_biogeo/' 'Bar_PorePairsReshuffle_Subj-' char(64+my_subject) '_zoneonly.png' ],'-r400','-dpng')
end
