%% FIGURE 3 SCRIPT


%% DATA

% Cluster info
load('2_snvs/cluster_names')
clusters_all_slst_super = cellfun(@(x) x(1), clusters_all_slst);
% Sample info
SampleNamesLong_all = load('data/cluster_step_variables','SampleNamesLong_final'); SampleNamesLong_all = SampleNamesLong_all.SampleNamesLong_final;
sample_clusters = cellfun(@(x) str2num(x(end-1:end)), SampleNamesLong_all );

% Metadata
subjects_all = load('data/cluster_step_variables','subjects_final'); subjects_all = subjects_all.subjects_final;
types_all = load('data/cluster_step_variables','types_final'); types_all = types_all.types_final;
zones_all = load('data/cluster_step_variables','zones_final'); zones_all = zones_all.zones_final;

% Keys
type_names={'Extract', 'Strip', 'Scrape'};
zone_list={{'Forehead','Fo','Tz'},{'Chin','Ch','Cn','Ce'},{'RightCheek','Rc'},{'LeftCheek','Lc'},{'Nose','No','no'},{'NearMouth','Mo','nc'},{'Back','Ba'},{'Neck','Nk'},{'Shoulder','Sh'}};
zone_names_face = arrayfun(@(x) zone_list{x}{1}, 1:1:5, 'UniformOutput', false);
zone_shortnames_face = arrayfun(@(x) zone_list{x}{2}, 1:1:5, 'UniformOutput', false);

% Load SLST colors
load( 'data/fig_format_slst_colors.mat' )


%% STRAIN COMPOSITION AT FACIAL SKIN SITES

% 3A: Compares relative abundance of lineages at different facial skin
% sites on a given subject

% Implementation notes:
% % Only facial sites (excluding mouth because not that many samples)
% % Unclustered colonies are gray
% % For now, treating all colonies as independent (even though later we will learn they are not)

% Loop through by subject
my_subjects_list = unique(subjects_all);
for s=1:numel(my_subjects_list)+2%1:numel(my_subjects_list)+2

    % Get data for this subject
    if s<=numel(my_subjects_list)
        my_subject = my_subjects_list(s);
        my_subject_char = char(my_subject+64);
        my_samples = ( subjects_all == my_subject ) & ... % this subject only
            ( zones_all <=5 ) ; % facial sites only
    elseif s==numel(my_subjects_list)+1
        my_subject_char = 'all';
        my_samples = ( zones_all <=5 ) ; % facial sites only
    elseif s==numel(my_subjects_list)+2
        my_subject_char = 'all-A';
        my_samples = ( zones_all <=5 ) & ... % facial sites only
            ( subjects_all ~= 1 ); % not subject A
    end
    my_samples_lineages = sample_clusters( my_samples );
    my_samples_zones = zones_all( my_samples );
    my_lineages = [setdiff( unique( my_samples_lineages ), [0]) 0]; % always includes unclustered and puts it at the end
    
    % Generate table: number of colonies in each zone from each lineage
    my_table_lineages_zone = zeros( numel(zone_names_face), numel(slst_list_short) );
    for c=1:numel(my_lineages)

        % Get SLST of the lineage
        next_lineage = my_lineages(c);
        if next_lineage == 0 % unclustered
            next_lineage_slst = 'XX';
        else % cluster
            next_lineage_slst = clusters_all_slst{ my_lineages(c) };
        end
        next_lineage_slst_short = next_lineage_slst(1);
        index_next_slst = find( ismember( slst_list_short, next_lineage_slst_short ) );
        
        % Update table for all zones
        for z=1:numel(zone_names_face)
            my_table_lineages_zone(z,index_next_slst) = my_table_lineages_zone(z,index_next_slst) + sum( my_samples_lineages == my_lineages(c) & my_samples_zones == z );
        end
        
    end

    % Bar chart: colony abundance of strains by facial site
    fig=figure(1);
    clf(1)
    fs = 32;
    hold on
    box on
    b = bar(my_table_lineages_zone, 'stacked', 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 2 );
    for k = 1:numel(slst_list_short)
        b(k).CData = slst_list_colors_short(k,:);
    end
    % Axes
    ax = gca;
    ax.XAxis.FontSize = fs-6;
    ax.YAxis.FontSize = fs-6;
    % X
    xticks( 1:1:numel(zone_names_face) )
    xticklabels( zone_shortnames_face )
    xlabel('facial skin region', 'FontSize', fs)
    % Y
    ylabel('number of colonies', 'FontSize', fs)
    hold off

    % Save fig
    print(fig,[ '3_straintype_biogeo' '/BarStrain_FacialZones_Subj-' my_subject_char '_raw.png' ],'-dpng')

    % Bar chart: relative abundance of strains by facial site
    fig=figure(2);
    clf(2)
    fs = 32;
    hold on
    box on
    b = bar(my_table_lineages_zone./sum(my_table_lineages_zone,2), 'stacked', 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 2 );
    for k = 1:numel(slst_list_short)
        b(k).CData = slst_list_colors_short(k,:);
    end
    % Axes
    ax = gca;
    ax.XAxis.FontSize = fs-6;
    ax.YAxis.FontSize = fs-6;
    % X
    xticks( 1:1:numel(zone_names_face) )
    xticklabels( zone_shortnames_face )
    xlabel('facial skin region','FontSize',fs)
    % Y
    ylabel('relative abundance','FontSize',fs)
    hold off

    % Save fig
    print(fig,[ '3_straintype_biogeo' '/BarStrain_FacialZones_Subj-' my_subject_char '_rel.png' ],'-dpng')
    
end


%% STRAIN COMPOSITION ON SURFACE VS IN PORES

% 3B: Compares relative abundance of lineages at different facial skin
% sites on a given subject

% Implementation notes:
% % Only facial sites
% % Unclustered colonies are gray
% % For now, treating all colonies as independent (even though later we will learn they are not)

% Loop through by subject
my_subjects_list = unique(subjects_all);
for s=numel(my_subjects_list)+1:numel(my_subjects_list)+2
    
    % Get data for this subject
    if s<=numel(my_subjects_list)
        my_subject = my_subjects_list(s);
        my_subject_char = char(my_subject+64);
        my_samples = ( subjects_all == my_subject ) & ... % this subject only
            ( zones_all <=5 ) ; % facial sites only
    elseif s==numel(my_subjects_list)+1
        my_subject_char = 'all';
        my_samples = ( zones_all <=5 ) ; % facial sites only
    elseif s==numel(my_subjects_list)+2
        my_subject_char = 'all-A';
        my_samples = ( zones_all <=5 ) & ... % facial sites only
            ( subjects_all ~= 1 ); % not subject A
    end
    my_samples_lineages = sample_clusters( my_samples );
    my_samples_types = types_all( my_samples );
    my_lineages = [setdiff( unique( my_samples_lineages ), [0]) 0]; % always includes unclustered and puts it at the end
    
    % Generate table: number of colonies in each zone from each lineage
    my_table_lineages_poresurf = zeros( 2, numel(slst_list_short) );
    for c=1:numel(my_lineages)
        
        % Get SLST of the lineage
        next_lineage = my_lineages(c);
        if next_lineage == 0 % unclustered
            next_lineage_slst = 'XX';
        else % cluster
            next_lineage_slst = clusters_all_slst{ my_lineages(c) };
        end
        next_lineage_slst_short = next_lineage_slst(1);
        index_next_slst = find( ismember( slst_list_short, next_lineage_slst_short ) );
        
        % Update table
        my_table_lineages_poresurf(1,index_next_slst) = my_table_lineages_poresurf(1,index_next_slst) + sum( my_samples_lineages == my_lineages(c) & ismember(my_samples_types,[ 3 ]) ); % surface
        my_table_lineages_poresurf(2,index_next_slst) = my_table_lineages_poresurf(2,index_next_slst) + sum( my_samples_lineages == my_lineages(c) & ismember(my_samples_types,[ 1 2 ]) ); % pores
    
    end
    
    if min(sum(my_table_lineages_poresurf,2)) < 10 % fewer than 10 colonies for either surface or pore
        continue
    end

    % Bar chart: colony abundance of strains on surface vs in pores
    fig=figure(3);
    clf(3)
    fs = 32;
    hold on
    box on
    b = bar(my_table_lineages_poresurf, 'stacked', 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 2 );
    for k = 1:numel(slst_list_short)
        b(k).CData = slst_list_colors_short(k,:);
    end
    % Axes
    ax = gca;
    ax.XAxis.FontSize = fs-6;
    ax.YAxis.FontSize = fs-6;
    % X
    xlim([0.25 2.75])
    xticks( 1:1:2 )
    xticklabels( {'scrapes', 'pores'} )
    xlabel('sample type','FontSize',fs)
    % Y
    ylabel('number of colonies','FontSize',fs)
    % Formatting
    hold off

    % Save fig
    print(fig,[ '3_straintype_biogeo' '/BarStrain_PoresScrapes_Subj-' my_subject_char '_abs_narrow.png' ],'-dpng')
    
    % Bar chart: relative abundance of strains on surface vs in pores
    fig=figure(4);
    clf(4)
    fs = 32;
    hold on
    box on
    b = bar(my_table_lineages_poresurf./sum(my_table_lineages_poresurf,2), 'stacked', 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 2 );
    for k = 1:numel(slst_list_short)
        b(k).CData = slst_list_colors_short(k,:);
    end
    % Axes
    ax = gca;
    ax.XAxis.FontSize = fs-6;
    ax.YAxis.FontSize = fs-6;
    % X
    xlim([0.25 2.75])
    xticks( 1:1:2 )
    xticklabels( {'scrapes', 'pores'} )
    xlabel('sample type','FontSize',fs)
    % Y
    ylabel('relative abundance','FontSize',fs)
    hold off

    % Save fig
    print(fig,[ '3_straintype_biogeo' '/BarStrain_PoresScrapes_Subj-' my_subject_char '_rel.png' ],'-dpng')
    
end

