%% SUMMARY

% Compare strain composition on the face vs the back

% Implementation notes
% % Ignoring unclustered colonies
% % Treating colonies as "independent"


%% Parameters

min_num_cols_on_zone = 5;


%% Load data

% Load cluster data (rename when necessary)
clusters_all_subjects = load('data/cluster_step_variables','clusters_final_subjects'); clusters_all_subjects = clusters_all_subjects.clusters_final_subjects;
SampleNames_all = load('data/cluster_step_variables','SampleNames_final'); SampleNames_all = SampleNames_all.SampleNames_final;
SampleNamesSimple_all = load('data/cluster_step_variables','SampleNamesSimple_final'); SampleNamesSimple_all = SampleNamesSimple_all.SampleNamesSimple_final;
SampleNamesLong_all = load('data/cluster_step_variables','SampleNamesLong_final'); SampleNamesLong_all = SampleNamesLong_all.SampleNamesLong_final;
unclustered_all = load('data/cluster_step_variables','unclustered_final'); unclustered_all = unclustered_all.unclustered_final;

% Load metadata
subjects_all = load('data/cluster_step_variables','subjects_final'); subjects_all = subjects_all.subjects_final;
zones_all = load('data/cluster_step_variables','zones_final'); zones_all = zones_all.zones_final;
types_all = load('data/cluster_step_variables','types_final'); types_all = types_all.types_final;
specimen_number_all = load('data/cluster_step_variables','specimen_number_final'); specimen_number_all = specimen_number_all.specimen_number_final;
slst_all = load('data/cluster_step_variables','slst_final'); slst_all = slst_all.slst_final; slst_all_original = slst_all; % save because will be updated with assembley info later

% Load SLST from assemblies
load('2_snvs/cluster_names.mat')
clusters_all_slst_super = cellfun(@(x) x(1), clusters_all_slst );
% Process
sample_clusters = cellfun(@(x) str2num(x(end-1:end)), SampleNamesLong_all );

% Data key
type_names={'Extract', 'Strip', 'Scrape'};
zone_list={{'Forehead','Fo','Tz'},{'Chin','Ch','Cn','Ce'},{'RightCheek','Rc'},{'LeftCheek','Lc'},{'Nose','No','no'},{'NearMouth','Mo','nc'},{'Back','Ba'},{'Neck','Nk'},{'Shoulder','Sh'}};
zone_names_face = arrayfun(@(x) zone_list{x}{1}, 1:1:5, 'UniformOutput', false);
zone_shortnames_face = arrayfun(@(x) zone_list{x}{2}, 1:1:5, 'UniformOutput', false);

% Load SLST colors
load( [ 'data/fig_format_slst_colors.mat' ] )


%% LINEAGE COMPOSITION ON FACE VS BACK

% SUPP: Compares relative abundance of lineages on face vs back by subject

% Implementation notes:
% % Did not collapse pores!!
% % Face = six sites; back = one site; ignoring neck/shoulders
% % Unclustered colonies are gray

% Loop through by subject
my_subjects_list = unique(subjects_all);
for s=1:numel(my_subjects_list) 

    % Get data for this subject
    my_subject = my_subjects_list(s);
    my_samples = ( subjects_all == my_subject );
    my_samples_lineages = sample_clusters( my_samples );
    my_samples_zones = zones_all( my_samples );
    my_samples_types = types_all( my_samples );
    my_lineages = setdiff( unique( my_samples_lineages ), [0]); % do not include unclustered 
    % Generate table: number of colonies from each strain type in each zone
    my_table_strains_faceback = zeros( 2, numel(slst_list_short) );
    for c=1:numel(my_lineages)
        next_lineage_strain = clusters_all_slst_super( my_lineages(c) );
        next_lineage_strain_index = find( next_lineage_strain == slst_list_short );
        my_table_strains_faceback(1,next_lineage_strain_index) = my_table_strains_faceback(1,next_lineage_strain_index) + sum( my_samples_lineages == my_lineages(c) & ismember(my_samples_zones,[ 1 2 3 4 5 6 ]) ); % face
        my_table_strains_faceback(2,next_lineage_strain_index) = my_table_strains_faceback(2,next_lineage_strain_index) + sum( my_samples_lineages == my_lineages(c) & ismember(my_samples_zones,[ 7 ]) ); % back
    end

    if min(sum(my_table_strains_faceback,2)) < min_num_cols_on_zone % fewer than XX colonies for either face or back
        continue
    end

    % Bar chart: relative abundance of lineages on face or back
    fig=figure(1);
    clf(1)
    fs = 32;
    hold on
    box on
    b = bar(my_table_strains_faceback./sum(my_table_strains_faceback,2), 'stacked', 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 2 );
    for k = 1:numel(slst_list_short)
        b(k).CData = slst_list_colors_short(k,:);
    end
    % X
    xlabel('body site')
    xticks( 1:1:2 )
    xticklabels( {'face', 'back'} )
    % Y
    ylabel('relative abundance')
    % Title
    title([ 'subject ' char(my_subject+64) ])
    % Formatting
    set(gca,'FontSize',fs)
    hold off

    % Save fig
    print(fig,[ '3_straintype_faceback/' 'Bar_FaceBack_Subj-' char(my_subject+64) '_rel.png' ],'-dpng')
    
    % Bar chart: absolute abundance of lineages on face or back
    fig=figure(2);
    clf(2)
    fs = 32;
    hold on
    box on
    b = bar(my_table_strains_faceback, 'stacked', 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 2 );
    for k = 1:numel(slst_list_short)
        b(k).CData = slst_list_colors_short(k,:);
    end
    % X
    xlabel('body site')
    xticks( 1:1:2 )
    xticklabels( {'face', 'back'} )
    % Y
    ylabel('number of colonies')
    % Title
    title([ 'subject ' char(my_subject+64) ])
    % Formatting
    set(gca,'FontSize',fs)
    hold off

    % Save fig
    print(fig,[ '3_straintype_faceback/' 'Bar_FaceBack_Subj-' char(my_subject+64) '_abs.png' ],'-dpng')
    

    %% P value (binomial test) for face vs back

    p_face = sum(my_table_strains_faceback(1,:))/sum(sum(my_table_strains_faceback)); % prob colony from face
    p_back = 1-p_face; % prob colony from back

    strains_present = find( sum( my_table_strains_faceback ) > 0 );
    num_strains = numel( strains_present );
    p_value_cutoff_100 = 0.01/num_strains; % Bonferroni: divide by number of strains on this subject
    p_value_cutoff_1000 = 0.001/num_strains; % Bonferroni: divide by number of strains on this subject

    % Probability that face vs non-face counts are higher than observed given p & N
    format long
    fprintf(1, [ 'Subject ' char(my_subject+64) '\n' ] )
    for s=1:num_strains
        this_strain = strains_present(s);
        this_strain_num_face = my_table_strains_faceback(1,this_strain);
        this_strain_num_back = my_table_strains_faceback(2,this_strain);
        p_value_backenrich = binocdf( this_strain_num_back-1, this_strain_num_face+this_strain_num_back, p_back, 'upper' );
        p_value_faceenrich = binocdf( this_strain_num_face-1, this_strain_num_face+this_strain_num_back, p_face, 'upper' );
        if p_value_backenrich < p_value_cutoff_100
            if p_value_backenrich < p_value_cutoff_1000
                fprintf( 1, [ 'Strain ' slst_list_short(this_strain) ' enriched on the back. p=' num2str(p_value_backenrich) ' (**).\n' ] )
            else
                fprintf( 1, [ 'Strain ' slst_list_short(this_strain) ' enriched on the back. p=' num2str(p_value_backenrich) ' (*).\n' ] )
            end
        end
        if p_value_faceenrich < p_value_cutoff_100
            if p_value_faceenrich < p_value_cutoff_1000
                fprintf( 1, [ 'Strain ' slst_list_short(this_strain) ' enriched on the face. p=' num2str(p_value_faceenrich) ' (**).\n' ] )
            else
                fprintf( 1, [ 'Strain ' slst_list_short(this_strain) ' enriched on the face. p=' num2str(p_value_faceenrich) ' (*).\n' ] )
            end
        end
    end
    % Note: imperfect stats since colonies aren't "independent"
    
    
end

% Subject A
% Strain A enriched on the back. p=4.128e-29 (**).
% Strain C enriched on the face. p=2.4024e-05 (**).
% Strain D enriched on the face. p=2.392e-07 (**).
% Subject B
% Strain H enriched on the back. p=0.000314 (*).



