%% MAIN SCRIPT FOR PARALLEL EVOLUTION ANALYSIS

% Last updated: Arolyn, 2021.04.06


%% SETUP

%% Set up directories and environment

% Directory for Lieberman Lab scripts:
dir_scripts_liebermanlab = '../lab_scripts';
path(dir_scripts_liebermanlab,path);

% Add my scripts
dir_my_scripts = [ pwd '/' 'scripts/myscripts_evo'];
path(dir_my_scripts,path);

% Where to find the reference genome on the Lieberman Lab Dropbox:
dir_ref_genome = [ pwd '/reference_genomes/Pacnes_C1'];
[~, GenomeLength, ~, ~] = genomestats(dir_ref_genome);


%%

%%%%%%%%%%%%%%%
%% LOAD DATA %%
%%%%%%%%%%%%%%%


%% Load data about all within-lineage mutations

% Load lineage names
dir_clusters = '2_snvs';
load( [ dir_clusters '/' 'cluster_names.mat' ],'cluster_names','clusters_all') % lineage names
lineage_names = cluster_names; clear cluster_names; 
cluster_sizes = cellfun(@(x) numel(x), clusters_all);

% List of path to files containing annotation_full for each lineage
paths_to_files = cellfun(@(lineage_name) [ dir_clusters '/' lineage_name '/' 'data_' lineage_name '.mat' ], lineage_names, 'UniformOutput', false );
file_tags = lineage_names;

% Get structure with all mutations
annotation_full_everything_all = concatenate_annotation_fulls( paths_to_files, file_tags );

% Number of mutations
% numel(annotation_full_everything_all) % 2813 all muts
% sum(round([annotation_full_everything_all.gene_num])==[annotation_full_everything_all.gene_num]) % 2445 coding regions

% Other nucS mutations?
% sum([annotation_full_everything_all.gene_num]==1223)


%% Load observed mutation spectrum

% Pick mutation spectrum based on observed SNPs
load([ '5_parallel_evo/mutspec/mut_spec_all_lineages.mat' ]) % mut_spec_prob; mut_spec_total; mut_spec_names={ 'AT/TA', 'AC/TG', 'AG/TC', 'GC/CG', 'GT/CA', 'GA/CT' }
% from mutation_spectrum_module via div_matrix2_6types 
% AT, TA % transversion
% AC, TG % transversion
% AG, TC % transition
% GC, CG % transversion
% GT, CA % transversion
% GA, CT % transition


%% Load positions masked in SNP calling

load( [ 'data/plasmid_pos_to_mask.mat' ], 'plasmid_pos_on_genome' )
positions_to_mask = plasmid_pos_on_genome;



%% GENERAL PARAMETERS FOR SUPP FIGS

% Find number of genes under parallel evolution detected at these false
% discovery rates:
fdr_to_test = 0:0.1:1;

% Compute dN/dS for this many genes with the lowest Poisson p-values
top_num_genes_to_test = [ 5, 10 ];

% Minimum number of mutations to evaluate set
min_num_muts = 100;

dir_supp = '5_parallel_evo/supp_figs';
if ~exist( dir_supp, 'dir' )
    mkdir( dir_supp )
end



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARALLEL EVO ANALYSIS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ANALYZE ALL MUTATIONS

% Everything
data_set_name = 'all-muts';
annotation_full_everything = annotation_full_everything_all;

% Call parallel evolution function
dir_save = [ '5_parallel_evo/figs_' data_set_name ];
[ dNdS_all, CIs_all, num_genes_detected_all, dNdS_top_genes_all, CIs_top_genes_all ] = evo_parallel_main( ...
    data_set_name, dir_save, annotation_full_everything, ...
    dir_ref_genome, GenomeLength, positions_to_mask, mut_spec_prob, ...
    fdr_to_test, top_num_genes_to_test );


%% Plots

for n=1:numel(top_num_genes_to_test)
    dNdS_for_plot = dNdS_top_genes_all(n);
    CI_for_plot = CIs_top_genes_all(n,:);
    plot_title = 'all within-lineage SNVs';
    x_axis_label = '';
    x_tick_labels = {['top ' num2str(top_num_genes_to_test(n))]};
    y_axis_factor = 4;
    dir_save = dir_supp;
    plot_file_name = [ 'supp_dNdS_all_top-' num2str(top_num_genes_to_test(n)) ];
    plot_dnds_bar_supp( dNdS_for_plot, CI_for_plot, ...
        plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
        dir_save, plot_file_name ) 
end
%%
dNdS_for_plot = [ dNdS_all(1); dNdS_top_genes_all ];
CI_for_plot = [ CIs_all(1,:); CIs_top_genes_all ];
plot_title = 'all within-lineage SNVs';
x_axis_label = 'genes';
x_tick_labels = { 'all', 'top 10', 'top 5'};
y_axis_factor = 4;
dir_save = dir_supp;
plot_file_name = 'supp_dNdS_all_top+all';
plot_dnds_bar_supp_all( dNdS_for_plot, CI_for_plot', ...
        plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
        dir_save, plot_file_name ) 
%%
dNdS_for_plot = dNdS_top_genes_all;
CI_for_plot = CIs_top_genes_all;
plot_title = 'within-lineage SNVs';
x_axis_label = 'genes';
x_tick_labels = {'top 10', 'top 5'};
y_axis_factor = 4;
dir_save = dir_supp;
plot_file_name = 'supp_dNdS_all_top-all';
plot_dnds_bar_supp_all( dNdS_for_plot, CI_for_plot', ...
        plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
        dir_save, plot_file_name ) 
%%
plot_file_name = 'supp_numgenes_all';
series_labels = {'all genes'};
plot_numgenesdetected_supp( 100*fdr_to_test, num_genes_detected_all, series_labels, plot_title, ...
    dir_save, plot_file_name )
%% MAIN FIG v1
data_set_name = 'all-muts';
load( [ '5_parallel_evo/figs_' data_set_name '/data_dNdS_by_density.mat' ])
dNdS_for_plot = [ dNdS_by_geneset(1:end-1), 1, dNdS_by_geneset(end) ] ;
%CI_for_plot = [ dNdS_lower_by_geneset; dNdS_upper_by_geneset ]';
CI_for_plot = [ [ dNdS_lower_by_geneset(1:end-1), 1, dNdS_lower_by_geneset(end) ]; [ dNdS_upper_by_geneset(1:end-1), 1, dNdS_upper_by_geneset(end) ] ]';
plot_title = 'all lineages';
x_axis_label = 'mutations per kb on gene             ';
x_tick_labels = arrayfun(@(x) { [ '\geq' num2str(x) ] }, list_density_thresholds ); x_tick_labels{end+1} = ''; x_tick_labels{end+1} = 'all';
y_axis_factor = 4;
dir_save = dir_supp;
plot_file_name = 'supp_dNdS_all_by-density';
plot_dnds_bar_supp_all_den( dNdS_for_plot, CI_for_plot', ...
        plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
        dir_save, plot_file_name ) 


%% ANALYZE MUTATIONS BY LINEAGE

dNdS_lineages = cell( numel(lineage_names),1 ); % initialize
CIs_lineages = cell( numel(lineage_names),1 ); % initialize
num_genes_detected_lineages = zeros( numel(lineage_names), numel(fdr_to_test) ); % initialize
dNdS_top_genes_lineages = zeros( numel(lineage_names), numel(top_num_genes_to_test) ); % initialize
CIs_top_genes_lineages = zeros( numel(lineage_names), numel(top_num_genes_to_test), 2 ); % initialize

for lin_index=1:numel(lineage_names)
    
    this_lineage = lineage_names{lin_index};
    annotation_full_everything = annotation_full_everything_all( [annotation_full_everything_all.data_index] == lin_index );
    data_set_name = this_lineage;

    if length( annotation_full_everything ) >= min_num_muts
        % Call parallel evolution function
        dir_save = [ '5_parallel_evo/figs_' data_set_name ];
        [ dNdS_lineage, CIs_lineage, num_genes_detected_lineage, dNdS_top_genes_lineage, CIs_top_genes_lineage ] ...
            = evo_parallel_main( data_set_name, dir_save, annotation_full_everything, ...
            dir_ref_genome, GenomeLength, positions_to_mask, mut_spec_prob, ...
            fdr_to_test, top_num_genes_to_test );
        dNdS_lineages{lin_index} = dNdS_lineage;
        CIs_lineages{lin_index} = CIs_lineage;
        num_genes_detected_lineages(lin_index,:) = num_genes_detected_lineage;
        dNdS_top_genes_lineages(lin_index,:) = dNdS_top_genes_lineage;
        CIs_top_genes_lineages(lin_index,:,:) = CIs_top_genes_lineage;
    else
        fprintf(1,['Insufficient mutations in Lineage ' this_lineage ' (n_muts=' num2str(length(annotation_full_everything)) ')' '\n'])
        num_genes_detected_lineages(lin_index,:) = nan;
        dNdS_top_genes_lineages(lin_index,:) = nan;
        CIs_top_genes_lineages(lin_index,:,:) = nan;
    end
    
end

%% Plot

dNdS_for_plot = dNdS_top_genes_lineages( ~isnan(dNdS_top_genes_lineages(:,1)),: );
CI_for_plot = CIs_top_genes_lineages( ~isnan(dNdS_top_genes_lineages(:,1)),:,: );
plot_title = 'by lineage';
x_axis_label = 'lineage';
x_tick_labels = find( ~isnan(dNdS_top_genes_lineages(:,1)) );
x_tick_labels = {'A-1','A-2','A-3','H-1','F-1','B-1'}; % manually; first six lineages
y_axis_factor = 4;
dir_save = dir_supp;
plot_file_name = 'supp_dNdS_by-lineage';
plot_dnds_bar_supp_multi( dNdS_for_plot, CI_for_plot, top_num_genes_to_test, ...
    plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
    dir_save, plot_file_name ) 
%%
plot_file_name = 'supp_numgenes_by-lineage';
series_labels = x_tick_labels;
num_genes_for_plot = num_genes_detected_lineages( ~isnan(dNdS_top_genes_lineages(:,1)),: );
plot_numgenesdetected_supp( 100*fdr_to_test, num_genes_for_plot, series_labels, plot_title, ...
    dir_save, plot_file_name )



%% ANALYZE MUTATIONS BY SUBJECT

subjects_list = unique( [annotation_full_everything_all.subject] );

dNdS_subjects = cell( numel(subjects_list),1 ); % initialize
CIs_subjects = cell( numel(subjects_list),1 ); % initialize
num_genes_detected_subjects = zeros( numel(subjects_list), numel(fdr_to_test) ); % initialize
dNdS_top_genes_subjects = zeros( numel(subjects_list), numel(top_num_genes_to_test), 1 ); % initialize
CIs_top_genes_subjects = zeros( numel(subjects_list), numel(top_num_genes_to_test), 2 ); % initialize

for lin_index=1:numel(subjects_list)
    
    this_subject = subjects_list(lin_index);
    annotation_full_everything = annotation_full_everything_all( [annotation_full_everything_all.subject] == this_subject );
    data_set_name = [ 'subj-' this_subject ];

    if length( annotation_full_everything ) >= min_num_muts
        % Call parallel evolution function
        dir_save = [ '5_parallel_evo/figs_' data_set_name ];
        [ dNdS_subject, CIs_subject, num_genes_detected_subject, dNdS_top_genes_subject, CIs_top_genes_subject ] ...
            = evo_parallel_main( data_set_name, dir_save, annotation_full_everything, ...
            dir_ref_genome, GenomeLength, positions_to_mask, mut_spec_prob, ...
            fdr_to_test, top_num_genes_to_test );
        dNdS_subjects{lin_index} = dNdS_subject;
        CIs_subjects{lin_index} = CIs_subject;
        num_genes_detected_subjects(lin_index,:) = num_genes_detected_subject;
        dNdS_top_genes_subjects(lin_index,:) = dNdS_top_genes_subject;
        CIs_top_genes_subjects(lin_index,:,:) = CIs_top_genes_subject;
    else
        fprintf(1,['Insufficient mutations in Subject ' this_subject ' (n_muts=' num2str(length(annotation_full_everything)) ')' '\n'])
        num_genes_detected_subjects(lin_index,:) = nan;
        dNdS_top_genes_subjects(lin_index,:) = nan;
        CIs_top_genes_subjects(lin_index,:,:) = nan;
    end
    
end


%%

dNdS_for_plot = dNdS_top_genes_subjects( ~isnan(dNdS_top_genes_subjects(:,1)),: );
CI_for_plot = CIs_top_genes_subjects( ~isnan(CIs_top_genes_subjects(:,1)),:,: );
plot_title = 'by subject';
x_axis_label = 'subject';
x_tick_labels = subjects_list( ~isnan(dNdS_top_genes_subjects(:,1)) ); x_tick_labels = arrayfun(@(x) {x}, x_tick_labels);
y_axis_factor = 4;
dir_save = dir_supp;
plot_file_name = 'supp_dNdS_by-subject';
plot_dnds_bar_supp_multi( dNdS_for_plot, CI_for_plot, top_num_genes_to_test, ...
    plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
    dir_save, plot_file_name ) 
%%
plot_file_name = 'supp_numgenes_by-subject';
series_labels = x_tick_labels;
num_genes_for_plot = num_genes_detected_subjects( ~isnan(num_genes_detected_subjects(:,1)),: );
plot_numgenesdetected_supp( 100*fdr_to_test, num_genes_for_plot, series_labels, plot_title, ...
    dir_save, plot_file_name )


%% ANALYZE MUTATIONS INFERRED TO HAVE OCCURRED INSIDE PORES

clusters_list = 1:1:numel(lineage_names); 

% Generate annotation_full 
for c=1:numel(clusters_list)

    % Next clade name
    lineage_next = lineage_names{clusters_list(c)};
    % Check if clade has enough SNPs for annotations; if so, load and append
    goodpos = load([ dir_clusters '/' lineage_next '/' 'data_' lineage_next '.mat' ], 'goodpos' );
    goodpos = goodpos.goodpos;
    
    if numel(goodpos) > 0

        % Split SNPs by within-pore vs not
        [ annotation_pores, annotation_other ] = ...
            filter_snps_by_pore( lineage_next, dir_clusters );
        
        % Fill in annotation_pores_all and annotation_other_all with SNPs from this lineage
        if c==1
            % Intialize annotation_pores_all and annotation_other_all
            annotation_pores_all = annotation_pores;
            annotation_other_all = annotation_other;
            fields_list = fieldnames(annotation_pores_all); % get all fields
            for i=1:length(annotation_pores_all)
                annotation_pores_all(i).lineage = lineage_next;
            end
            for i=1:length(annotation_other_all)
                annotation_other_all(i).lineage = lineage_next;
            end
        else
            % Add rows to annotation_pores_all
            old_length = length(annotation_pores_all);
            fields_next = fieldnames(annotation_pores);
            for snp=1:length(annotation_pores) % loop through eacn new SNP
                for f=1:length(fields_next) % update each field
                    annotation_pores_all(old_length+snp).(fields_next{f}) = annotation_pores(snp).(fields_next{f});
                    annotation_pores_all(old_length+snp).lineage = lineage_next;
                end
            end
            % Add rows to annotation_other_all
            old_length = length(annotation_other_all);
            fields_next = fieldnames(annotation_other);
            for snp=1:length(annotation_other) % loop through eacn new SNP
                for f=1:length(fields_next) % update each field
                    annotation_other_all(old_length+snp).(fields_next{f}) = annotation_other(snp).(fields_next{f});
                    annotation_other_all(old_length+snp).lineage = lineage_next;
                end
            end
        end
        
    end
end


%% Save pore data

% data_set_name = 'data_within-pore';
% if ~exist( data_set_name, 'dir' )
%     mkdir( data_set_name )
% end
% save( [ data_set_name '/' 'data_pores.mat' ], 'annotation_pores_all', 'annotation_other_all' )


%%

% Call parallel evolution function
top_num_genes_to_test_pores = [3];
data_set_name = 'within-pore';
dir_save = [ '5_parallel_evo/figs_' data_set_name ];
[ dNdS_withinpore, CIs_withinpnore, num_genes_detected_withinpore, dNdS_top_genes_withinpore, CIs_top_genes_withinpore ] = ...
    evo_parallel_main( data_set_name, dir_save, annotation_pores_all, ...
    dir_ref_genome, GenomeLength, positions_to_mask, mut_spec_prob, ...
    fdr_to_test, top_num_genes_to_test_pores );

%%
dNdS_for_plot = [ dNdS_withinpore(1), dNdS_top_genes_withinpore ];
CI_for_plot = [ CIs_withinpnore(1,:); CIs_top_genes_withinpore ]';
plot_title = 'within pore SNVs';
x_axis_label = 'number of mutations on gene';
x_tick_labels = { 'm>0', 'm>1' };
y_axis_factor = 4;
dir_save = dir_supp;
plot_file_name = 'supp_dNdS_pores';
dNdS_for_plot( isinf( dNdS_for_plot ) ) = y_axis_factor;
CI_for_plot( isinf( CI_for_plot ) ) = y_axis_factor;
plot_dnds_bar_supp( dNdS_for_plot, CI_for_plot, ...
    plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
    dir_save, plot_file_name ) 
%%
plot_title = '';
plot_file_name = 'supp_dNdS_pores_notitle';
plot_dnds_bar_supp( dNdS_for_plot, CI_for_plot, ...
    plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
    dir_save, plot_file_name ) 
%%
% save( [ dir_results '/' 'data_old_method.mat' ], ...
%     'filter_min_gene_length', 'min_density_list', 'max_num_muts_per_gene', ...
%     'record_obs_num_genes_with_n_muts', 'record_sim_num_genes_with_n_muts_mean', 'record_sim_num_genes_with_n_muts_std' )
load( [ '5_parallel_evo/figs_' data_set_name '/' 'data_old_method.mat' ] )
plot_title = 'within pore SNVs';
x_axis_label = '\geq m mutations';
plot_file_name = 'supp_doublebar_pores';
make_inset = false; inset_start = 0;
plot_parallelevo_doublebar_supp( record_obs_num_genes_with_n_muts(1,:), ...
    record_sim_num_genes_with_n_muts_mean(1,:), record_sim_num_genes_with_n_muts_std(1,:), ...
    plot_title, x_axis_label, dir_supp, plot_file_name, make_inset, inset_start ) 
%%
plot_title = 'within pore SNVs';
x_axis_label = '\geq m mutations';
plot_file_name = 'supp_doublebar_pores_short';
make_inset = false; inset_start = 0;
plot_parallelevo_doublebar_supp( record_obs_num_genes_with_n_muts(1,:), ...
    record_sim_num_genes_with_n_muts_mean(1,:), record_sim_num_genes_with_n_muts_std(1,:), ...
    plot_title, x_axis_label, dir_supp, plot_file_name, make_inset, inset_start )

%% MAIN FIG v2
data_set_name = 'all-muts';
load( [ '5_parallel_evo/figs_' data_set_name '/data_dNdS_by_density.mat' ])
dNdS_for_plot = [ dNdS_by_geneset(end), 1, dNdS_by_geneset(1:end-1), 1, dNdS_withinpore(1) ] ;
%CI_for_plot = [ dNdS_lower_by_geneset; dNdS_upper_by_geneset ]';
CI_for_plot = [ [ dNdS_lower_by_geneset(end), 1, dNdS_lower_by_geneset(1:end-1), 1, CIs_withinpnore(1,1) ]; [ dNdS_upper_by_geneset(end), 1, dNdS_upper_by_geneset(1:end-1), 1, CIs_withinpnore(1,2) ] ]';
plot_title = 'all lineages';
x_axis_label = 'mutations per kb on gene';
x_tick_labels = arrayfun(@(x) { [ '\geq' num2str(x) ] }, list_density_thresholds ); x_tick_labels{end+1} = ''; x_tick_labels{end+1} = 'intrapore'; x_tick_labels = horzcat( {'all', ''}, x_tick_labels );
y_axis_factor = 2;
dir_save = dir_supp;
plot_file_name = 'supp_dNdS_all_by-density_plus-pores';
plot_dnds_bar_supp_all_den_w_pores( dNdS_for_plot, CI_for_plot', ...
    plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
    dir_save, plot_file_name ) 


%%

%%%%%%%%%%%%%%%%%%%%
%% KEGG FUNCTIONS %%
%%%%%%%%%%%%%%%%%%%%


%% Load KEGG functions

% Load KEGG functional info (provided by Alex Poret)
load( [ 'data' '/' 'functional_analysis_revised.mat' ] )
% functional_analysis_structure
% overarching_cat
% sub_cat


%% Grab gene numbers for KEGG functions

% Coarse KEGG categories
kegg_cat_coarse_num = length(overarching_cat);
kegg_cat_coarse_names = cell( kegg_cat_coarse_num,1 ); % initialize
kegg_cat_coarse_genenums = cell( kegg_cat_coarse_num,1 ); % initialize
for i=1:kegg_cat_coarse_num
    kegg_cat_coarse_names{i} = overarching_cat{i,1};
    temp_function_indices = overarching_cat{i,2};
    temp_genenums = [];
    for j=1:numel(temp_function_indices)
        temp_genenums = [ temp_genenums, functional_analysis_structure{ temp_function_indices(j),5 } ];
    end
    kegg_cat_coarse_genenums{i} = temp_genenums( temp_genenums > 0 ); % remove nans for KEGG genes that aren't in this reference genome
end


%% NEUTRAL MODEL %%


%% Compute neutral model probability of nonsynonymous mutation for each coarse KEGG category

% Compute neutral model for genes in each KEGG category
dnds_expected_kegg_coarse = zeros( kegg_cat_coarse_num,1 ); % initialize

for i=1:kegg_cat_coarse_num

    % Compute
    probnonsyn_expected = compute_expected_dnds( dir_ref_genome, mut_spec_prob, kegg_cat_coarse_genenums{i} );
    
    % Save
    dnds_expected_kegg_coarse(i) = probnonsyn_expected/(1-probnonsyn_expected); % N/S expected from neutral model
    
end


%% OBSERVED MUTATIONS %%


%% Compute observed dN/dS for each coarse KEGG category

% Compute observed N/S and dN/dS for genes in each KEGG category
dnds_observed_kegg_coarse = zeros( kegg_cat_coarse_num,1 ); % initialize
dNdS_kegg_coarse = zeros( kegg_cat_coarse_num,1 ); % initialize
CI_kegg_coarse_upper = zeros( kegg_cat_coarse_num,1 ); % initialize
CI_kegg_coarse_lower = zeros( kegg_cat_coarse_num,1 ); % initialize
nummuts_observed_kegg_coarse = zeros( kegg_cat_coarse_num,1 ); % initialize

for i=1:kegg_cat_coarse_num

    % Compute observed N/S
    [ p_nonsyn, CI_nonsyn, num_muts_N, num_muts_S ] = compute_observed_dnds( annotation_full_everything_all, kegg_cat_coarse_genenums{i} );
    dnds_observed = (p_nonsyn/(1-p_nonsyn)); % N/S observed
    dNdS = dnds_observed/dnds_expected_kegg_coarse(i); % relative to neutral model for this KEGG category
    CI_lower = ( CI_nonsyn(1)/(1-CI_nonsyn(1)) ) / dnds_expected_kegg_coarse(i); 
    CI_upper = ( CI_nonsyn(2)/(1-CI_nonsyn(2)) ) / dnds_expected_kegg_coarse(i);

    % Save
    dnds_observed_kegg_coarse(i) = dnds_observed;
    dNdS_kegg_coarse(i) = dNdS; % initialize
    CI_kegg_coarse_lower(i) = CI_lower; % initialize
    CI_kegg_coarse_upper(i) = CI_upper; % initialize
    nummuts_observed_kegg_coarse(i) = num_muts_N + num_muts_S;

end

%% Figure

% Make a figure
dNdS_for_plot = dNdS_kegg_coarse;
CI_for_plot = [ CI_kegg_coarse_lower, CI_kegg_coarse_upper ]';
plot_title = [ 'kegg-coarse' ];
x_axis_label = 'coarse KEGG function';
x_tick_labels = kegg_cat_coarse_names;
dir_save = '5_parallel_evo/figs_dnds_kegg'; 
plot_file_name = [ 'bar_dnds_' plot_title ];
y_axis_factor = 4;
text_labels = nummuts_observed_kegg_coarse;
% Call plot function
plot_dnds_bar( dNdS_for_plot, CI_for_plot, plot_title, x_axis_label, x_tick_labels, y_axis_factor, dir_save, plot_file_name, text_labels ) 


%% Supp figure

dNdS_for_plot = dNdS_kegg_coarse;
CI_for_plot = [ CI_kegg_coarse_lower, CI_kegg_coarse_upper ]';
plot_title = 'by function';
x_axis_label = 'KEGG function';
x_tick_labels = kegg_cat_coarse_names;
y_axis_factor = 4;
dir_save = 'supp_figs';
plot_file_name = 'supp_dNdS_kegg-coarse';
plot_dnds_bar_supp( dNdS_for_plot, CI_for_plot, ...
    plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
    dir_save, plot_file_name ) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLUSTER 2 MUTATIONAL AGE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compute neutral model probability of nonsynonymous mutation

cdsindex_proteins = get_coding_genenums( dir_ref_genome, positions_to_mask ); % avoid genes in regions masked in SNP calling

probnonsyn_expected = compute_expected_dnds( dir_ref_genome, mut_spec_prob, cdsindex_proteins );
dnds_expected = probnonsyn_expected/(1-probnonsyn_expected); % N/S expected from neutral model


%% OBSERVED MUTATIONS %%


%% Load lineage A-2 SNVs

% Which lineage
this_lineage = 2; 
this_lineage_name = lineage_names{ this_lineage };

load([ dir_clusters '/' this_lineage_name '/' 'data_' this_lineage_name '.mat' ] )

Calls_lineage = Calls_for_analysis(:,~outgroup_isolates); % exclude outgroup isolates
names_lineage = SampleNames(~outgroup_isolates)';
num_samples_in_lineage = sum(~outgroup_isolates);

% Boolean showing which positions differ from the inferred lineage MRCA for each colony in that lineage
diff_mrca = ( ...
    ( Calls_lineage ~= repmat(anc_nti_goodpos,1,num_samples_in_lineage) ) ...
    & Calls_lineage > 0 ...
    );

%% Load lineage A-2 tree info

lin2_table = readtable('data/lineage2_subcladeinfo.csv');
lin2_names = lin2_table.sample_name;
lin2_subclades = lin2_table.sub_clade_num;
num_subclades = max( lin2_subclades );

subclade_indices = arrayfun(@(x) ismember( names_lineage, lin2_names( x==lin2_subclades ) ), 1:1:num_subclades, 'UniformOutput', false );


%% Determine which mutations are on long branches

is_on_long_branch = zeros( size( goodpos ), 'logical' );
for i=1:numel(goodpos)
    mut_var_in_subclade = cell2mat( cellfun(@(x) numel( unique( diff_mrca(i,x) ) ) > 1, subclade_indices, 'UniformOutput', false ) );
    is_on_long_branch(i) = ( sum( mut_var_in_subclade ) == 0 );
end


%% Evaluate mutation types and get gene numbers

muts_genenums = [ annotation_full.gene_num ];
genenums_old = [ annotation_full(is_on_long_branch).gene_num ];
genenums_new = [ annotation_full(~is_on_long_branch).gene_num ];

types_old = [annotation_full(is_on_long_branch).type]; % long branches
types_new = [annotation_full(~is_on_long_branch).type]; % within subclades

genenums_old_coding = genenums_old( types_old == 'N' | types_old == 'S' );
genenums_new_coding = genenums_new( types_new == 'N' | types_new == 'S' );


%% Compute dN/dS for old mutations

% Compute observed N/S
next_num_N = sum( types_old=='N' );
next_num_S = sum( types_old=='S' );
x = next_num_N;
n = next_num_N + next_num_S;
alpha = 0.05;
[ p_nonsyn, CI_nonsyn ] = binofit( x, n, alpha );
dnds_observed = (p_nonsyn/(1-p_nonsyn)); % N/S observed

% Compute expected N/S
next_genenums = unique( genenums_old_coding );
probnonsyn_expected = compute_expected_dnds( dir_ref_genome, mut_spec_prob, next_genenums );
dnds_expected_local = probnonsyn_expected/(1-probnonsyn_expected); % N/S expected from neutral model

% Compute dN/dS
dNdS_old = dnds_observed/dnds_expected_local; % relative to neutral model
CI_lower_old = ( CI_nonsyn(1)/(1-CI_nonsyn(1)) ) / dnds_expected_local;
CI_upper_old = ( CI_nonsyn(2)/(1-CI_nonsyn(2)) ) / dnds_expected_local;



%% Compute dN/dS for new mutations

% Compute observed N/S
next_num_N = sum( types_new=='N' );
next_num_S = sum( types_new=='S' );
x = next_num_N;
n = next_num_N + next_num_S;
alpha = 0.05;
[ p_nonsyn, CI_nonsyn ] = binofit( x, n, alpha );
dnds_observed = (p_nonsyn/(1-p_nonsyn)); % N/S observed

% Compute expected N/S
next_genenums = unique( genenums_new_coding );
probnonsyn_expected = compute_expected_dnds( dir_ref_genome, mut_spec_prob, next_genenums );
dnds_expected_local = probnonsyn_expected/(1-probnonsyn_expected); % N/S expected from neutral model

% Compute dN/dS
dNdS_new = dnds_observed/dnds_expected_local; % relative to neutral model
CI_lower_new = ( CI_nonsyn(1)/(1-CI_nonsyn(1)) ) / dnds_expected_local;
CI_upper_new = ( CI_nonsyn(2)/(1-CI_nonsyn(2)) ) / dnds_expected_local;


%% Plot

dNdS_for_plot = [ dNdS_old, dNdS_new ];
CI_for_plot = [ CI_lower_old, CI_lower_new ; CI_upper_old, CI_upper_new ];
plot_title = 'lineage A-2';
x_axis_label = 'mutation age';
%x_tick_labels = {'between subclade','within subclade'};
x_tick_labels = {'old','new'};
y_axis_factor = 4;
dir_save = dir_supp;
plot_file_name = 'supp_dNdS_lineage2mrca';
plot_dnds_bar_supp( dNdS_for_plot, CI_for_plot, ...
    plot_title, x_axis_label, x_tick_labels, y_axis_factor, ...
    dir_save, plot_file_name ) 



%% McDonald-Kreitman Test

num_muts_N_old = sum( types_old=='N' );
num_muts_N_new = sum( types_new=='N' );
num_muts_S_old = sum( types_old=='S' );
num_muts_S_new = sum( types_new=='S' );

odds_MK = (num_muts_N_old/num_muts_S_old)/(num_muts_N_new/num_muts_S_new)
table_MK = table([num_muts_N_new;num_muts_S_new],[num_muts_N_old;num_muts_S_old],'VariableNames',{'New','Old'},'RowNames',{'N','S'})
[h,p,stats] = fishertest(table_MK,'Alpha',0.05)


% OLD

% odds_MK =
% 
%     0.8157
% 
% 
% table_MK =
% 
%   2×2 table
% 
%          New    Old
%          ___    ___
% 
%     N    136    50 
%     S     71    32 
% 
% 
% h =
% 
%   logical
% 
%    0
% 
% 
% p =
% 
%     0.4964
% 
% 
% stats = 
% 
%   struct with fields:
% 
%              OddsRatio: 1.2259
%     ConfidenceInterval: [0.7227 2.0795]