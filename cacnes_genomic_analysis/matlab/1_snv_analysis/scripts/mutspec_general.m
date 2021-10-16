%% SUMMARY

% This script computes the observed mutation spectrum.


%% Version History

% Arolyn, 2021.02.03: combined old scripts that did various things relating
% to the mutation spectrum



%% SETUP


%% Directory setup

% Directory for Lieberman Lab scripts:
dir_scripts_liebermanlab = '../lab_scripts';
path(dir_scripts_liebermanlab,path);

% Add my scripts
dir_my_scripts = [ pwd '/' 'scripts/myscripts_evo'];
path(dir_my_scripts,path);

% Where to find cluster data
dir_lineages = '2_snvs';

% Names of lineages
load( [ dir_lineages '/' 'cluster_names' ] )
num_lineages = numel(cluster_names);


%% Basic info

% Nucleotides: 1=A, 2=T, 3=C, 4=G
NTs = 'ATCG';



%% COMPUTE MUTATION SPECTRUM FOR EACH LINEAGE

% Directory to store info
dir_mut_spec = '5_parallel_evo/mutspec';
if ~exist([dir_mut_spec],'dir')
    mkdir([dir_mut_spec])
end

% Initialize
lineage_slst = cell( num_lineages,1 ); % save SLST
lineage_mut_specs = zeros( num_lineages,6 ); % save mutation spectrum 

% Loop through lineages
for c=1:num_lineages

    next_lineage = cluster_names{c};
    fprintf(1, [ 'Analyzing mutation spectrum for lineage ' next_lineage '...\n' ] )

    % Load data
    load([ dir_lineages '/' next_lineage '/' 'data_' next_lineage '.mat' ])
    lineage_slst{c} = unique(slst(~outgroup_isolates));

    %% Remove hypermutator mutations in lineage 1
    if contains( next_lineage, 'Cluster-1-A' )
        fprintf(1,'Removing hypermutator SNPs from lineage 1...\n')
        is_hyper = ( specimen_numbers==31 | specimen_numbers==41 );
        calls_hypers=Calls_for_analysis(:,is_hyper);
        diff_mrca=((calls_hypers~=repmat(anc_nti_goodpos,1,sum(is_hyper))) & calls_hypers>0);
        hyper_muts = sum(diff_mrca,2)>0;
        fprintf(1,['Ignoring ' num2str(sum(hyper_muts)) ' mutations...\n'])
        goodpos = goodpos( ~hyper_muts );
    end

    % Analyze mutations
    save_data = 1;
    [ mutationmatrix, mut_observed, typecounts, prob_nonsyn ] = ...
        mutation_spectrum_module( goodpos, anc_nti_goodpos, Calls_for_analysis, annotation_full, p, ...
        next_lineage, dir_mut_spec, save_data );

    % Record
    lineage_mut_specs(c,:) = mut_observed;

end



%% COMPUTE AGGREGATE MUTATION SPECTRUM ACROSS ALL LINEAGES

mut_spec_total = sum( lineage_mut_specs );
mut_spec_prob = mut_spec_total/sum(mut_spec_total);

% mut_observed
%AT, TA % transversion
%AC, TG % transversion
%AG, TC % transition
%GC, CG % transversion
%GT, CA % transversion
%GA, CT % transition
% from mutation_spectrum_module via div_matrix2_6types 
mut_types_names = { 'AT/TA', 'AC/TG', 'AG/TC', 'GC/CG', 'GT/CA', 'GA/CT' };


%% Bar chart

figure(1)
clf(1)
hold on
box on
color_face = 0.75*[ 1 1 1 ];
bar( mut_spec_total, 'FaceColor', color_face )
xticks(1:1:6)
xticklabels( mut_types_names )
xlabel('mutation type')
ylabel('number of mutations')
title('mutation spectrum (all SNVs)')
set(gca, 'FontSize', 16)
hold off

% Save image
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
print([ dir_mut_spec '/Bar_MutSpec_all.png' ],'-dpng')


%% Save mutation spectrum

save( [ dir_mut_spec '/mut_spec_all_lineages' ], 'mut_spec_total', 'mut_spec_prob', 'mut_types_names')


%% COMPARE MUTATION SPECTRA ACROSS STRAIN TYPES


%% Compute mutation spectra for each strain type

% Get superSLST for each lineage
lineage_slst_super = char(num_lineages);
for c=1:num_lineages
    next_slst_list = lineage_slst{c};
    next_slst_super = unique(cellfun(@(x) x(1), next_slst_list ));
    lineage_slst_super(c) = next_slst_super;
end

% List of unique superSLSTs
slst_super_list = unique( lineage_slst_super );
num_slst_super = numel( slst_super_list );

% Get mutation spectrum for each superSLST
mut_spec_superslst = zeros( num_slst_super,6 );
for s=1:num_slst_super
    next_slst_super = slst_super_list(s);
    if sum(next_slst_super==lineage_slst_super)==1
        next_mut_spec = lineage_mut_specs( next_slst_super==lineage_slst_super,: );
    else
        next_mut_spec = sum( lineage_mut_specs( next_slst_super==lineage_slst_super,: ) );
    end
    mut_spec_superslst(s,:) = next_mut_spec;
end

% Convert to fractions
mut_spec_superslst_prob = mut_spec_superslst./sum(mut_spec_superslst,2);


%% Bar chart

% colormap
% write numbers

figure(2)
clf(2)
hold on
box on
% Heatmap
imagesc( flipud(mut_spec_superslst_prob) )
% Colormap
my_colormap = ones( 100,3 ); my_colormap_scale = fliplr(0.1:(1-0.1)/99:1);
my_colormap = my_colormap_scale.*my_colormap'; my_colormap = my_colormap';
caxis([0 1])
colormap( my_colormap )
colorbar( 'Ticks', 0:0.25:1 )
% Numbers
mut_spec_superslst_flipped = flipud( mut_spec_superslst );
for i=1:numel(slst_super_list)
    for j=1:6
        text( j, i, num2str(mut_spec_superslst_flipped(i,j)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 16 )
    end
end
% Grids
for i = 1:numel(slst_super_list)+1
    lw = 1.25;
    plot([.5,6.5],[i-.5,i-.5],'k-','LineWidth',lw);
    plot([i-.5,i-.5],[.5,numel(slst_super_list)+0.5],'k-','LineWidth',lw);
end
% Labels
xlabel('mutation type')
xticks(1:1:6)
xticklabels( mut_types_names )
xlim( [0.5 6.5] )
ylabel('superSLST')
yticks(1:1:numel(slst_super_list))
yticklabels( fliplr(arrayfun(@(x) {x}, slst_super_list)) )
ylim( [0.5 numel(slst_super_list)+0.5] )
title('mutation spectrum by superSLST')
set(gca, 'FontSize', 16)
hold off

% Save image
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
print([ dir_mut_spec '/Heatmap_MutSpec_superSLST.png'],'-dpng')


%% Compare mutation spectra across strain types

% Chi-squared test across all strain types
pval_all = chi_squared_test( mut_spec_superslst );

% pval_all =
% 
%    2.1394e-08


%%

% Chi-squared test comparing each strain type to across all strain types
pval_strain = zeros( 1,num_slst_super );
for s=1:num_slst_super
    pval_strain(s) = chi_squared_test( [ mut_spec_total; mut_spec_superslst(s,:) ] );
end

% pval_strain =
% 
%     0.2050    0.0031    0.0000    0.2991    0.2973    0.0407    0.5146






