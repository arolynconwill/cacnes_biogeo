%% FIGURE 5 SCRIPT

% For plots with dMRCAs

% NOTE: interpore comparison is interpore mean dMRCA
% i.e. infer ancestor of each pore
% infer ancestor between pair of pores
% take average of distance between ancestor and pore A and distance
% between ancestor and pore B


%% ENVIRONMENT

dir_scripts_aro_1 = [ 'scripts/myscripts_denovomuts' ];
dir_scripts_aro_2 = [ 'scripts/myscripts_biogeo' ];
path(path,dir_scripts_aro_1);
path(path,dir_scripts_aro_2);


%% DATA

% Cluster info
load('2_snvs/cluster_names')
SampleNames_all = load('data/cluster_step_variables','SampleNames_final'); SampleNames_all = SampleNames_all.SampleNames_final;

% Metadata
types_all = load('data/cluster_step_variables','types_final'); types_all = types_all.types_final;
specimen_number_all = load('data/cluster_step_variables','specimen_number_final'); specimen_number_all = specimen_number_all.specimen_number_final;
subjects_all = load('data/cluster_step_variables','subjects_final'); subjects_all = subjects_all.subjects_final;

% Load mutliplicity of pore specimens (strips and extracts)
load('data/spec_mult.mat')


%% GET dMRCAS

%% Calculate metrics for clade by clade

% Filters
min_num_pores = 5; % clade must have more than min_num_pores in order to be analyzed on its own; pores in small clades combined together

% For storing dMRCAs for each clade
clades_cladenums = [];
clades_dmrcas_max = [];
clades_dmrcas_median = [];
clades_dmrcas_median_all = []; % for all clades
clades_intrapore = {};
clades_interpore_pairs = {};
clades_interpore_min = {};
clades_numpores = [];
clades_porespecnums = {};

% For temporarily storing dMRCAs for pores in small clades
temp_dmrca_clade_max = [];
temp_dmrca_intrapore = [];
temp_dmrca_interpore_pairs = [];
temp_dmrca_interpore_min = [];
temp_num_pores = 0;
temp_pore_specnum = [];
temp_pore_cladenum = [];

% Import data and compute
for c=1:numel(clusters_all) % only clusters that have at least 10 samples
    
    this_cluster_name = cluster_names{c};
    fprintf(1,[this_cluster_name '\n'])

    % Compute info for clade
    [ dmrca_clade_max, dmrca_clade_median, ...
        dmrca_intrapore, ...
        dmrca_interpore_pairs_list, dmrca_interpore_min_list, ...
        num_pores, pore_specnum ] = ...
        get_clade_dmrcas( this_cluster_name, spec_mult_specnums, spec_mult_porenums );
    clades_dmrcas_median_all(c) = dmrca_clade_median;
    % excludes hypermutator specimens 31 and 41
    % excludes pore specimens that came from multiple pores
    % returns NaN's if no data for clade
    
    % Save
    if num_pores > min_num_pores
        clades_cladenums(end+1) = c;
        clades_dmrcas_max(end+1) = dmrca_clade_max;
        clades_dmrcas_median(end+1) = dmrca_clade_median;
        clades_intrapore{end+1} = dmrca_intrapore;
        clades_interpore_pairs{end+1} = dmrca_interpore_pairs_list;
        clades_interpore_min{end+1} = dmrca_interpore_min_list;
        clades_numpores(end+1) = num_pores;
        clades_porespecnums{end+1} = pore_specnum;
    elseif num_pores > 1 % need pore pair to have at least one interpore dMRCA
        temp_dmrca_clade_max = [ temp_dmrca_clade_max; dmrca_clade_max ];
        temp_dmrca_intrapore = [ temp_dmrca_intrapore; dmrca_intrapore ];
        temp_dmrca_interpore_pairs = [ temp_dmrca_interpore_pairs; dmrca_interpore_pairs_list ];
        temp_dmrca_interpore_min = [ temp_dmrca_interpore_min; dmrca_interpore_min_list ];
        temp_num_pores = temp_num_pores + num_pores;
        temp_pore_specnum = [ temp_pore_specnum, pore_specnum ];
        temp_pore_cladenum = [ temp_pore_cladenum, c*ones(1,num_pores) ];
    elseif num_pores == 1 % just keep intrapore
        temp_dmrca_clade_max = [ temp_dmrca_clade_max; dmrca_clade_max ];
        temp_dmrca_intrapore = [ temp_dmrca_intrapore; dmrca_intrapore ];
        temp_num_pores = temp_num_pores + num_pores;
        temp_pore_specnum = [ temp_pore_specnum, pore_specnum ];
        temp_pore_cladenum = [ temp_pore_cladenum, c ];
    end
    
end

% Add data from small clades aggregated as clade "0"
if temp_num_pores>1
    clades_cladenums(end+1) = 0;
    clades_dmrcas_max(end+1) = max(temp_dmrca_clade_max);
    clades_intrapore{end+1} = temp_dmrca_intrapore;
    clades_interpore_pairs{end+1} = temp_dmrca_interpore_pairs;
    clades_interpore_min{end+1} = temp_dmrca_interpore_min;
    clades_numpores(end+1) = temp_num_pores;
    clades_porespecnums{end+1} = temp_pore_specnum;
end


%% Get some numbers for intrapore dMRCAs across all pores
% includes all pores that come from a single follicle for which there are
% at least 2 colonies

clades_intrapore_list_all = cell2mat( cellfun(@(x) x', clades_intrapore, 'UniformOutput', false ));
median(clades_intrapore_list_all) % O SNVs
mean(clades_intrapore_list_all) % 1.3 SNVs
quantile( clades_intrapore_list_all, [0.25 0.75] ) % 0 to 0.6 SNVs


%% Get some numbers for intrapore dMRCAs in Lineage A-1

clade_1_intrapore = clades_intrapore{1}; % non hypermutators
median(clade_1_intrapore) % O SNVs
mean(clade_1_intrapore) % 1.4 SNVs
quantile( clade_1_intrapore, [0.25 0.75] ) % 0 to 1.1 SNVs

clade_1_interpore = clades_interpore_pairs{1};
median( clade_1_interpore ) % 6 SNVs
mean( clade_1_interpore ) % 6.2 SNVs
quantile( clade_1_interpore, [0.25 0.75] ) % 4-8.5 SNVs


%% PLOT DMRCAS FOR EACH CLADE

for this_clade = 1:numel(clades_cladenums)
    
    this_clade_num = clades_cladenums(this_clade);

    % Make figure
    figure(1)
    clf(1)
    % Intrapore
    subplot(1,3,1)
    boxplot(clades_intrapore{this_clade}, ...
        'Notch','on', 'OutlierSize',10, 'Colors','kkkk', 'Symbol','k+')
    xlim([.85 1.15])
    xticks([1])
    xticklabels({'Intrapore'})
    ylim([ -0.5 clades_dmrcas_max(this_clade) ])
    ylabel('Mean dMRCA')
    set(gca,'FontSize',20)
    % Interpore: pairwise
    subplot(1,3,2)
    boxplot(clades_interpore_pairs{this_clade}, ...
            'Notch','on', 'OutlierSize',10, 'Colors','kkkk', 'Symbol','k+')
    xlim([.8 1.2])
    xticks([1])
    xticklabels({'Interpore'})
    ylim([ -0.5 clades_dmrcas_max(this_clade) ])
    ylabel('Min dMRCA')
    set(gca,'FontSize',20)
    % Interpore: minimum
    subplot(1,3,3)
    boxplot(clades_interpore_min{this_clade}, ...
            'Notch','on', 'OutlierSize',10, 'Colors','kkkk', 'Symbol','k+')
    xlim([.8 1.2])
    xticks([1])
    xticklabels({'Interpore'})
    ylim([ -0.5 clades_dmrcas_max(this_clade) ])
    ylabel('Pairwise dMRCA')
    set(gca,'FontSize',20)
    % Title
    suptitle([ 'Clade ' num2str(this_clade_num) ...
        ': Npores=' num2str(clades_numpores(this_clade)) ] )
    
    % Save figure
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [8 6]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 8 6]);
    print([ '4_pore_biogeo/' 'ExtraPlot_Clade-' num2str(this_clade_num) '_intra-v-inter.png'],'-dpng')

end


%% Wilcoxon rank-sum test

min_pores = 5;
clades_keep = [];

% Wilcoxon rank-sum test with Bonferroni correction
bonferroni_factor = sum(clades_numpores>=min_pores);
wilcoxon_pval_by_clade = [];
wilcoxon_pval_exponent_by_clade_with_bonferroni = [];

% dMRCA data
list_all_intrapore_raw = [];

% Naming schemes
% For clade assignments
list_all_intraclade = {};
list_all_interclade = {};
% For clade assignments with number of pores included
list_all_intracladenum = {};
list_all_intercladenum = {};
% For clade assignments with number of pores and inter/intra included
list_all_intracladenumtag = {};
list_all_intercladenumtag = {};
% For combined inter and intra
list_all_dmrcas = [];
list_all_dmrcas_raw = [];
list_all_cladenumtag = {};
% Just a list of clades
list_all_clades = {}; 
list_all_clades_wnum = {}; 

for this_clade = 1:numel(clades_cladenums)

    this_clade_num = clades_cladenums(this_clade);
    
    if clades_numpores(this_clade) < min_pores
        continue
    end
    
    % Clade short name
    if this_clade_num > 0
        this_clade_name = [ clusters_all_subjects{this_clade_num} '-' num2str(this_clade) ];
    else
        this_clade_name = 'Other';
    end
    fprintf(1,['Clade: ' this_clade_name '\n'])
    clades_keep(end+1) = this_clade;

    % Get dMRCA info
    next_clade_intra = clades_intrapore{this_clade};
    next_clade_inter_pairs = clades_interpore_pairs{this_clade};
    next_clade_inter_min = clades_interpore_min{this_clade};
    
    % Do Wilcoxon rank-sum test
    [wilcoxon_pval,~,~] = ranksum(next_clade_intra,next_clade_inter_pairs,'tail','left');
    wilcoxon_pval_by_clade(end+1) = wilcoxon_pval;
    wilcoxon_pval_exponent_by_clade_with_bonferroni(end+1) = ceil(log10(bonferroni_factor*wilcoxon_pval));
    
    % Append dMRCA info; NORMALIZE TO CLADE MAX
%     list_all_intrapore = [ list_all_intrapore; next_clade_intra/clades_dmrcas_max(this_clade) ];
%     list_all_interpore_pairs = [ list_all_interpore_pairs; next_clade_inter_pairs/clades_dmrcas_max(this_clade) ];
%     list_all_interpore_min = [ list_all_interpore_min; next_clade_inter_min/clades_dmrcas_max(this_clade) ];
%     list_all_dmrcas = [ list_all_dmrcas; next_clade_intra/clades_dmrcas_max(this_clade) ];
%     list_all_dmrcas = [ list_all_dmrcas; next_clade_inter_pairs/clades_dmrcas_max(this_clade) ];
    
    % Append dMRCA info; NOT NORMALIZED TO CLADE MAX
%     list_all_intrapore_raw = [ list_all_intrapore_raw; next_clade_intra ];
%     list_all_interpore_raw_pairs = [ list_all_interpore_raw_pairs; next_clade_inter_pairs ];
%     list_all_interpore_raw_min = [ list_all_interpore_raw_min; next_clade_inter_min ];
    list_all_dmrcas_raw = [ list_all_dmrcas_raw; next_clade_intra ];
    list_all_dmrcas_raw = [ list_all_dmrcas_raw; next_clade_inter_pairs ];

    % which clade
    list_all_clades{end+1} = this_clade_name;
    list_all_clades_wnum{end+1} = [ this_clade_name ' (n=' num2str(clades_numpores(this_clade)) ')']; 
    for i=1:numel(next_clade_intra)
        list_all_intraclade{end+1} = this_clade_name;
    end
    for i=1:numel(next_clade_inter_pairs)
        list_all_interclade{end+1} = this_clade_name;
    end
    % also with number of pores
    for i=1:numel(next_clade_intra)
        list_all_intracladenum{end+1} = [this_clade_name ' (n=' num2str(numel(next_clade_intra)) ')'];
    end
    for i=1:numel(next_clade_inter_pairs)
        list_all_intercladenum{end+1} = [this_clade_name ' (n=' num2str(clades_numpores(this_clade)) ')'];
    end
    % also with number of pores and inter/intra
    for i=1:numel(next_clade_intra)
        list_all_intracladenumtag{end+1} = [this_clade_name '-Intra (n=' num2str(numel(next_clade_intra)) ')'];
        list_all_cladenumtag{end+1} = [this_clade_name '-Intra (n=' num2str(numel(next_clade_intra)) ')'];
    end
    for i=1:numel(next_clade_inter_pairs)
        list_all_intercladenumtag{end+1} = [this_clade_name '-Inter (n=' num2str(clades_numpores(this_clade)) ')'];
        list_all_cladenumtag{end+1} = [this_clade_name '-Inter (n=' num2str(clades_numpores(this_clade)) ')'];
    end

end


%% PLOT HISTOGRAM OF INTRAPORE DMRCAS BY LINEAGE

% Binning
my_bin_edges = -0.5:1:22.5;
my_bin_centers = 0:1:22;
my_counts = [];
for i=1:numel(clades_intrapore)
    my_counts(end+1,:) = histcounts( clades_intrapore{i}, my_bin_edges);
end

% Appearance
fs = 18;
my_colors = [ ...
    239,243,255; ...
    189,215,231; ...
    107,174,214; ...
    49,130,189; ...
    8,81,156 ]/256; % https://colorbrewer2.org/#type=sequential&scheme=Blues&n=5

% Figure
fig=figure(3);
clf(3)
hold on
box on
b=bar(my_counts', 'stacked', 'FaceColor', 'flat' );
for i=1:numel(clades_intrapore)
    b(i).CData = my_colors(5-i+1,:);
end
xticks(1:1:numel(my_bin_centers))
xticklabels(my_bin_centers)
xlabel('intrapore dMRCA')
ylabel('number of pores')
lgd=legend( {'1a','1b','3a','3b','other'} );
lgdtit = get(lgd,'Title');
set(lgdtit,'String','lineage')
set(gca,'LineWidth',1)
set(gca,'FontName','Helvetica','FontSize',fs)
hold off

% Save
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [10 5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 5]);
print([ '4_pore_biogeo/' 'Hist_Intrapore_split.png'],'-r400','-dpng')


%% PLOT COMPARING INTRA- AND INTER- PORE DMRCAS 

legend_bool = false;

fig=figure(5);
clf(5)
box on
fs = 20;
% Intrapore and Interpore together
b=boxplot(...
    list_all_dmrcas_raw, list_all_cladenumtag, ... % 'Orientation', 'horizontal',
    'Colors','kkkk', 'Widths',0.65, 'Symbol','ko', 'Jitter',0.5 );
set(b,'LineWidth',1.5);
set(gca,'FontSize',fs)
% % X axis
xticks(1.5:2:2*numel(list_all_clades)-0.5)
list_all_clades_wnum_manual = { '1a', '1b', '3a', '3b', 'other' }; % clades_numpores
%list_all_clades_wnum_manual = { 'A-1', 'A-2', 'F-1', 'F-2', 'other' }; % clades_numpores
xticklabels(list_all_clades_wnum_manual)
xlabel('lineage','FontSize',fs)
% % Y axis
ylim([-1 38])
ax = gca;
ax.YAxis.FontSize = fs-4;
ylabel('dMRCA (#SNVs)','FontSize',fs)
yticks(0:8:40)
% % Boxes
color_plots = [.1 .4 .5];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    if rem(j,2)==0 % evens are inter
        patch(get(h(j),'XData'),get(h(j),'YData'),color_plots,'FaceAlpha',.85,'LineWidth',1.5); % all gray
    else % odds are intra 
        patch(get(h(j),'XData'),get(h(j),'YData'),color_plots,'FaceAlpha',.15,'LineWidth',1.5); % all gray
    end
end
% % Legend % for some reason has to be here to get line to show up
if legend_bool
    lgd=legend('interpore pairwise dMRCA','intrapore mean dMRCA','Location','northoutside');
    lgd.FontSize = fs-4;
    %lgd.Box = 'off';
end
% % Add KS test results
wilcoxon_height = 36;
wilcoxon_fontsize = fs-6;
for ci=1:numel(list_all_clades)
    if wilcoxon_pval_exponent_by_clade_with_bonferroni(ci)<=-2
        text(2*ci-0.5,wilcoxon_height,[' p < 10^{' num2str(max(-10,wilcoxon_pval_exponent_by_clade_with_bonferroni(ci))) '}'],'FontSize',wilcoxon_fontsize,'HorizontalAlignment','center')
    else
        text(2*ci-0.5,wilcoxon_height-0.33,['ns'],'FontSize',wilcoxon_fontsize,'HorizontalAlignment','center')
    end
end
wilcoxon_line_height = wilcoxon_height-2;
lw_pval = 1.5;
for ci=1:numel(list_all_clades)
    line( [2*ci-1, 2*ci], [wilcoxon_line_height, wilcoxon_line_height], ...
        'HandleVisibility', 'off', 'Color', 'k', 'LineWidth', lw_pval )
    line( [2*ci-1, 2*ci-1], [wilcoxon_line_height-1, wilcoxon_line_height], ...
        'HandleVisibility', 'off', 'Color', 'k', 'LineWidth', lw_pval)
    line( [2*ci, 2*ci], [wilcoxon_line_height-1, wilcoxon_line_height], ...
        'HandleVisibility', 'off', 'Color', 'k', 'LineWidth', lw_pval)
end
set(gca,'LineWidth',1)
set(gca,'FontName','Helvetica')

% Save figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 6 5]);
if legend_bool
    print([ '4_pore_biogeo/' 'BoxPlot_IntraPore-v-InterPore_lgd.png'],'-r400','-dpng')
else
    print([ '4_pore_biogeo/' 'BoxPlot_IntraPore-v-InterPore.png'],'-r400','-dpng')
end


%% PLOT DISTANCE TO CLOSEST RELATED PORE

legend_bool = false;

fig=figure(8);
clf(8)
box on
fs = 20;
% Histogram for each set
%titles_manual = { 'A-1', 'A-2', {'F-1','lineage'}, 'F-2', 'other' };
titles_manual = { '1a', '1b', '3a', '3b', 'other' };
for s=1:numel(titles_manual)
    subplot(1,numel(titles_manual),s)
    box on
    bin_width = 2;
    bins = 0:bin_width:40; 
    counts = histcounts(2*clades_interpore_min{s},bins); % factor of 2 to convert joint dMRCA to distance
    b=barh( counts, 'FaceColor', color_plots, 'FaceAlpha', 0.5, 'LineWidth', 1 );
    set(b,'LineWidth',1.5);
    % y axis
    ylim([0.5 numel(bins)-0.5])
    yticks(0.5:1:numel(bins)+1)
    if s==1
        bins_labels = arrayfun(@(x) num2str(x), bins, 'UniformOutput', false);
        for i=1:numel(bins)
            if rem(i-1,4) ~=0
                bins_labels{i} = '';
            end
        end
        yticklabels(bins_labels)
    else
        yticklabels('')
    end
    % x axis
    xlim([ 0 1.1*max(counts) ])
    xticks(0:1:max(counts))
    xticks_labels = arrayfun(@(x) num2str(x), 0:1:max(counts), 'UniformOutput', false);
    for i=1:numel(xticks_labels)
        if i>1 && i<max(counts)+1
            xticks_labels{i} = '';
        end
    end
    xticklabels(xticks_labels)
    if s<numel(titles_manual)
        line( [0 max(xlim)], [2*median(clades_interpore_pairs{s})/bin_width+0.5, 2*median(clades_interpore_pairs{s})/bin_width+0.5], 'Color', color_plots, 'LineStyle', '-', 'LineWidth', 2 )
    end
    % titles
    t=title(titles_manual{s},'FontWeight','normal','FontSize',fs);
    % labeling
    set(gca,'FontSize',fs)
    ax = gca;
    ax.TitleFontSizeMultiplier=1;
    ax.YAxis.FontSize = fs-4;
    ax.XAxis.FontSize = fs-4;
    if s==1%numel(titles_manual)
        ylabel({'#SNVs to most closely related pore'}, 'FontSize', fs)
    end
    if s==3
        xlabel({'number of pores'},'FontSize',fs)
    else
        xlabel('','FontSize',fs)
    end
    set(gca,'LineWidth',1)
    set(gca,'FontName','Helvetica')
    % legend
    if legend_bool
        if s~=4
            lgd = legend({'',''},'Location','northoutside');
            lgd.FontSize = fs-4;
        else
            lgd = legend({'interpore minimum distance','median interpore distance (#SNVs)'},'Location','northoutside');
            lgd.FontSize = fs-4;
        end
    end
end

% Save figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 6 5]);
if legend_bool
    print([ '4_pore_biogeo/' 'Hists_InterPore_Min-sideways_lgd.png'],'-r400','-dpng')
else
    print([ '4_pore_biogeo/' 'Hists_InterPore_Min-sideways.png'],'-r400','-dpng')
end



%% LINEAGE MEDIAN dMRCA HISTOGRAM

figure(10)
clf(10)
fs = 26;
lw = 1.5;
hold on
box on
histogram( clades_dmrcas_median_all, 0:2.5:40,'FaceColor', [.5,.5,.5], 'EdgeColor', 'k', 'LineWidth', lw  )
% X
xlabel('median lineage dMRCA (SNVs)')
xticks( 0:5:40 )
% Y
ylabel('number of lineages')
% Formatting
set(gca,'FontSize',fs)
hold off

% Save figure
print([ '4_pore_biogeo/' 'Hist_LineageMRCA.png'],'-dpng')


%%

min(clades_dmrcas_median_all)
max(clades_dmrcas_median_all)

% ans =
% 
%      0
% 
% 
% ans =
% 
%     26


