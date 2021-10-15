%% SUMMARY

% Looking at data tables from Oh et al 2014 in the context of strain-type
% (super SLST)


%% SETUP

% Add matlab functions to path
path(path,'miniscripts')

% Folder for storing figures
mkdir('figures')


%% IMPORT DATA

% All data from Oh Segre biogeography paper

% Import metadata
xls_meta = readtable('data/ohsegre_biogeo_s1.xls'); % issue with xlsx; converted

% Import abundances
xls_taxa = readtable('data/ohsegre_biogeo_s11.xls'); % issue with xlsx; converted
% Because excel/matlab issue prevents more than 256 columns from being
% imported 
xls_taxa_2 = readtable('data/ohsegre_biogeo_s11_256.xls'); % issue with xlsx; converted
% FMI: https://www.mathworks.com/matlabcentral/answers/132062-how-to-import-csv-file-of-200-rows-by-4096-columns-without-truncation

% Import strain key (via Figure 1 in Scholz et al 2014)
xls_strainkey = readtable('data/strain_key.xlsx'); % for conversion to SLST


%% PROCESS METADATA

% Grab metadata
sample_names = xls_meta.SAMPLE;
sample_subjects = xls_meta.PatientID; % ignore 'MOCK' and 'SH' (not adults)
sample_sites = xls_meta.Site_Symmetry;
sample_characteristic = xls_meta.SiteCharacteristic;

% Filter metadata 
% Remove mock communities and children and 'Ax' site
samples_to_keep  = ~cell2mat( cellfun(@(x) contains(x,'MOCK'), sample_subjects, 'UniformOutput', false) ) ...
    & ~cell2mat( cellfun(@(x) contains(x,'Mock'), sample_subjects, 'UniformOutput', false) ) ...
    & ~cell2mat( cellfun(@(x) contains(x,'SH'), sample_subjects, 'UniformOutput', false) ) ...
    & ~cell2mat( cellfun(@(x) contains(x,'Ax'), sample_sites, 'UniformOutput', false) );
% Update variables
sample_names = sample_names( samples_to_keep );
sample_subjects = sample_subjects( samples_to_keep );
sample_sites = sample_sites( samples_to_keep );
sample_characteristic = sample_characteristic( samples_to_keep );

% Clean up metadata
% Remove symmetry information from sample sites
sample_sites = cellfun(@(x) strtok(x,'-'), sample_sites, 'UniformOutput', false); % only keep info before dash
% Get information about unique sample sites
key_sample_sites_names = unique( sample_sites );
key_sample_sites_characteristics = cellfun(@(x) unique(sample_characteristic(ismember( sample_sites, x))), key_sample_sites_names);


%% PROCESS STRAIN DATA

% Get information from strain key
key_names = xls_strainkey.StrainName;
key_slst = xls_strainkey.SLST;
key_tradtype = xls_strainkey.Traditional;

% Get strain names from taxa table
xls_taxa_names_0 = xls_taxa.Taxa;
% Only keep substring after '_'
xls_taxa_names = cellfun(@(x) reverse(strtok(reverse(x),'_')), xls_taxa_names_0, 'UniformOutput', false);
% Fix names that got messed up
xls_taxa_names{7} = 'HL037PA1'; % originally annoted in xls as 'HL037PA2' but tree and SLST paper have 'HL037PA1'...?
xls_taxa_names{38} = 'ATCC-11828'; % fix name parsing
xls_taxa_names{63} = 'PRP-38'; % fix name parsing
xls_taxa_names{76} = 'FZ1-2-0'; % fix name parsing 
% Check if these all match one name in the key
% thing=cellfun(@(x) ismember(x,key_names), xls_taxa_names )
% xls_taxa_names(~thing)
% Remove header lines

% Extract taxa table
taxa_table_1 = table2array( xls_taxa(4:end,2:end)); 
taxa_table_names = xls_taxa_names(4:end); % rows
taxa_table_samples_1 = xls_taxa.Properties.VariableNames(2:end); % columns
% And the rest of the taxa table
taxa_table_2 = table2array( xls_taxa_2(4:end,2:end)); 
taxa_table_samples_2 = xls_taxa_2.Properties.VariableNames(2:end); % columns
% Combine
taxa_table = [ taxa_table_1, taxa_table_2 ];
taxa_table_samples = horzcat( taxa_table_samples_1, taxa_table_samples_2 );

% Extract number of reads
sample_reads = [ table2array( xls_taxa(1,2:end) ), table2array( xls_taxa_2(1,2:end) ) ];

%%

% Remove taxon that doesn't match anything in the key and also has zero
% abundance
% 'HL044PA1'
%taxa_table_names{21-3};
%sum(taxa_table(21-3,:))
taxa_to_remove = zeros( size(taxa_table_names) );
taxa_to_remove(21-3) = 1;
% Downsize
taxa_table = taxa_table( ~taxa_to_remove,: );
taxa_table_names = taxa_table_names( ~taxa_to_remove );

% Reorder according to SLST tree
% Get order
key_indices_1 = cell2mat(cellfun(@(x) find(ismember(key_names,x)), taxa_table_names, 'UniformOutput', false ));
[ ~,key_indices_2 ] = sort(key_indices_1);
% Do reordering
taxa_table = taxa_table( key_indices_2,: );
taxa_table_names = taxa_table_names( key_indices_2 );
% Also get order typing for taxa_table_names
taxa_table_names_slst = key_slst(key_indices_1);
taxa_table_names_slst = taxa_table_names_slst(key_indices_2);
taxa_table_names_tradtype = key_tradtype(key_indices_1);
taxa_table_names_tradtype = taxa_table_names_tradtype(key_indices_2);



%% CONVERT ABUNDANCES TO SLSTs

% Unique SLSTs
list_slsts = unique( taxa_table_names_slst );
taxa_table_by_slst = zeros( numel(list_slsts), numel(taxa_table_samples) );
for i=1:numel(list_slsts)
    next_slst = list_slsts{i};
    next_slst_indices = ismember( taxa_table_names_slst, next_slst );
    taxa_table_names_slst(next_slst_indices);
    if sum(next_slst_indices) == 1
        next_row = taxa_table( next_slst_indices,: );
    else
        next_row = sum( taxa_table( next_slst_indices,: ) );
    end
    taxa_table_by_slst(i,:) = next_row;
end

% Unique super SLSTs
taxa_table_names_slst_super = cellfun(@(x) x(1), taxa_table_names_slst );
list_slsts_super = unique( taxa_table_names_slst_super );
taxa_table_by_slst_super = zeros( numel(list_slsts_super), numel(taxa_table_samples) );

for i=1:numel(list_slsts_super)
    next_slst_super = list_slsts_super(i);
    next_slst_indices = ismember( taxa_table_names_slst_super, next_slst_super );
    taxa_table_names_slst_super(next_slst_indices);
    if sum(next_slst_indices) == 1
        next_row = taxa_table( next_slst_indices,: );
    else
        next_row = sum( taxa_table( next_slst_indices,: ) );
    end
    taxa_table_by_slst_super(i,:) = next_row;
end



%% FIGURES


%% Define super SLST colors

colors_super = [ ... % seven letters
    [ 215,48,39 ]; ... % A
    [ 244,109,67 ]; ... % B *new
    [ 253,174,97 ]; %[ 252,141,89]; ... % C
    [ 254,224,144]; ... % D
    [ 255,255,191 ]; ... % E
    [ 224,243,248 ]; ... % F
    [ 171,217,233 ]; ... % G *new
    [ 116,173,209]; %[ 145,191,219 ]; ... % H
    [ 69,117,180 ]; ... % K
    [ 200,200,200 ]; % unknown
    ];
colors_super = colors_super/256;


%% Make a big plot for all subjects, separating skin sites by characteristic, manual order
% Only sebaceous sites

% Get lists of all subjects and all sample sites
list_subjects = unique( sample_subjects );
[ list_sites, temp_indices ] = unique( sample_sites );
list_characteristics = sample_characteristic( temp_indices );
% horzcat( list_sites, list_characteristics ) % check site characteristic matching

% Reorder in terms of characteristic (dry, moist, sebaceous, toenail)
% [ list_characteristics, temp_indices ] = sort( list_characteristics );
% list_sites = list_sites( temp_indices );

% Reorder manually
list_sites = { 'Gb', 'Ch', 'Al', 'Ba' }; % oily only
%list_sites_long = { 'Gb (face)', 'Ch (face)', 'Al (face)', 'Ba (back)' }; % oily
list_sites_long = { 'Gb (face)', 'Ch (face)', 'Al (face)', 'Ba (back)' }; % oily

% Initialize figure
fig3=figure(3);
clf(3)

for s=1:numel(list_subjects)
    
    this_subject = list_subjects{s}; 

    % Get this subject's sample info from metadata
    this_subject_samples = sample_names( cellfun(@(x) isequal(x,this_subject), sample_subjects) );
    this_subject_samples_sites = sample_sites( cellfun(@(x) isequal(x,this_subject), sample_subjects) );
    this_subject_samples_characteristics = sample_characteristic( cellfun(@(x) isequal(x,this_subject), sample_subjects) );

    % Find these samples in the taxa table
    %taxa_table_indices = cellfun(@(x) find(ismember(taxa_table_samples,x)), this_subject_samples);
    taxa_table_indices = cellfun(@(x) find(ismember(taxa_table_samples,x)), this_subject_samples, 'UniformOutput', false);
    taxa_table_indices_keep = cellfun(@(x) numel(x),taxa_table_indices)>0;
    taxa_table_indices = cell2mat(taxa_table_indices(taxa_table_indices_keep)); % in case can't find some of the samples
    % Update other variables if missing some samples
    this_subject_samples_sites = this_subject_samples_sites(taxa_table_indices_keep);
    this_subject_samples_characteristics = this_subject_samples_characteristics(taxa_table_indices_keep);
    
    % Extract taxa table for this subject
    this_subject_taxa_table = taxa_table_by_slst_super( :,taxa_table_indices );

    % Add gray bar for leftover relative abundance
    this_subject_taxa_table = [ this_subject_taxa_table; max( 0, 100 - sum( this_subject_taxa_table ) ) ];
    
    % Fill in gaps for missing sites
    this_subject_taxa_table_full = zeros( numel(list_slsts_super)+1,numel(list_sites) );
    for i=1:numel(list_sites)
        next_site = list_sites{i};
        next_site_index = find(ismember( this_subject_samples_sites, next_site ));
        if numel(next_site_index) == 1
            this_subject_taxa_table_full(:,i) = this_subject_taxa_table( :,next_site_index );
        elseif numel(next_site_index) > 1 % if more than one sample for same site
            next_site_index = next_site_index(1); % just choose the first
            this_subject_taxa_table_full(:,i) = this_subject_taxa_table( :,next_site_index );
        else
            this_subject_taxa_table_full(:,i) = zeros( numel(list_slsts_super)+1,1 );
        end
    end
    
    % Add rows of zeros between site types
    this_subject_taxa_table_full_gaps = ...
        [ this_subject_taxa_table_full(:,1:3), ...
        zeros(numel(list_slsts_super)+1,1) ...
        this_subject_taxa_table_full(:,4) ];

    % Make a bar chart
    %[~,reorder] = sort(list_sites); % alphabetical order
    reorder = 1:1:numel(list_sites)+1; % use what was defined above
    subplot(3,5,s)
    box on
    hold on
    b=bar(this_subject_taxa_table_full_gaps(:,reorder)','stacked','LineWidth',1.5);
    for i=1:numel(list_slsts_super)+1
        b(i).FaceColor = colors_super(i,:);
    end
    % lines
    line_width = 1;
    set(gca,'LineWidth',line_width)
    % title
    title(['Subject: ' this_subject])
    % x axis
    xticks([1 2 3 5]); %xticks(1:1:numel(list_sites)+1)
    xticklabels(list_sites)
    xtickangle(0)
    xlabel('       face               back')
%    xlabel('Sample Location')
    % y axis
    ylim([0 100])
    if ismember( s, [1,6,11] )
        ylabel('relative abundance')
    end
    % font
    set(gca,'FontSize',14)
    % legend
%    legend(b,list_slsts,'Location','eastoutside')
    hold off

end

suptitle('Cutibacterium acnes abundance data from Oh et al 2014 (sebaceous only)')

% Save
print(fig3,[ 'figures' '/' 'BarCharts_SiteSebaceous.png' ],'-dpng')


%% BRAY-CURTIS for FACE/BACK

%% Compute Bray-Curtis between face vs face-back

% Get read counts by super SLST for each sample
taxa_table_by_slst_super_reads = round( taxa_table_by_slst_super.*sample_reads );
num_samples = size( taxa_table_by_slst_super_reads,2 );

% Compute Bray-Curtis
bray_curtis_mat = zeros( num_samples, num_samples );
for i=1:num_samples
    for j=1:num_samples
        bray_curtis_mat(i,j) = compute_bray_curtis( taxa_table_by_slst_super_reads(:,i), taxa_table_by_slst_super_reads(:,j) );
    end
end

% Booleans: face/back
list_sites_face = { 'Gb', 'Ch', 'Al' }; 
face_samples = sample_names( ismember( sample_sites, list_sites_face ) );
is_face = ismember( taxa_table_samples, face_samples );
list_sites_back = { 'Ba' }; 
back_samples = sample_names( ismember( sample_sites, list_sites_back ) );
is_back = ismember( taxa_table_samples, back_samples );
% Booleans: same/diff subject
is_same_subject = zeros( num_samples, num_samples, 'logical' );
for i=1:num_samples
    for j=1:num_samples
        samp_i = taxa_table_samples{i};
        samp_j = taxa_table_samples{j};
        if ismember( samp_i, sample_names ) && ismember( samp_j, sample_names )
        subj_i = sample_subjects( ismember(sample_names,samp_i) );
        subj_j = sample_subjects( ismember(sample_names,samp_j) );
        is_same_subject(i,j) = isequal( subj_i, subj_j );
        else
            is_same_subject(i,j) = false;
        end
    end
end

%%

% Get Bray-Curtis for same subject, different face sites
bc_same_subj_face_sites = bray_curtis_mat( ...
    is_same_subject ...
    & ( is_face & is_face' ) ...
    & 1:1:num_samples < [1:1:num_samples]' ...
    );
% Get Bray-Curtis for same subject, face+back sites
bc_same_subj_face_back = bray_curtis_mat( ...
    is_same_subject ...
    & ( is_face & is_back' ) ...
    & 1:1:num_samples < [1:1:num_samples]' ...
    );

% Get Bray-Curtis for different subject, face sites
bc_diff_subj_face_face = bray_curtis_mat( ...
    ~is_same_subject ...
    & ( is_face & is_face' ) ...
    & 1:1:num_samples < [1:1:num_samples]' ...
    );
% Get Bray-Curtis for different subject, back sites
bc_diff_subj_back_back = bray_curtis_mat( ...
    ~is_same_subject ...
    & ( is_back & is_back' ) ...
    & 1:1:num_samples < [1:1:num_samples]' ...
    );
% Get Bray-Curtis for different subject, face+back sites
bc_diff_subj_face_back = bray_curtis_mat( ...
    ~is_same_subject ...
    & ( is_face & is_back' ) ...
    & 1:1:num_samples < [1:1:num_samples]' ...
    );

ranksum( bc_same_subj_face_sites, bc_same_subj_face_back ) % 0.0097
ranksum( bc_diff_subj_face_face, bc_diff_subj_back_back ) % 0.2435
ranksum( bc_diff_subj_face_back, bc_diff_subj_face_face ) % 0.8549
ranksum( bc_diff_subj_face_back, bc_diff_subj_back_back ) % 0.3341

%%

% Make plot
figure(1); clf(1);
lw = 1.5;
fs = 20;
bar_color = 0.75*[ 1 1 1 ];
%
subplot(2,1,1)
box on
hold on
bins= 0:0.1:1;
histogram( bc_same_subj_face_sites, bins, 'LineWidth', lw, 'FaceColor', bar_color )
ylim([0 12])
%xlabel('Bray-Curtis dissimilarity')
ylabel('# pairs of sites')
title('face sites on the same subject')
set(gca,'LineWidth',line_width,'FontSize',fs)
hold off
%
subplot(2,1,2)
box on
hold on
bins= 0:0.1:1;
histogram( bc_same_subj_face_back, bins, 'LineWidth', lw, 'FaceColor', bar_color  )
ylim([0 7])
xlabel('Bray-Curtis dissimilarity')
ylabel('# pairs of sites')
title('face vs back sites on the same subject')
set(gca,'LineWidth',line_width,'FontSize',fs)
hold off

print([ 'figures' '/' 'BrayCurtis_same-subj.png' ],'-dpng')


% Make plot
figure(2); clf(2);
lw = 1.5;
fs = 20;
bar_color = 0.75*[ 1 1 1 ];
%
subplot(3,1,1)
box on
hold on
bins= 0:0.1:1;
histogram( bc_diff_subj_face_face, bins, 'LineWidth', lw, 'FaceColor', bar_color )
ylim([0 400])
%xlabel('Bray-Curtis dissimilarity')
%ylabel('# pairs of sites')
title('face sites on different subject')
set(gca,'LineWidth',line_width,'FontSize',fs)
hold off
%
subplot(3,1,2)
box on
hold on
bins= 0:0.1:1;
histogram( bc_diff_subj_back_back, bins, 'LineWidth', lw, 'FaceColor', bar_color )
ylim([0 25])
%xlabel('Bray-Curtis dissimilarity')
ylabel('# pairs of sites')
title('back sites on different subject')
set(gca,'LineWidth',line_width,'FontSize',fs)
hold off
%
subplot(3,1,3)
box on
hold on
bins= 0:0.1:1;
histogram( bc_diff_subj_face_back, bins, 'LineWidth', lw, 'FaceColor', bar_color  )
ylim([0 110])
xlabel('Bray-Curtis dissimilarity')
%ylabel('# pairs of sites')
title('face vs back sites on different subject')
set(gca,'LineWidth',line_width,'FontSize',fs)
hold off

print([ 'figures' '/' 'BrayCurtis_diff-subj.png' ],'-dpng')
