%% Script to make relative abundance from 16S amplicon sequencing data


%% Setup

path( path, 'matlab_functions' )


%% Parameters 

% Sample filtering
min_num_reads = 100; 

% Taxa to show individually
min_abundance_to_show_species = 0.05; % average over all samples
min_abundance_to_show_genus = 0.33; % max over all samples


%% Processing 16S data

% Import data
rawdata = readtable('data/final_OTU_table.xlsx');

% Sample names
column_names = rawdata.Properties.VariableNames;
sample_names = string(column_names(5:end));

% Taxon names and taxon counts
taxons = string(rawdata{:,4});
taxon_counts = table2array(rawdata(:,5:end));

% Filter and condense taxons
% 1. Combine repeate taxon names
unique_taxons = unique(taxons,'stable');
condensed_counts = zeros(length(unique_taxons), length(sample_names));
for ii = 1:length(unique_taxons)
    condensed_counts(ii,:) = sum(taxon_counts(taxons==unique_taxons(ii),:),1);
end

% 2. Remove unidenitified bactera and spike ins (Caulobacter) 
taxa_to_remove = ( unique_taxons == "Bacteria" ) | ( unique_taxons == "Unassigned" ) | ( contains(unique_taxons,"Caulobacter") );
condensed_counts(taxa_to_remove,:) = []; 
unique_taxons(taxa_to_remove) = []; 

% Filter samples
valid_patient = contains(sample_names, ["A_","F_"]) ... % ,"G_"
    & ( ~contains(sample_names, {'posctrl','blank'}) & ~contains( sample_names, {'F_Sc_16_H11','F_St_G08'} ) ); % more blanks
enough_coverage = sum(condensed_counts)>min_num_reads; % has at least nnn reads in a sample
sample_names = sample_names(valid_patient & enough_coverage);
condensed_counts = condensed_counts(:,valid_patient & enough_coverage);

% Relative abundance in each sample
condensed_counts_normalized = condensed_counts./sum(condensed_counts); % normalize by sample


%% Get info for sample_names_new

% Load data
namedata = readtable('data/Specimen_log_simple.csv');
table_specnum = namedata.SpecimenNumber;
table_spec = namedata.Specimen;
table_subject = namedata.Subject;
table_time = namedata.SamplingTime;

% Update names
sample_names_old_short = cell( size(sample_names) ); % initialize
sample_types = cell( size(sample_names) ); % initialize
sample_names_new_long = cell( size(sample_names) ); % initialize
G_tally = 0; % since G doesn't have specimen numbers in the log
for n=1:numel(sample_names)
    old_name = sample_names{n};
    this_subject = old_name(1);
    if contains( old_name, 'Ex' ) || contains( old_name, 'Sc' )
        this_tag = old_name(end-5:end-4);
        if this_tag(1)=='0'
            this_tag = this_tag(2);
        end
    elseif contains( old_name, 'St' )
        this_tag = old_name(end-2:end);
    end
    this_bool = ismember( table_spec, this_tag ) & ismember( table_subject, this_subject ) & cellfun(@(x) contains(x,'2016'), table_time );
    if sum(this_bool)~=1 % unique match not found
        fprintf(1, [ old_name ':' num2str(sum(this_bool)) '\n'] )
        next_name_old_short = old_name;
        next_name_new_long = old_name;
        next_name_type = old_name;
    else
        next_name_old_short = [ old_name(1:5) num2str(table_specnum(this_bool)) ];
        next_name_new_long = [ old_name(1:5) num2str(table_specnum(this_bool)) ];
        next_name_type = old_name;
    end
    % Get sample type
    if contains( next_name_type, 'Sc_' )
        next_name_type = 'skin_scrape';
    elseif contains( next_name_type, 'St_' )
        next_name_type = 'pore_strip';
    elseif contains( next_name_type, 'Ex_' )
        next_name_type = 'pore_extract';
    else
        next_name_type = 'unrecognized_type';
    end
    % Update to subject numbers
    next_name_new_long = strrep( next_name_new_long, 'A_', 'subj-1_' );
    next_name_new_long = strrep( next_name_new_long, 'F_', 'subj-3_' );
    next_name_new_long = strrep( next_name_new_long, 'G_', 'subj-17_' );
    % Update sample descriptions
    next_name_new_long = strrep( next_name_new_long, '_Sc_', '_scrape-' );
    next_name_new_long = strrep( next_name_new_long, '_St_', '_pore-' );
    next_name_new_long = strrep( next_name_new_long, '_Ex_', '_pore-' );
    % Add 16S to long nname
    next_name_new_long = [ next_name_new_long '_16S' ];
    % Save new name
    sample_types{n} = next_name_type;
    sample_names_old_short{n} = next_name_old_short;
    sample_names_new_long{n} = next_name_new_long;
end

% Match back to original filenames (needed for SRA upload)
filedata = readtable('data/MANIFEST.csv','Delimiter',',');
filedata_filenames = filedata.filename;
file_index = [];
for i=1:numel(sample_names)
    next_name_tag = strrep( sample_names{i},'_','-' );
    next_index = find( cellfun(@(x) contains(x,next_name_tag), filedata_filenames ) );
    file_index = [ file_index next_index ];
end

% Write csv with key to different sample nomenclature
fid = fopen( 'sample_names_key.csv', 'w' );
fprintf(fid, [ 'sample_names_raw,sample_names_old_short,sample_names_new_long,sample_type,filename,' '\n' ] );
for i=1:numel(sample_names)
    fprintf(fid, [ sample_names{i} ',' ...
        sample_names_old_short{i} ',' ...
        sample_names_new_long{i} ',' ...
        sample_types{i} ',' ...
        filedata_filenames{file_index(i)} ',' '\n' ] );
end
fclose(fid);


%% Get classifications at each taxonomic levels

% Parse taxa labels to get classification at each taxonomic level available
split_by_classification = cellfun(@(x) split(x, ';'), unique_taxons, 'UniformOutput', false);
taxonomic_classes_split = repmat({strings(size(split_by_classification))}, 7,1); % initialize
for ii=1:7 % check all taxonomic levels
    is_classified_in_class = cellfun(@length, split_by_classification) >=ii; % boolean for which samples have this level of classification
    taxonomic_classes_split{ii}(is_classified_in_class) = cellfun(@(x) x(ii), split_by_classification(is_classified_in_class));
end
taxonomic_classes_split{7} = replace(taxonomic_classes_split{7}, '__', ' ');

% Get species with mean relative abundance over min_abundance_to_show_species
list_species = taxonomic_classes_split{7};
unique_species = unique(taxonomic_classes_split{7}, "stable"); unique_species(unique_species == "") = [];
unique_species_mean_abundance = zeros( numel(unique_species), 1 );
for i=1:numel(unique_species)
    bool_taxon = ismember( list_species, unique_species{i} );
    if sum(bool_taxon)==1
        unique_species_mean_abundance(i) = mean( condensed_counts_normalized( bool_taxon, : ) );
    else
        unique_species_mean_abundance(i) = mean( sum( condensed_counts_normalized( bool_taxon, : ) ) );
    end
end
species_to_graph = unique_species( unique_species_mean_abundance >= min_abundance_to_show_species );

% Get genuses with mean relative abundance over min_abundance_to_show_genus
list_genuses = taxonomic_classes_split{6};
unique_genuses = unique(taxonomic_classes_split{6}, "stable"); unique_genuses(unique_genuses == "") = [];
unique_genus_mean_abundance = zeros( numel(unique_genuses), 1 );
unique_genus_max_abundance = zeros( numel(unique_genuses), 1 );
for i=1:numel(unique_genuses)
    bool_taxon = ismember( list_genuses, unique_genuses{i} ) & ~ismember( list_species, species_to_graph );
    if sum(bool_taxon)==1
        unique_genus_mean_abundance(i) = mean( condensed_counts_normalized( bool_taxon, : ) );
        unique_genus_max_abundance(i) = max( condensed_counts_normalized( bool_taxon, : ) );
    else
        unique_genus_mean_abundance(i) = mean( sum( condensed_counts_normalized( bool_taxon, : ) ) );
        unique_genus_max_abundance(i) = max( sum( condensed_counts_normalized( bool_taxon, : ) ) );
    end
end
%genuses_to_graph = unique_genuses( unique_genus_mean_abundance >= min_abundance_to_show_genus )
genuses_to_graph = unique_genuses( unique_genus_max_abundance >= min_abundance_to_show_genus );
%     "Lawsonella"
%     "Staphylococcus"
%     "Anaerococcus"
%     "Corynebacterium"
%     "Peptoniphilus"
%     "Delftia"


%% Create array for graphing with only genus level classifications

% NOTE: the ordering of the below keywords excludes any species implicated
% in a search before it - ex. "staphylococcus" is behind "aureus" and will
% therefore not count any "aureus" species.

% List of taxa to plot
graph_labels_taxa = vertcat( species_to_graph, genuses_to_graph );
number_catagories_to_show = length(graph_labels_taxa) + 1; % include a spot for "other"
species_seperated_counts_normalized = zeros(number_catagories_to_show, size(condensed_counts_normalized,2)); % initialize

% Fill in table for specified taxa
species_not_to_count = zeros(length(unique_taxons),1);
for ii=1:length(graph_labels_taxa)
    contains_spec_name = contains(unique_taxons, graph_labels_taxa(ii)) & ~species_not_to_count;
    species_not_to_count = species_not_to_count | contains_spec_name;
    species_seperated_counts_normalized(ii,:) = sum(condensed_counts_normalized(contains_spec_name,:),1);
end

% Include rest in "other" category
graph_labels = vertcat( graph_labels_taxa, {"Other"} );
species_seperated_counts_normalized(end,:) = sum(condensed_counts_normalized(~species_not_to_count,:),1);


%% Categorize samples for plotting

%  Sort by C acnes abundance 
[~,acnes_sorting] = sort(species_seperated_counts_normalized(1,:),'descend');

% By subject
subjectA = contains(sample_names, "A_");
subjectF = contains(sample_names, "F_");
%subjectG = contains(sample_names, "G_");

% By sample type
strip = contains(sample_names, "St");
scrape = contains(sample_names, "Sc");
extract = contains(sample_names, "Ex");
sampletype = {extract,strip,scrape};
sampletype_names = { 'pore extract', 'pore strip', 'scrape' };

% Make arrays containing all subject A/F samples sorted in order of C. acnes abundance
% The first cell contains strips, the next scrapes, and the third extracts
A_samples_to_show = {[];[];[]}; 
F_samples_to_show = {[];[];[]}; 
%G_samples_to_show = {[];[];[]}; 
for ii=1:numel(sampletype)
    A_samples_to_show{ii} = acnes_sorting(ismember(acnes_sorting,find(sampletype{ii} & subjectA)));
    F_samples_to_show{ii} = acnes_sorting(ismember(acnes_sorting,find(sampletype{ii} & subjectF)));
%    G_samples_to_show{ii} = acnes_sorting(ismember(acnes_sorting,find(sampletype{ii} & subjectG)));
end


%% Plotting

dir_plots = 'figs';
if ~exist( dir_plots, 'dir')
    mkdir( dir_plots )
end

legend_bool = false;

% Subject A / Subject 1
for i=1:numel(sampletype_names)
    
    figure(1)
    clf(1)
    these_samples = A_samples_to_show{i};
    these_samples_name = [ 'subject 1: ' sampletype_names{i} ];
    bar_temp = [species_seperated_counts_normalized(:,these_samples)];
    x_labels = sample_names_old_short(these_samples);
    
    if legend_bool
        fig_filename = [ dir_plots '/' '16S_subj-1_' sampletype_names{i} '-wleg.png'];
    else
        fig_filename = [ dir_plots '/' '16S_subj-1_' sampletype_names{i} '.png'];
    end
    make_miniplot( these_samples_name, bar_temp, x_labels, graph_labels, legend_bool, fig_filename )
    
end

% Subject F / Subject 3
for i=1:numel(sampletype_names)
    
    figure(1)
    clf(1)
    these_samples = F_samples_to_show{i};
    these_samples_name = [ 'subject 3: ' sampletype_names{i} ];
    bar_temp = [species_seperated_counts_normalized(:,these_samples)];
    x_labels = sample_names_old_short(these_samples);
    
    if legend_bool
        fig_filename = [ dir_plots '/' '16S_subj-3_' sampletype_names{i} '-wleg.png'];
    else
        fig_filename = [ dir_plots '/' '16S_subj-3_' sampletype_names{i} '.png'];
    end
    make_miniplot( these_samples_name, bar_temp, x_labels, graph_labels, legend_bool, fig_filename )

end

% % Subject G % skipping bcs we do not have WGS colony data for Subject G


