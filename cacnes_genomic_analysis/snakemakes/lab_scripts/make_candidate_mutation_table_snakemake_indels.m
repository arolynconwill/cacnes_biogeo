function make_candidate_mutation_table_snakemake( path_to_p_file, path_to_sample_names_file, path_to_outgroup_boolean_file, path_to_list_of_quals_files, path_to_list_of_diversity_files, path_to_candidate_mutation_table )
% Gather everything together for candidate_mutation_table

% Inputs:
% % path_to_p_file: where to find all_positions.mat
% % path_to_sample_names_file: where to find text file with sample names
% (space delimited)
% % path_to_outgroup_boolean_file: where to find text file with outgroup
% booleans (space delimited, 1=outgroup, 0=not)
% % path_to_list_of_quals_files: where to find text file with list of
% quals.mat files for each sample (space delimited)
% % path_to_list_of_diversity_files: where to find text file with list of
% diversity.mat files for each sample (space delimited)
% Output:
% % path_candidate_mutation_table: where to write
% candidate_mutation_table.mat, ex. results/candidate_mutation_table.mat

% Note: All paths should be relative to pwd!


%% Version history

% % This is adapted from TDL's build_mutation_table_master_smaller_file_size_backup.m
% % Arolyn, 2018.12.19: This script was written as part of the transition
% to snakemake. It performs the part of the case step that gathers data for
% Quals and counts and saves candidate_mutation_table.mat
% % Arolyn, 2019.02.12: Added another matlab variable that stores indel
% statistics; called 'indel_counter'.


%% p: positions on genome that are candidate SNPs
fprintf(1,'Processing candidate SNP positions...\n');

load(path_to_p_file,'p') % from previous step, should include variable called p

fprintf(1,['Total number of positions: ' num2str(length(p)) '\n']);


%% SampleNames: list of names of all samples
fprintf(1,'Processing sample names...\n');

% Input is space separated text file??
fname = [ pwd '/' path_to_sample_names_file ];
fid = fopen( fname, 'r' );
fstring = fgetl(fid); % from echo, just one line
SampleNames = strsplit(fstring,' '); 
fclose(fid);

numSamples = length( SampleNames ); % save number of samples

fprintf(1,['Total number of samples: ' num2str(numSamples) '\n']);


%% in_outgroup: booleans for whether or not each sample is in the outgroup
fprintf(1,'Processing outgroup booleans...\n');

% Input is space separated text file??
fname = [ pwd '/' path_to_outgroup_boolean_file ];
fid = fopen( fname, 'r' );
in_outgroup_string = fgetl(fid); % from echo, just one line
fclose(fid);
in_outgroup_cell = strsplit(in_outgroup_string,' ');
in_outgroup = cellfun(@(x) str2num(x), in_outgroup_cell);
% check if this is oriented in the right direction??


%% Quals: quality score (relating to sample purity) at each position for all samples
fprintf(1,'Gathering quality scores at each candidate position...\n');

% Import list of directories for where to quals.mat for each sample
fname = [ pwd '/' path_to_list_of_quals_files ];
fid = fopen( fname, 'r' );
fstring = fgetl(fid); % from echo, just one line
paths_to_quals_files = strsplit(fstring,' ');
fclose(fid);

% Make Quals
Quals = zeros(length(p), numSamples,'int16'); % initialize
for i=1:numSamples
    fprintf(1,'Loading quals matrix for sample: %g  \n',i) ;
    fprintf(1,['Filename: ' paths_to_quals_files{i} '\n']) ;
    load(paths_to_quals_files{i});
    Quals(:,i)=quals(p);
end


%% counts: counts for each base from forward and reverse reads at each candidate position for all samples
fprintf(1,'Gathering counts data at each candidate position...\n');

% Import list of directories for where to diversity.mat for each sample
fname = [ pwd '/' path_to_list_of_diversity_files ];
fid = fopen( fname, 'r' );
fstring = fgetl(fid); % from echo, just one line
paths_to_diversity_files = strsplit(fstring,' ');
fclose(fid);

% Make counts
counts = zeros(8, length(p), numSamples,'uint16'); % initialize
indel_counter=zeros(2, length(p), numSamples, 'uint16'); % Added 2019.02.12
for i=1:numSamples
    fprintf(1,'Loading counts matrix for sample: %g  \n',i) ;
    fprintf(1,['Filename: ' paths_to_diversity_files{i} '\n']) ;
    load(paths_to_diversity_files{i});
    counts(:,:,i)=data(1:8,p);
    indel_counter(:,:,i)=data(39:40,p); % Added 2019.02.12; reads supporting indels and reads supporting deletions
end

%fprintf(1,'Getting all the coverage information...\n');
%[all_coverage_per_bp, ~, all_maf_per_bp] = get_all_coverage(SampleInfo, GenomeLength);


%% Save everything!
fprintf(1,'Saving everything...\n');

save([ pwd '/' path_to_candidate_mutation_table ], 'SampleNames', 'p', 'counts', 'Quals', 'in_outgroup', 'indel_counter', '-v7.3') ;
fprintf(1,['Filename: ' pwd '/' path_to_candidate_mutation_table '\n']) ;

fprintf(1,'WARNING! This version does not generate a coverage matrix.\n');
%save('coveragematrix', 'all_coverage_per_bp', '-v7.3'); 

