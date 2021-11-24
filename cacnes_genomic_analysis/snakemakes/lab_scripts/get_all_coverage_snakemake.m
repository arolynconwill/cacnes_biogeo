function [all_coverage_per_bp, cov_modes, all_maf_per_bp] = ...
    get_all_coverage_snakemake( path_to_sample_names_file, path_to_list_of_diversity_files, path_to_coverage_matrix_file ) 

%% Summary

% INPUT: 
% path_to_sample_names_file: text file with list of sample names
% path_to_list_of_diversity_files: text file with list of paths to diversity.mat's
% path_to_coverage_matrix_file: string with the path and filename for
% saving the coverage matrix relative to snakemake's working directory

% OUTPUT: where n is number of samples 
% all_coverage_per_bp: [n x GenomeLength] coverage at each bp
% cov_modes: [n x 1] coverage mode for each sample  

% SNAKEMAKE: example rule for implementing this as part of the standard
% Lieberman Lab case step snakemake
%
% # generate coverage matrix (matlab version)
% rule coverage_matrix_mat:
%     input: 
%         string_sampleID_names = "1-temp_pos/string_sampleID_names.txt",
%         string_diversity_mat = "1-temp_pos/string_diversity_mat.txt",
%     output:
%         coverage_mat = "2-candidate_mutation_table/coverage_matrix.mat",
%     log:
%         'logs/coverage_matrix.log',
%     shell:
%         """
%         module add mit/matlab/2015b; 
%         matlab -r "path('{SCRIPTS_DIRECTORY}',path); get_all_coverage_snakemake( '{input.string_sampleID_names}', '{input.string_diversity_mat}', '{output.coverage_mat}' )" > {log} 2>&1
%         """


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


%% Files: diversity.mat for each sample
fprintf(1,'Processing diversity.mat filenames...\n');

% Import list of directories for where to quals.mat for each sample
fname = [ pwd '/' path_to_list_of_diversity_files ];
fid = fopen( fname, 'r' );
fstring = fgetl(fid); % from echo, just one line
paths_to_diversity_files = strsplit(fstring,' ');
fclose(fid);


%% Coverage matrix
fprintf(1,'Generating coverage matix...\n');

% Loop through samples to fill in matrix rows
for i = 1:numSamples
    
    filename = paths_to_diversity_files{i};
    
    % Get coverage for this sample
    fprintf(1,['Getting coverage for isolate ' num2str(i) ' of ' num2str(numSamples) '...\n'])
    [coverage, coverage_mode, maf] = get_sample_coverage_snakemake( filename );  
    
    % Initialize variables (if not done yet)
    if i==1
        genome_size = length(coverage);
        all_coverage_per_bp = uint16(zeros( numSamples, genome_size ));
        all_maf_per_bp = uint16(zeros( numSamples, genome_size ));
        cov_modes = zeros( numSamples,1 ); 
    end
    
    % Record data for this sample
    all_coverage_per_bp(i,:) = coverage;
    all_maf_per_bp(i,:)=1000*maf;
    cov_modes(i) = coverage_mode(1); 
    
end


%% Save data
fprintf(1,'Saving coverage matix...\n');

save([ pwd '/' path_to_coverage_matrix_file ], 'SampleNames', 'all_coverage_per_bp', 'all_maf_per_bp', 'cov_modes', '-v7.3') ;
fprintf(1,['Filename: ' pwd '/' path_to_coverage_matrix_file '\n']) ;

clear all


end