clear all; close all;

fid=fopen('dna-sequences.fasta','r');

alternate_textscan = textscan(fid,'%s','delimiter','\n');
alternate_textscan_decompressed = alternate_textscan{:};
number_lines = length(alternate_textscan_decompressed);

genuses_to_clean = {{'cutibacterium', 'propionibacterium'}, 'staphylococcus', {'corynebacterium','lawsonella'},'neisseria'};
genuses_to_clean_file_demon = ["cutibacterium", "staphylococcus","corynebacterium_lawsonella", "neisseria"];
genus_fasta_files = repmat({zeros(number_lines,1)},size(genuses_to_clean));

NOT_specified_genus_fasta_temp = zeros(number_lines,1);
confusing_genus_fasta = zeros(number_lines,1);

confusing_species_index = cellfun(@(x) sum(contains_genus(lower(x),genuses_to_clean))>1, alternate_textscan_decompressed);

parfor ii = 1:length(genuses_to_clean)
    contains_genus_name = contains(lower(alternate_textscan_decompressed), genuses_to_clean{ii});
    contains_genus_name(confusing_species_index) = 0;
    indexes_of_sequence = find(contains_genus_name) + 1;
    contains_genus_name(indexes_of_sequence) = 1;
    genus_fasta_files(ii) = {contains_genus_name};
    NOT_specified_genus_fasta_temp = NOT_specified_genus_fasta_temp + contains_genus_name;
end

fasta_as_string = cellfun(@string, alternate_textscan_decompressed, 'UniformOutput',false);
NOT_specified_genus_fasta = NOT_specified_genus_fasta_temp==0;
NOT_specified_genus_fasta(confusing_species_index) = 0;
NOT_specified_genus_fasta(find(confusing_species_index)+1) = 0;


for ii = 1:length(genuses_to_clean)
    newfilename = "SILVA_v132_" + genuses_to_clean_file_demon(ii) + "_unmodified.fasta";
    fid = fopen(newfilename, 'wt');    
    fprintf(fid, '%s\n', [fasta_as_string{genus_fasta_files{ii}}]);
    fclose(fid);
end

fid = fopen('SILVA_v132_not_specified_genus_unmodified.fasta', 'wt');
fprintf(fid, '%s\n', [fasta_as_string{NOT_specified_genus_fasta}]);
fclose(fid);

confusing_species_index(find(confusing_species_index)+1) = 1;
fid = fopen('SILVA_v132_confusing_genus_unmodified.fasta', 'wt');
fprintf(fid, '%s\n', [fasta_as_string{confusing_species_index}]);
fclose(fid);
