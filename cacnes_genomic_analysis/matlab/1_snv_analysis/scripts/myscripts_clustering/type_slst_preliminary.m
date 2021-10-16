function [ slst_all, SampleNamesLong_all ] = type_slst_preliminary( ...
    dir_ref_genome, p, maNT_all, NTs, SampleNames_all )


%% Read in genome information
fprintf(1,'Reading in reference genome...\n')

[ChrStarts, GenomeLength, ~, ScafNames]= genomestats(dir_ref_genome);
refnt_all = extract_outgroup_mutation_positions(dir_ref_genome, p2chrpos(p,ChrStarts));

[refnti_all,~]=ismember(refnt_all,NTs);


%% Extract SLST sequences and print to a file
fprintf(1,'Making file for SLST typing...\n')
% SLST = single locus sequencing type
% slst_types.txt is the file generated to use on the web-based tool 

% Note: SLST sequence can vary in length; some of these won't match the
% reference SLST set because they are forced to be the same length as the
% SLST sequence in the reference genome (even though they aren't). These
% will be updated later with data from genome assemblies.

% Start/end of SLST sequence on reference genome
p_slst_start = 1469806;
p_slst_end = 1470289;

% SLST in reference genome
reference_slst = extract_outgroup_mutation_positions(dir_ref_genome, p2chrpos(p_slst_start:p_slst_end,ChrStarts));
[~,ancestra_slsti]=ismember(reference_slst,NTs);
anc_nuc =NTs(ancestra_slsti);

% SLST positions in p
in_slst = find(p>=p_slst_start & p<=p_slst_end);
p_slst = p(in_slst);
char_slst = p_slst-1469805;

% Extract SLST sequences
% Note: fills in non-candidate-SNP positions with reference genome allele
slst_variable_nuc = cell( size(SampleNames_all) ); % initialize
SLST_sequences = cell( size(SampleNames_all) ); % initialize
for i=1:numel(SampleNames_all)
    slst_variable_nuc{i} = NTs(maNT_all(in_slst,i));
    SLST_sequences{i} = anc_nuc;
    SLST_sequences{i}(char_slst) = slst_variable_nuc{i};
end

% Save SLST sequences in a file
fid=fopen('data/slst_sample_seqs.txt','w');
for i=1:numel(SampleNames_all)
    fprintf(fid,['>' SampleNames_all{i} '\n' SLST_sequences{i} '\n']);
end
fclose(fid);


%% Pause so that user can get SLST types from web tool
% Use webpage for typing: http://medbac.dk/slst/pacnes

fprintf(1,'Used http://medbac.dk/slst/pacnes to get preliminary sample SLSTs...\n')


%% Load SLST types

% Load SLST (obtained via web tool)
slst_all = readtable('data/slst_sample_types.csv','Delimiter',',');
slst_all = slst_all.Type;


%% Update SampleNames to include SLST type

SampleNamesFormattedSLST = cell( size(SampleNames_all) );
for i=1:numel(SampleNames_all)
    SampleNamesFormattedSLST{i} = [ SampleNames_all{i} '_SLST-' slst_all{i} ];
end
SampleNamesLong_all = SampleNamesFormattedSLST;


end