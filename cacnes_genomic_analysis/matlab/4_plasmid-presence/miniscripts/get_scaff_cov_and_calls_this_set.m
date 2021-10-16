function [ Calls, coverage, p, counts, Quals ] = get_scaff_cov_and_calls_this_set( filename, SampleNames_this_set, ...
    min_qual_for_call, min_maf_for_call, min_cov_each_strand_for_call )

% Load data
load( filename )

% Find samples from this set in candidate mutation table
[ in_set, index_in_table ] = ismember( SampleNames_this_set, SampleNames );

% Downsize to samples with this scaffold type
counts = counts(:,:,index_in_table);
Quals = Quals(:,index_in_table);
SampleNames = SampleNames(index_in_table);
% ignoring in_outgroup and indel_counter because not using currently

% Get coverage from counts
coverage = squeeze( sum( counts ) );
cov_fwd_strand=squeeze(sum(counts(1:4,:,:)));
cov_rev_strand=squeeze(sum(counts(5:8,:,:)));

% Compute major allele and major allele frequency
[maf, maNT, ~, ~] = div_major_allele_freq( counts );

% Switch sign of quals
Quals = -1*Quals;

% Get loose calls for each position using filters provided in input
Calls = maNT;
Calls( Quals < min_qual_for_call ...
    | maf < min_maf_for_call ...
    | cov_fwd_strand < min_cov_each_strand_for_call ...
    | cov_rev_strand < min_cov_each_strand_for_call ...
    ) = 0; % set to N

end