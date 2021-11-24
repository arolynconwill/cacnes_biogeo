%% Parameters  -- IMPORTANT TO ADJUST !!!!
clear all; close all;
% Much of this will need to vary according to your reference genome,
% coverage, and particular samples
global refgenome

%{
% A:017 subtree creation
load('A_017_candidate_mutation_table_full.mat')
refgenome='Cutibacterium_granulosum_A_017_assembled';
samples_with_contamination = SampleNames(~contains(SampleNames,["A:017-Ext-crMoB-4_t00", "A:017-Ext-crMoB-3_t00", "A:017-Ext-crMoB-2_t00", "A:280-Ext-crMoW-1_t08", "A:028-Ext-tlMoWd-3_t00"])); 
%}

% A:018 subtree creation
load('A_018_candidate_mutation_table_full.mat')
refgenome='Cutibacterium_granulosum_A_018_assembled';
samples_with_contamination = SampleNames(~contains(SampleNames,["A:288-Scr-clCh-4_t08", "A:018-Ext-crMoBl-4_t00", "A:018-Ext-crMoBl-3_t00", "A:018-Ext-crMoBl-2_t00", "A:018-Ext-crMoBl-1_t00", "A:259-Stp-Lc82-2_t08", "A:023-Ext-tcChBp-4_t00", "A:258-Stp-Lc82-1_t08"])); 

%{
% B:087 subtree creation
refgenome='Cutibacterium_granulosum_B_087_assembled';
load('B_087_candidate_mutation_table_full.mat')
samples_with_contamination = SampleNames(~contains(SampleNames,['B:087',"B:076-Scr-blBa-2_t00"])); 
%}

%% Presets
min_average_coverage_to_include_sample = 5;

%Filter on a per-call basis
min_maf_for_call = .75; % maf = major allele frequency
min_cov_per_strand_for_call = 3;% cov = coverage
min_qual_for_call = 30;  % qual = quality from samtools

%Filter on a by-position basis (across samples)
max_fraction_ambigious_samples = .34; % .2 % ambiguous = a sample that failed one filter at that position
min_median_coverage_position = 12; % note: related to min_average_coverage_to_include_sample
promotersize=250;

%filters for using once we have a trusted position
min_maf_for_analysis=.67;
dirty_sample_cutoff = .05; %maximum number of N's in a sample after all other filtering
% corresponds to about .05 of the relevant positions

% Remove all samples that do not contain a desired clade - this is a
% holdover from the larger tree 
%samples_with_contamination = ["A:023-Ext-tcChBp-1_t00",'A:281-Ext-crMoW-1_t08','B:088-Ext-ccChBs-4_t00','B:303-Scr-Rc-3_t11','E:096-Scr-Ch-3_t00','J:327-Scr-Ch-4_t11',"J:327-Scr-Ch-6_t11"];
disp('Params set')
%% Enviornment set up -- probably won't need to change

workingdir=char(pwd);
REFGENOMEFOLDER=['~/Dropbox (MIT)/Lieberman Lab/Reference_Genomes/' refgenome];
SCRIPTSDIRECTORY = ['~/Dropbox (MIT)/Lieberman Lab/scripts'];
path(SCRIPTSDIRECTORY,path);

NTs='ATCG';

%% Remove undesired samples based on name and/or coverage
coverage=squeeze(sum(counts(1:8,:,:)));

goodsamples = mean(coverage) > min_average_coverage_to_include_sample & ~contains(SampleNames,samples_with_contamination);

SampleNames_orig = SampleNames;
SampleNames=SampleNames(goodsamples);
counts=counts(:,:,goodsamples);
Quals=Quals(:,goodsamples);
coverage=coverage(:,goodsamples);
num_samples=numel(SampleNames);

coverage_forward_strand=squeeze(sum(counts(1:4,:,:)));
coverage_reverse_strand=squeeze(sum(counts(5:8,:,:)));

Quals = -1*Quals; %use -Quals because this way higher numbers are more confident

%% Read in genome information

[ChrStarts, GenomeLength, ~, ScafNames]= genomestats(REFGENOMEFOLDER);
refnt = extract_outgroup_mutation_positions(REFGENOMEFOLDER, p2chrpos(p,ChrStarts));

% Is be more complicated if refnt ~= ancnt
ancnt=refnt;
[~,ancnti]=ismember(refnt,NTs);
ancnti_m=repmat(ancnti,1,num_samples);

%% Make some basic structures for finding mutations

contig_positions=p2chrpos(p, ChrStarts);
[maf, maNT, minorNT, minorAF] = div_major_allele_freq(counts);
[~,refnti]=ismember(refnt,NTs);

mutantAF=zeros(size(maNT));
mutantAF(maNT~=ancnti_m)=maf(maNT~=ancnti_m);
mutantAF(minorNT~=ancnti_m)=mutantAF(minorNT~=ancnti_m)+minorAF(minorNT~=ancnti_m); %this construction allows for positions with two different mutations

disp('Basic structures made')
%% Find positions with polymorphic mutations

%diversemutation=div_test_thresholds(counts,STRICT_PARAMETERS, coveragethresholds);
%Calls(diversemutation>0)=-1;  %-1 records that a diverse mutation was called
diversemutation=zeros(size(maNT));

%% Find positions with fixed mutations

% Find the nucleotide identity at each position.
Calls=maNT;
Calls(abs(Quals) < min_qual_for_call | maf< min_maf_for_call | coverage_forward_strand < min_cov_per_strand_for_call | coverage_reverse_strand < min_cov_per_strand_for_call)=0; %Quals < min_qual_for_call |

Calls(sum(Calls<1,2)>=(num_samples*max_fraction_ambigious_samples) | median(coverage,2)<min_median_coverage_position,:)=0;

[MutQual, MutQualIsolates] = ana_mutation_quality(Calls,Quals) ; 

fixedmutation=((maNT~=repmat(ancnti,1,num_samples)) & maNT>0 & repmat(MutQual,1,num_samples)>=1 & maf> min_maf_for_analysis);

hasmutation= fixedmutation | diversemutation;
goodpos=find(sum(fixedmutation,2)>0);

disp('Finding SNPs')

%% Display table (with or without annotation)

positions_to_show_in_table=goodpos;
samples_to_show=1:numel(SampleNames);
num_contigs=max(contig_positions(:,1));
contig_lengths=[ChrStarts(2:end) GenomeLength]-ChrStarts;

% Whether or not to sort the clickable table by quality
QualSort=0; % toggle to sort the table by quality or not
SpecificSamples = 0; %manually input a sample list to observe s
SpecificSite = 0; %manualy select a site NOT WORKING 

sitetoparse = positions_to_show_in_table;
site_to_view_logic = (calls_for_treei(:,26)~=calls_for_treei(:,27)) & (calls_for_treei(:,27)~=0) & (calls_for_treei(:,26)~=0);
sitetoparse = positions_to_show_in_table(site_to_view_logic);
%{
% %if SpecificSite
%     %Z = find(822703)
%     %sitetoparse = Z{1}
% %end

%randomorder=randperm(numel(SampleNames));
%samples_to_show=randomorder(1:10);
%samples_to_show=[samples_to_show'; find(ismember(SampleNames,{'P02LT10504';'P02LT-3903';'P02LT-2412';'P02LT-2407';'P02LT-2403';'P02LT-3513';'P02LT-3510';'P02LT-0602'}))];
%}

% Uncomment the following three lines if there is an annotated genome

%annotations = annotate_mutations_gb(p2chrpos(p(positions_to_show_in_table),ChrStarts),REFGENOMEFOLDER) ;
%annotation_full= append_annotations(annotations, ancnti(positions_to_show_in_table), Calls(positions_to_show_in_table,:), counts(:,positions_to_show_in_table,:), hasmutation(positions_to_show_in_table,:), promotersize) ; %% adds information about particular mutations observed, based on
%clickable_snp_table(annotation_full(site_to_view_logic), Calls(sitetoparse,samples_to_show), counts(:,sitetoparse,samples_to_show), SampleNames(samples_to_show), ScafNames, MutQual(positions_to_show_in_table), QualSort);
%clickable_snp_table(annotation_full(site_to_view_logic), calls_for_treei(sitetoparse,samples_to_show), counts(:,sitetoparse,samples_to_show), SampleNames(samples_to_show), ScafNames, MutQual(positions_to_show_in_table), QualSort);

hasmutation=hasmutation(goodpos,:);

chrpos=p2chrpos(p,ChrStarts);

%Creating table to check for chromosome/plasmid ejection
%{
chromosome_filtered_counts = counts(:,positions_to_show_in_table,samples_to_show);
[size_fil_z,size_fil_x,size_fil_y] = size(chromosome_filtered_counts);
compressed_counts = zeros(size_fil_x,size_fil_y);

for z=1:size_fil_z
    for x = 1:size_fil_x
        for y = 1:size_fil_y
            compressed_counts(x,y) = compressed_counts(x,y) + chromosome_filtered_counts(z,x,y);
        end
    end
end

%filtered_positions = p(positions_to_show_in_table);
%Filteredcalls = Calls(positions_to_show_in_table,samples_to_show);
%bar(filtered_positions,compressed_counts)

%}
%% Finding potential contaminants

maf_good=maf(positions_to_show_in_table,:);
maf_mean=mean(maf_good);
maf_outliers=sum(maf_good<.9);

 %% parsimony tree

disp('Making tree')

samplestoplot=1:numel(SampleNames);
 
quality_positions=goodpos;
 
calls_for_treei=maNT(quality_positions,samplestoplot);
maf2 = maf(quality_positions,samplestoplot);
calls_for_treei(maf2 < min_maf_for_analysis)=0;

%checks for samples with too many N's
badsamples = sum(calls_for_treei==0);
indices_dirty_samples = find(badsamples>dirty_sample_cutoff*size(maf2,1));

samplestoplot(indices_dirty_samples) = [];
%reruns code to reform calls_for_treei

calls_for_treei=maNT(quality_positions,samplestoplot);
maf2 = maf(quality_positions,samplestoplot);
calls_for_treei(maf2 < min_maf_for_analysis)=0;

calls_for_tree=zeros(size(calls_for_treei));
calls_for_tree(calls_for_treei > 0)=NTs(calls_for_treei(calls_for_treei>0));
calls_for_tree(calls_for_treei < 1)='?';

outgroup_nts = extract_outgroup_mutation_positions(REFGENOMEFOLDER, p2chrpos(p(quality_positions),ChrStarts));
 
TreeSampleNames= {SampleNames{samplestoplot}};
 
calls_for_tree = [outgroup_nts, calls_for_tree];
TreeSampleNames={'Reference' SampleNames{samplestoplot}};
 

% parsimony
[treefilename, tempnames] =generate_parsimony_tree_rename(calls_for_tree, TreeSampleNames, refgenome);
 
fprintf(1,'\nDone with tree\n');
 
%relabeltree(treefilename,newtreefilename,tempnames, newnames);
 
%% Make a tree for each SNP location
%{
%samplestoplot=samplestoplot(~in_outgroup);
 
samplestoplot=1:numel(SampleNames);
 
if exist('tree_counting','dir')
    rmdir('tree_counting','s') 
end
mkdir('tree_counting')
cd('tree_counting')
 
fid=fopen('for_tree_labeling.csv','w');
fprintf(fid,'chr,pos');
for i=1:numel(TreeSampleNames)
    fprintf(fid,[',' TreeSampleNames{i}]);
end
for i=1:numel(goodpos)
    fprintf(fid,['\n' num2str(contig_positions(positions_to_show_in_table(i),1)) ',' num2str(contig_positions(positions_to_show_in_table(i),2))]);
    for j=1:size(calls_for_tree,2)
        fprintf(fid,[',' calls_for_tree(i,j)]);
    end
end
fprintf(fid,'\n');
fclose(fid);
eval(['! python2.7 '  SCRIPTSDIRECTORY '/countMutations.py ../' treefilename ' for_tree_labeling.csv'])
cd('..')
%}