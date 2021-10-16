function read_gb_modified(RefGenomeDirectory)


% May need revisions for your reference genome

%%
%Edited by TDL 2018_08_13 to not save genenum -- can be generated quickly
% with genomic_position_all
% FMK: CDS generation using genbankFile.Sequence if exists rather than full
% genome. Previously led to wrong seq/translation for chr2 and onwards.
% 2019/11 JSB and TDL changed extractFeature sub-function of
% featurelocation div_add_in_nonCDS

%% Intialize

% Load fasta
fastaFilename = [RefGenomeDirectory '/genome.fasta'];
fastaFile = fastaread(fastaFilename);

nts='atcg';
rc='tagc';

[ChrStarts, glength, ~, ScafNames]= genomestats(RefGenomeDirectory);

%%
tic; fprintf(1,'Reading in gb file...\n ');

clear CDS
for i=1:length(ScafNames)
    
    ScafNames_i=ScafNames{i};
    
    fprintf(1,[ScafNames_i '...\n']);
    
    % GRAB HEADER FOR GB FILE
    f=find(ScafNames_i=='|',2,'last');
    if f > 1
        fn = ScafNames_i(f(1)+1:f(2)-1) ;
    else
        fn=ScafNames_i; % if header doesn't have NCBI formatted header
    end
    
    % FIND ANNOTATION SOURCE FOR SCAFFOLD SEQUENCE (GB or FASTA)
    gbfilename = [RefGenomeDirectory '/' fn '.gb'];
    if ~exist(gbfilename, 'file')
        gbfilename = [RefGenomeDirectory '/' fn '.gb']; % add by FMK. When many contigs/chr store gb files in separate folder called: single_contig_anno
        %gbfilename = [RefGenomeDirectory '/single_contig_anno/' fn '.gb']; % add by FMK. When many contigs/chr store gb files in separate folder called: single_contig_anno
       
    end
    if exist(gbfilename, 'file')
        % LOAD GENBANK
        genbankFile = genbankread(gbfilename);
        if ~isfield(genbankFile,'Sequence')
            fprintf(1,'Warning: Genbank file may have no CDS/Feature or first line of the genbank file may not be at least 79 characters including whitespaces\n');
            genbankFile.Sequence = lower(fastaFile(i).Sequence);
        elseif isempty(genbankFile.Sequence)% check if sequence is empty
            genbankFile.Sequence = lower(fastaFile(i).Sequence);
        end
        scafSeq = genbankFile.Sequence;
        
        
    else
        error(['Could not find a genebank file named ' gbfilename]);
    end
    
    if isfield(genbankFile,'CDS') | isfield(genbankFile,'Features')
        genes = locustag_from_text(genbankFile.CDS) ;
        genes = div_add_in_nonCDS(genes, genbankFile.Features); %fixes many things that matlab's genbank reader fails to import
        if isfield(genbankFile, 'Sequence')
            CDS{i} = parse_all_locations_gb(genes, genbankFile.Sequence);
        else
            CDS{i} = parse_all_locations_gb(genes, char([fastaFile.Sequence]+32)) ;  %also reverses strands in this function, sorts by position
        end
        %sort by position on genome
        if ~isempty(CDS{i})
            [~,sortedpositions]=sort([CDS{i}.loc1]);
            CDS{i}=CDS{i}(sortedpositions);
        end
    end
    
    
end



%% save

save([RefGenomeDirectory '/cds_sorted'],'CDS')

return
end