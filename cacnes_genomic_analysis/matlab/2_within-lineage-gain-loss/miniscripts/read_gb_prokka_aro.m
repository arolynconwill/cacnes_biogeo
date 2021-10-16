function read_gb_prokka_aro(RefGenomeDirectory)

% May need revisions for your reference genome

%% Revision Historys

%Edited by TDL 2018_08_13 to not save genenum -- can be generated quickly
% with genomic_position_all
% FMK: CDS generation using genbankFile.Sequence if exists rather than full
% genome. Previously led to wrong seq/translation for chr2 and onwards.
% 2019/11 JSB and TDL changed extractFeature sub-function of
% featurelocation div_add_in_nonCDS
% JSB winter 20/21 -- made to work with prokka


%% Intialize

% Load fasta
fastaFile = [RefGenomeDirectory '/genome.fasta'];
fastaFile = fastaread(fastaFile);

genbankAll = genbankread([RefGenomeDirectory '/prokka_out.gbk']); % ARO: added /

[~, ~, ~, ScafNames]= genomestats(RefGenomeDirectory);


%%

fprintf(1,'Reading in gb file...\n'); % ARO: removed tic; no toc

clear CDS

for i=1:length(ScafNames)
    if numel(genbankAll(:)) ~= length(ScafNames)
        error('scaffold number must equal the number of scaffolds in the genbank file. ')
    end
    
    genbankFile = genbankAll(i);
    
    % LOAD GENBANK
    if ~isfield(genbankFile,'Sequence')
        fprintf(1,'Warning: Genbank file may have no CDS/Feature or first line of the genbank file may not be at least 79 characters including whitespaces\n');
        genbankFile.Sequence = lower(fastaFile(i).Sequence);
    elseif isempty(genbankFile.Sequence)% check if sequence is empty
        genbankFile.Sequence = lower(fastaFile(i).Sequence);
    end
    
    if isfield(genbankFile,'CDS') || isfield(genbankFile,'Features')
        genes = locustag_from_text(genbankFile.CDS) ;
        genes = div_add_in_nonCDS(genes, genbankFile.Features);
        
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


%% Save CDS

save([RefGenomeDirectory '/cds_sorted'],'CDS')


end
