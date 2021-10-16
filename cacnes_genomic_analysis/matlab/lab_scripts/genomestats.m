function [CStarts, GLength, CIndicator, ScafNames, Sequences]= genomestats(REF_GENOME_DIRECTORY)

% Comment: Note that indexing for genome positions in CStarts begins at 0,
% NOT 1!!!!

fr = fastaread([REF_GENOME_DIRECTORY '/genome.fasta']) ;

GLength=0;
CStarts=[];
Sequences={};

ScafNames = {fr.Header} ;
for i=1:length(ScafNames)
    f=find(ScafNames{i}==' ',1) ;
    if f>0
        ScafNames{i} = ScafNames{i}(1:f-1) ;
    end
    CStarts(end+1)=GLength;
    GLength=GLength+numel(fr(i).Sequence);
    Sequences{end+1}=fr(i).Sequence;

CIndicator = find(fr(1).Header=='.',1) - 1; %Assumes fewer than 10 chromosomes

end

end
