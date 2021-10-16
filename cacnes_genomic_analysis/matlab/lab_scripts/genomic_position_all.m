function [cds_indices,locus_tag_nums] = genomic_position_all(cds_struct,genomelength, ChrStarts)
%

% cds_indices: the index number of the corresponding; returns 1.5 for a 
%   position between CDS{1}(1) and CDS{1}(2) (this is the old behavior)
% locus_tag_nums: the tag number of CDS (e.g. '5' for PAG1_00005); returns 0.5
%   for any intergenic position

%% Version history

% Edited 2019_01 by TDL and FMK to fix error handling multiple contigs with
% overwriting and to initizalize cds_indicies as 0.5 instead of 0

% Edited by TDL 2018_08_13 to return two different arguments-- cds_num was
% only previous output. Also now requires ChrStarts as input

%%
if nargin < 3 & size(cds_struct) > 1
    error('More than 1 chromosome/contig detected in cds_struct but ChrStarts not provided')
elseif size(cds_struct) == 1
    ChrStarts=0;
end

locus_tag_nums=0.5*ones(genomelength,1) ;
cds_indices=0.5*ones(genomelength,1) ;

for c=1:numel(cds_struct)
    
    this_chr=cds_struct{c};
    
    if isempty(this_chr) %Skip contigs w/o CDS annotation.
        continue
    end
    
    cds_num=div_get_gene_numbers(this_chr);
    
    genestarts=ChrStarts(c)+[this_chr.loc1];
    geneends=ChrStarts(c)+[this_chr.loc2];
    
    for i=1:(numel(genestarts)-1)
        locus_tag_nums(genestarts(i):geneends(i))=cds_num(i);
        locus_tag_nums(geneends(i)+1:genestarts(i+1)-1)=0.5;
        
        cds_indices(genestarts(i):geneends(i))=i;
        cds_indices(geneends(i)+1:genestarts(i+1)-1)=i+0.5;
    end

    
    if geneends(i+1) < genestarts(i+1)
       %if last gene overlaps with origin, discount part from origin to
       %gene end
        if c < numel(ChrStarts)
            geneends(i+1)=ChrStarts(c+1);  
        else
            geneends(i+1)=genomelength; 
        end
    end
    
    locus_tag_nums(genestarts(i+1):geneends(i+1))=cds_num(i)+1;
    cds_indices(genestarts(i+1):geneends(i+1))=i+1;
    
    if c < numel(ChrStarts)
        cds_indices(geneends(i+1):ChrStarts(c+1))=i+1.5;
    else
        cds_indices(geneends(i+1):genomelength)=i+1.5;
    end
        
end



return
end