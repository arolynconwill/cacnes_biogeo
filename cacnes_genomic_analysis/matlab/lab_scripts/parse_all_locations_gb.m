function gb = parse_all_locations_gb(old_gb, Sequence) 

%% Version History

% Arolyn, 2021.02.24: using bacterial genetic code when populating translation2 field. 

%% Code

gb = old_gb ;
notgenes=[];


for i=1:length(gb)
    f = parse_location_gb(gb(i).location) ;
    if isempty(f.loc1)
        notgenes(end+1)=i;
    else
        gb(i).loc1 = f.loc1 ;
        gb(i).loc2 = f.loc2 ;
        %there seems to be a mess up somewhere above, reverse strands here
        if f.strand
            gb(i).strand = 0;
            gb(i).Sequence = Sequence(f.loc1:f.loc2) ;
        else
            gb(i).strand = 1;
            gb(i).Sequence = seqrcomplement(Sequence(f.loc1:f.loc2)) ;
        end
    end
    gb(i).translation2=nt2aa(gb(i).Sequence,'ACGTOnly',0,'GeneticCode',11); % Arolyn, 2021.02.24: specified bacterial genetic code
end

gb(notgenes)=[];

return 
end