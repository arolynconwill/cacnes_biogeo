function ref_ID  = extract_outgroup_mutation_positions(ref_folder, positions)
    %
    % ref_ID  = extract_outgroup_mutation_positions(ref_folder, positions)
    %
    % Looks up positions in the reference genome saved in ref_folder as a
    % fasta file. Extracts the entry, be it a nucleic acid or amino acid,
    % and passes it to ref_ID.
    
    fastafile = [ref_folder '/genome.fasta']; 
    
    fr = fastaread(fastafile) ;
    
    ref_ID=zeros(size(positions,1),1);
    for i=1:size(ref_ID)
        ref_ID(i) = fr(positions(i,1)).Sequence(positions(i,2)); 
    end
    
end