function [coverage, coverage_mode, MAF] = get_sample_coverage_snakemake( filename ) 

    % Summary: this script loads diversity.mat and returns information on
    % coverage and major allele frequency (MAF)

    diversity_file = strcat(filename); 
    
    if exist(diversity_file, 'file')
        
        diversity = load(diversity_file);
        countsdata = diversity.data;
    
        % calculate coverage
        coverage = sum(countsdata(1:8,:));
        average_coverage = mean(coverage); 
        coverage_mode=mode(coverage);
        
        % calculate major allele frequency
        c=countsdata(1:4,:)+countsdata(5:8,:); % counts for each base at each position
        [sorted, sortedpositions] = sort(c,1);
        maxcount = sorted(end,:,:); % number of counts for major allele at each position
        MAF=double(maxcount)./sum(c,1);
        MAF(isnan(MAF))=0; %set to 0 to indicate no data

    else
        
        error(['No diversity.mat file found at :' diversity_file])
    
    end
   
end