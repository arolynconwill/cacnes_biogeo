function [MutQual, MutQualIsolates] = ana_mutation_quality(Call,Qual)
%
% [MutQual, MutQualIsolates] = ana_mutation_quality(Call,Qual)
%
% This function calls mutations within the data itself, instead of w.r.t. 
% the reference genome and calculates the new quality of those calls. It 
% takes as input the called nucleotides (Call) and the quality of the call 
% (Qual) and outputs only the mutation call quality (MutQual) calculated as 
% max_i(min_j(Qual_i, Qual_j)) for i,j over all samples with different Call
% entries. It also outputs the samples giving this maxmin quality 
% (MutQualIsolates).
% 
% If there is no within data mutation in a given position (all nucleotides
% are equal, but different from the reference), MutQual returns 0 at that
% position.

%% Version History
% November 2012: added variable strains so that an outgroup can be removed
% Arolyn,   2019.01.17: Added comments, no changes to code
% GTV,      2020.03.24: Added description of function

%%

[Nmuts, NStrain] = size(Call) ;



MutQual = zeros(Nmuts,1) ; 
MutQualIsolates = zeros(Nmuts,2); 

for k=1:Nmuts
    if (length(unique([Call(k,:), 0]))<=2) %if there is only one type of non-N call, skip this location
        MutQual(k) = nan ;
        MutQualIsolates(k,:) = 0; 
    else
        c=Call(k,:) ; c1=c(ones(1,NStrain),:) ; c2=c1' ;
        q=Qual(k,:) ; q1=q(ones(1,NStrain),:) ; q2=q1' ;
        g = (c1~=c2 & c1~=0 & c2~=0) ; % boolean matrix identifying pairs of samples where calls disagree (and are not N) at this position
        positive_pos = find(g); 
        
        % get MutQual + logical index for where this occurred
        [MutQual(k), MutQualIndex] = max(min(q1(g),q2(g))) ;%min(q1(g),g2(g)) gives lower qual for each disagreeing pair of calls, we then find the best of these
        % store which strains were used to call this mutation 
        [strain_i, strain_j] = ind2sub(size(g), positive_pos(MutQualIndex));
        MutQualIsolates(k,:) = [strain_i, strain_j]; 
        
    end
end

MutQual(isnan(MutQual))=0;

return
end
