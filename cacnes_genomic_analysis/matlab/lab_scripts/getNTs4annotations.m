function [NTs4annotations] = getNTs4annotations(Calls,p_all,Nsample,goodpos,outgroup_isolates,refnt_all)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JSB and ARO 18 May 2020
% this function generates a list of positions and
% corresponding nucleotides of the anscetral alleles for a given cluster
% in order to properly call mutations with annotate_mutations_gb.m
% instead of identifying a mutation with respect to a reference genome,
% these ps and NTs will update the reference sequence such that a mutation
% is called from the genome of the MRCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check to make sure the right version of Calls is passed in


if numel(Calls(:,1)) ~= numel(p_all) || numel(Calls(1,:)) ~= Nsample
    error('Calls should have the same number of rows as p_all and the same number of columns as Nsample. Ensure that you are passing the version of calls which includes non-variable positions i.e. before Calls=Calls(~nonvariablep,:);');
end


%% Initialize output arrays


[NTs4annotations,isgoodpos] = deal(zeros(numel(p_all),1));

isgoodpos(goodpos) = 1; % turn goodpos into a boolean array with the same size as p_all

%% Determine which positions are nonvariable within the cluster


nonvariableNTs=(sum(Calls==1,2)==Nsample | sum(Calls==2,2)==Nsample | sum(Calls==3,2)==Nsample | sum(Calls==4,2)==Nsample);

% add these Nts to NTs4annotations using the first column of calls, since they
% are invariant across samples

NTs4annotations(nonvariableNTs) = Calls(nonvariableNTs,1);


%% Define ancestor for positions calling a function which treats positions in goodpos and ~goodpos differently

% for positions not in goodpos, the function will start with the mode of
% the ingroup calls. If it is not in goodpos, it will start with the mode
% of the outgroup calls

%calls the function, called anscestor
NTsvariablegoodpos = anscestor(isgoodpos&~nonvariableNTs);
NTsvariableNOTgoodpos = anscestor(~isgoodpos&~nonvariableNTs);

%assigns the nucleotides of the output based on the function called above
NTs4annotations(isgoodpos&~nonvariableNTs)= NTsvariablegoodpos(isgoodpos&~nonvariableNTs);
NTs4annotations(~isgoodpos&~nonvariableNTs) = NTsvariableNOTgoodpos(~isgoodpos&~nonvariableNTs);

%%% this function has a different order of operations for goodpos positions
%%% and those which are not goodpos
    function anc_nti = anscestor(bool)
    % initialize output array
    anc_nti = nan(numel(bool),1); 
    
    outgroup_calls = Calls(:,outgroup_isolates);
    outgroup_calls( outgroup_calls==0 ) = NaN;  %mode takes the lowest value when there is a tie, setting to nan avoids this
    outgroup_calls = mode(outgroup_calls,2); 
    
    ingroup_calls = Calls(:,~outgroup_isolates);
    ingroup_calls(ingroup_calls==0) = NaN;  %mode takes the lowest value when there is a tie, setting to nan avoids this
    ingroup_calls = mode(ingroup_calls,2);
    
    if any(bool&isgoodpos) % this checks if you are looking at goodpos
        % 1 outgroup calls
        anc_nti = outgroup_calls;
        
        % 2 reference genome calls
        anc_nti(isnan(anc_nti)) = refnt_all(isnan(anc_nti));
        
        % 3 ingroup calls
        %anc_nti(anc_nti==0) = ingroup_calls(anc_nti==0);
        % not implemented since reference should always have a call

    else 
        % 1 ingroup calls
        anc_nti = ingroup_calls; 

        % 2 outgroup calls
        anc_nti(isnan(anc_nti)) = outgroup_calls(isnan(anc_nti));

        % 3 replace with reference if nothing else
        anc_nti(isnan(anc_nti)) = refnt_all(isnan(anc_nti)); % Fixed (Aro, 2020.07.02)

    end


    end
end

