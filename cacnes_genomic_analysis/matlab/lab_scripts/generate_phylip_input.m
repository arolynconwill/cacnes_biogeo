function tempnames = generate_phylip_input(calls, names, fileout)

%Prepares input file. Changes names temporarily if names are too long

%Print to screen if output file not specified
if nargin < 3
    fid=1;
else
    fid = fopen(fileout,'w');
end


%Short temporary name
tempnames=names;
if max(cellfun(@length,names)) > 9
    for i=1:numel(names)
        tempnames(i)={['T^' num2str(i)]};
    end
end



NStrain=size(calls,2);

%phylip header
fprintf(fid, ['   ' num2str(NStrain) '   ' num2str(size(calls,1)) '\n']);

%phylipinput=char(10+size(printCalls,1),size(printCalls,2));

for i=1:NStrain
    name=tempnames{i}; 
    
    for j=1:(10-numel(name))
        name=[name ' '];
    end
    
    fprintf(fid, [name '   ' calls(:,i)' '\n']);
end
end