function [newtreefilename, newnames] = generate_parsimony_tree_old(calls, names, prefix_for_labeled_tree)

%% Revision history:
% 2018.08.15, Arolyn: Updated section at end that relabels trees; added
% optional third argument ('prefix_for_labeled_tree') that prepends a
% string to the labeled tree filename
%%

newnames=names;
%Removes any ':' from final names, which screws up tree visualization
for i=1:numel(newnames)
    if ~isempty(find(newnames{i}==':',1))
        newnames{i}(newnames{i}==':')='^';
    end
    if ~isempty(find(newnames{i}=='|',1))
        newnames{i}(newnames{i}=='|')='_';
    end
    if ~isempty(find(newnames{i}=='/',1))
        newnames{i}(newnames{i}=='/')='_';
    end
    if ~isempty(find(newnames{i}=='(',1))
        newnames{i}(newnames{i}=='(')='_';
    end
end



timestamp=datestr(now, 'yyyy-mm-dd-HH-MM-SS');

%write input file
tempnames = generate_phylip_input_old(calls, newnames, [timestamp '_infile.txt']);

%write parameter file
fid = fopen([timestamp '_optionfile.txt'],'w');
fprintf(fid, [timestamp '_infile.txt\n']);
fprintf(fid, 'f\n');
fprintf(fid, [timestamp '_out.txt\n']);
fprintf(fid, '5\n');
fprintf(fid, 'V\n');
fprintf(fid, '1\n');
fprintf(fid, 'y\n');
fprintf(fid, 'f\n');
fprintf(fid, [timestamp '_out.tree\n\n']);
fclose(fid);

treefilename=[timestamp '_out.tree'];

%run
fid = fopen('temp.sh','w');
fprintf(fid, ['! ../dnapars < ' timestamp '_optionfile.txt > ' timestamp '_outfile.txt\n']);
fclose(fid);

if ~exist('outfile','file')
    fid = fopen('outfile','w');
    fprintf(fid, 'foo');
    fclose(fid);
    fid = fopen('outtree','w');
    fprintf(fid, 'foo');
    fclose(fid);
end


! chmod +x temp.sh
! ./temp.sh

%this relabeles the tree because the function that generates the input file
%uses temporary names
if exist('prefix_for_labeled_tree','var')
    newtreefilename=[prefix_for_labeled_tree '_' treefilename(1:end-5) '.tree'];
else
    newtreefilename=['labeledtree_' treefilename(1:end-5) '.tree'];
end
relabeltree_old(treefilename,newtreefilename,tempnames,newnames);
