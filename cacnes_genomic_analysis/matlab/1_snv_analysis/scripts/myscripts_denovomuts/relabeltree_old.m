function relabeltree_old(treefile,newtreefile,oldnames,newnames)

%% Version history:
% Arolyn, 2018.10.05: 
% % need char to write file... fprintf(fout,char(newtree));
% % fixed weird spacing after warning ": " from ": "
% Arolyn, 2018.11.15:
% % added wait time in case treefile doesn't exist yet...
%%

%%  Purpose: To relabel the tips of trees in a Newick-formatted file.
%Especially for using with phylip programs, which require tip labels which
%are <=10 characters long

%Tami Lieberman 2016

%%

% Check to make sure file exists; otherwise wait
% wait = true;
% while wait
%     if exist(treefile)
%         wait = false;
%     else
%         fprintf(1,['Waiting 5 seconds for ' treefile '\n']);
%         pause(5); % wait 5 seconds
%     end
% end

fin=fopen(treefile,'r');
fout=fopen(newtreefile,'w');

x=fread(fin)';
newtree=x;

for i=1:numel(oldnames)
    j=strfind(char(newtree),[oldnames{i} ':']);
    oldnamelength=numel(oldnames{i});
    if ~isempty(j)
        if numel(j)==1 %& j < numel(x)-30
            newtree=[newtree(1:j-1) newnames{i} newtree((j+oldnamelength):end)];
        else
            error(['two instances of oldname: ' oldnames{i}]);
        end
    else
        fprintf(1,['Warning: No sample found in treefile with name: ' oldnames{i} '\n']);
    end
end

fprintf(fout,char(newtree));

fclose(fin);
fclose(fout);
        
            
        

