%% Cleaning staph species
fid=fopen('dna-sequences_noMislabel_StaphSpecLvl.taxonomy','r');
staph_clean_felix = textscan(fid,'%s','delimiter','\n');
staph_clean_felix_decompressed =  staph_clean_felix{:};

staph_lines_to_keep = contains(lower(staph_clean_felix_decompressed), 'staphylococcus;staphylococcus');
staph_ids_taxa_to_keep  = split(staph_clean_felix_decompressed(staph_lines_to_keep), '	');
staph_ids_to_keep = staph_ids_taxa_to_keep(:,1);
staph_taxa_to_keep = staph_ids_taxa_to_keep(:,2);

fid = fopen("final_SATIVA_cleaned_files/final_SATIVA_cleaned_staph.tax", 'wt');
fprintf(fid, '%s\n', string(staph_clean_felix_decompressed(staph_lines_to_keep)));
fclose(fid);

fid2=fopen('SILVA_v132_staphylococcus_unmodified.fasta','r');
staph_clean_Alex = textscan(fid2,'%s','delimiter','\n');
staph_clean_Alex_decompressed =  staph_clean_Alex{:};
headerlines = contains(staph_clean_Alex_decompressed, '>');

staph_Alex_seqs = staph_clean_Alex_decompressed(~headerlines);
staph_Alex_seqs = cellfun(@(x) x(1:180), staph_Alex_seqs, 'UniformOutput', false);
staph_clean_Alex_headers = staph_clean_Alex_decompressed(headerlines);

staph_Alex_IDs = split(staph_clean_Alex_headers, {' Bacteria;',' Eukaryota'});
staph_Alex_ID_only = erase(staph_Alex_IDs(:,1), '>');

sequences_to_pull = cellfun(@(x) find(strcmp(staph_Alex_ID_only, x)), staph_ids_to_keep);
sequences_to_pull_order = sort(sequences_to_pull);

final_fasta_file = strings(length(sequences_to_pull_order)*2,1);
final_fasta_file(1:2:end) = ">" + string(staph_Alex_ID_only(sequences_to_pull_order));
final_fasta_file(2:2:end) = staph_Alex_seqs(sequences_to_pull_order);

fid = fopen("final_SATIVA_cleaned_files/final_SATIVA_cleaned_staph.fasta", 'wt');
fprintf(fid, '%s\n', final_fasta_file);
fclose(fid);

%% Cleaning Cutibacterium Species

clear all; close all;
fid=fopen('SILVA_v132_cutibacterium_CLEAN.fasta','r');
fasta = textscan(fid,'%s','delimiter','\n');
fasta = fasta{:};
fasta_clean = fasta;

fid2=fopen('SILVA_v132_cutibacterium_CLEAN.tax','r');
tax = textscan(fid2,'%s','delimiter','\n');
tax = tax{:};
tax_clean = tax;

fid3=fopen('sativa_outputs/SILVA_v132_cutibacterium_CLEAN.mis.tsv','r');
SATIVA_res = textscan(fid3,'%s','delimiter','\n');
SATIVA_res = SATIVA_res{:};
SATIVA_res_split = split(SATIVA_res(6:end,:),'	');

deletable_indexes = (str2double(SATIVA_res_split(:,5))<=.9);
label_swap_indexes = (str2double(SATIVA_res_split(:,5))>.9);

% Swap labels
for ii=1:length(label_swap_indexes)
    old_label = SATIVA_res_split((ii),6);
    new_label = SATIVA_res_split((ii),7);
    taxline = tax(contains(tax,SATIVA_res_split((ii),1)));
    newtaxline = replace(taxline, old_label, new_label);
    tax_clean(contains(tax,SATIVA_res_split((ii),1))) = newtaxline;
end

% Delete nessesary samps 
tax_clean(contains(tax, SATIVA_res_split(deletable_indexes,1))) = [];
delete_fasta_indexes = find(contains(fasta, SATIVA_res_split(deletable_indexes,1)));
delete_fasta_indexes = [delete_fasta_indexes; delete_fasta_indexes+1];
fasta_clean(delete_fasta_indexes) = [];

% Deals with another propinibacterium genus labelling issue 
tax_clean = cellfun(@(x) replace(x, "Propionibacterium propionicum" ,"Pseudopropionibacterium propionicum"), tax_clean, 'UniformOutput', false);

fid = fopen('final_SATIVA_cleaned_files/final_SATIVA_cleaned_cuti.fasta', 'wt');
fprintf(fid, '%s\n', string(fasta_clean));
fclose(fid);

fid = fopen('final_SATIVA_cleaned_files/final_SATIVA_cleaned_cuti.tax', 'wt');
fprintf(fid, '%s\n', string(tax_clean));
fclose(fid);

%% Cleaning all corynebacterium species
clear all; close all;
prefix = "SILVA_v132_corynebacterium_lawsonella_CLEAN";
fid_log=fopen("sativa_outputs/" + prefix + ".log",'r');
log_file = textscan(fid_log,'%s','delimiter','\n');
log_file = log_file{:};

taxa_to_replace = {};
taxa_replacement = strings(0);
all_taxa_combined_clusters = find(contains(log_file, ["WARNING: Following taxa share >60% indentical sequences und thus considered indistinguishable:","indistinguishable taxa"]));
for ii=1:length(all_taxa_combined_clusters)-1
    block_to_search = log_file(all_taxa_combined_clusters(ii):all_taxa_combined_clusters(ii+1));
    taxa_to_replace = [taxa_to_replace, {string(block_to_search(contains(block_to_search, [";Corynebacterium ", ";Lawsonella "])))}];
    taxa_replacement = [taxa_replacement, block_to_search(contains(block_to_search, ["TAXCLUSTER"]))];
end

fid=fopen(prefix + ".fasta",'r');
fasta = textscan(fid,'%s','delimiter','\n');
fasta = fasta{:};
fasta_clean = fasta;

fid2=fopen(prefix + ".tax",'r');
tax = textscan(fid2,'%s','delimiter','\n');
tax = tax{:};
tax_clean = tax;

fid3=fopen('sativa_outputs/' + prefix + ".mis.tsv",'r');
SATIVA_res = textscan(fid3,'%s','delimiter','\n');
SATIVA_res = SATIVA_res{:};
SATIVA_res_split = split(SATIVA_res(6:end,:),'	');

deletable_indexes = (str2double(SATIVA_res_split(:,5))<=.9);
label_swap_indexes = (str2double(SATIVA_res_split(:,5))>.9);

% Swap labels
for ii=1:length(label_swap_indexes)
    old_label = SATIVA_res_split((ii),6);
    new_label = SATIVA_res_split((ii),7);
    taxline = tax(contains(tax,SATIVA_res_split((ii),1)));
    newtaxline = replace(taxline, old_label, new_label);
    tax_clean(contains(tax,SATIVA_res_split((ii),1))) = newtaxline;
end

% Delete nessesary samps 
tax_clean(contains(tax, SATIVA_res_split(deletable_indexes,1))) = [];
delete_fasta_indexes = find(contains(fasta, SATIVA_res_split(deletable_indexes,1)));
delete_fasta_indexes = [delete_fasta_indexes; delete_fasta_indexes+1];
fasta_clean(delete_fasta_indexes) = [];

% Deals with another propinibacterium genus labelling issue 
for ii = 1:length(taxa_replacement)
    tax_clean = cellfun(@(x) replace(x, taxa_to_replace{ii} ,taxa_replacement(ii)), tax_clean, 'UniformOutput', false);
end

fid = fopen('final_SATIVA_cleaned_files/final_SATIVA_cleaned_coryne.fasta', 'wt');
fprintf(fid, '%s\n', string(fasta_clean));
fclose(fid);

fid = fopen('final_SATIVA_cleaned_files/final_SATIVA_cleaned_coryne.tax', 'wt');
fprintf(fid, '%s\n', string(tax_clean));
fclose(fid);

%% Cleaning nessereia species
clear all; close all;
prefix = "SILVA_v132_neisseria_CLEAN";
fid_log=fopen("sativa_outputs/" + prefix + ".log",'r');
log_file = textscan(fid_log,'%s','delimiter','\n');
log_file = log_file{:};

Neisseria_genuses = (";" + ["Alysiella", "Amantichitinum", "Aquaphilus", "Aquella", "Bergeriella", "Chitinolyticbacter", "Conchiformibius", "Crenobacter", "Eikenella", "Kingella", "Morococcus", "Neisseria", "Prolinoborus", "Rivicola", "Simonsiella", "Snodgrassella", "Stenoxybacter", "Uruburuella", "Vitreoscilla"] + " ");
taxa_to_replace = {};
taxa_replacement = strings(0);
all_taxa_combined_clusters = find(contains(log_file, ["WARNING: Following taxa share >60% indentical sequences und thus considered indistinguishable:","indistinguishable taxa"]));
for ii=1:length(all_taxa_combined_clusters)-1
    block_to_search = log_file(all_taxa_combined_clusters(ii):all_taxa_combined_clusters(ii+1));
    taxa_to_replace = [taxa_to_replace, {string(block_to_search(contains(block_to_search, Neisseria_genuses)))}];
    taxa_replacement = [taxa_replacement, block_to_search(contains(block_to_search, ["TAXCLUSTER"]))];
end

fid=fopen(prefix + ".fasta",'r');
fasta = textscan(fid,'%s','delimiter','\n');
fasta = fasta{:};
fasta_clean = fasta;

fid2=fopen(prefix + ".tax",'r');
tax = textscan(fid2,'%s','delimiter','\n');
tax = tax{:};
tax_clean = tax;

fid3=fopen('sativa_outputs/' + prefix + ".mis.tsv",'r');
SATIVA_res = textscan(fid3,'%s','delimiter','\n');
SATIVA_res = SATIVA_res{:};
SATIVA_res_split = split(SATIVA_res(6:end,:),'	');

deletable_indexes = (str2double(SATIVA_res_split(:,5))<=.9);
label_swap_indexes = (str2double(SATIVA_res_split(:,5))>.9);

% Swap labels
for ii=1:length(label_swap_indexes)
    old_label = SATIVA_res_split((ii),6);
    new_label = SATIVA_res_split((ii),7);
    taxline = tax(contains(tax,SATIVA_res_split((ii),1)));
    newtaxline = replace(taxline, old_label, new_label);
    tax_clean(contains(tax,SATIVA_res_split((ii),1))) = newtaxline;
end

% Delete nessesary samps 
tax_clean(contains(tax, SATIVA_res_split(deletable_indexes,1))) = [];
delete_fasta_indexes = find(contains(fasta, SATIVA_res_split(deletable_indexes,1)));
delete_fasta_indexes = [delete_fasta_indexes; delete_fasta_indexes+1];
fasta_clean(delete_fasta_indexes) = [];

% Deals with another propinibacterium genus labelling issue 
for ii = 1:length(taxa_replacement)
    tax_clean = cellfun(@(x) replace(x, taxa_to_replace{ii} ,taxa_replacement(ii)), tax_clean, 'UniformOutput', false);
end

fid = fopen('final_SATIVA_cleaned_files/final_SATIVA_cleaned_neisseria.fasta', 'wt');
fprintf(fid, '%s\n', string(fasta_clean));
fclose(fid);

fid = fopen('final_SATIVA_cleaned_files/final_SATIVA_cleaned_neisseria.tax', 'wt');
fprintf(fid, '%s\n', string(tax_clean));
fclose(fid);


%% Truncate the rest of the unmodified genes
clear all; close all; 
fid=fopen("SILVA_v132_not_specified_genus_unmodified.fasta",'r');
unmodfile = textscan(fid,'%s','delimiter','\n');
unmodfile = unmodfile{:};

taxaline = contains(unmodfile, '>');

split_labels = cellfun(@(x) split(x, ' '), unmodfile, 'UniformOutput', false);
new_fasta_file = cellfun(@(x) x(1), split_labels, 'UniformOutput', false);
%new_fasta_file(~taxaline) = cellfun(@(x) x(1:180), new_fasta_file(~taxaline), 'UniformOutput', false);

new_taxa_file = cellfun(@(x) erase(x, '>'), unmodfile(taxaline), 'UniformOutput', false);
new_taxa_file = cellfun(@(x) replace(x,' Bacteria;', '	Bacteria;'), new_taxa_file, 'UniformOutput', false);
new_taxa_file = cellfun(@(x) replace(x,' Eukaryota;', '	Eukaryota;'), new_taxa_file, 'UniformOutput', false);
new_taxa_file = cellfun(@(x) replace(x,' Archaea;', '	Archaea;'), new_taxa_file, 'UniformOutput', false);

split_by_tax_class = cellfun(@(x) split(x, ';'),new_taxa_file,'UniformOutput',false);
number_of_tax_class = cellfun(@length, split_by_tax_class);
badtax = number_of_tax_class~=7;
new_taxa_file(badtax) = [];
badtax_index_fasta = [find(badtax)*2-1; find(badtax)*2];
new_fasta_file(badtax_index_fasta) = [];

fid = fopen('final_SATIVA_cleaned_files/final_SATIVA_other_species.fasta', 'wt');
fprintf(fid, '%s\n', string(new_fasta_file));
fclose(fid);

fid = fopen('final_SATIVA_cleaned_files/final_SATIVA_other_species.tax', 'wt');
fprintf(fid, '%s\n', string(new_taxa_file));
fclose(fid);

%% 