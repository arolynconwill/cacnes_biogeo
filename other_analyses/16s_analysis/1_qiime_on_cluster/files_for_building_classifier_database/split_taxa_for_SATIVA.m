%% SPLIT TAXA 
% Cleaning steps: 
% 1. check to insure the correct number of genus mentions -- too little
% means no species classification
% 2. remove species that are "uncultured" or "unidentified" or contain
% "sp." (a poor species classification)
% 3. Remove any phrase after "subsp","Type" - 16S generally can't
% distinguish on a strain level basis
% 4. Change propionibacterium/cutibacterium mixups 

toggle_cuti_clean = 1;
toggle_coryne = 1;
genuses_to_clean_file_demon = ["cutibacterium", "staphylococcus","corynebacterium_lawsonella","neisseria"];

Neisseria_genuses = lower(["Alysiella", "Amantichitinum", "Aquaphilus", "Aquella", "Bergeriella", "Chitinolyticbacter", "Conchiformibius", "Crenobacter", "Eikenella", "Kingella", "Morococcus", "Neisseria", "Prolinoborus", "Rivicola", "Simonsiella", "Snodgrassella", "Stenoxybacter", "Uruburuella", "Vitreoscilla"]);
all_searched_for_genuses = ["cutibacterium", 'propionibacterium', 'staphylococcus', 'corynebacterium','lawsonella',Neisseria_genuses];

for ij = 1:length(genuses_to_clean_file_demon)
    newfilename = "dividing_up_raw_reads/SILVA_v132_" + genuses_to_clean_file_demon(ij) + "_unmodified.fasta";

    fid=fopen(newfilename,'r');
    alternate_textscan = textscan(fid,'%s','delimiter','\n');
    alternate_textscan_decompressed = alternate_textscan{:};
    % Because of intricacies within Neisseriaceae containing the genus Neisseria
    alternate_textscan_decompressed = replace(alternate_textscan_decompressed, "Neisseriaceae", "Neis8seriaceae");
    
    number_lines = length(alternate_textscan_decompressed);

    fasta_lines = contains(alternate_textscan_decompressed,'>');
    fasta_lines_text = alternate_textscan_decompressed(fasta_lines);
    sequence_lines = alternate_textscan_decompressed(~fasta_lines);

    split_line = cellfun(@(x) strsplit(x,[" Bacteria;", " Eukaryota;"]), fasta_lines_text, 'UniformOutput', false);
    fasta_names = cellfun(@(x) x(1), split_line);

    taxanames = string(cellfun(@(x) "Bacteria;" + x{2:end}, split_line, 'UniformOutput', false));

    % Check for number of genus occurances 
    number_of_genus_occurances = cellfun(@(x) sum(contains(split(lower(x), [" ",";"]), all_searched_for_genuses)), taxanames);
    correct_genus_number = number_of_genus_occurances == 2; 
    % Remove forbidden terms
    contains_forbidden_term = contains(taxanames, ["uncultured", "unidentified","sp."]);
    % Removes below length sequences 
    removed_sequence_length = cellfun(@length, sequence_lines)<180;

    % Remove subspecies and strain types
    length_split_remaining_labels = cellfun(@(x) length(split(x, ' ')), taxanames)<2; % This catches the case where the species designation is "genus;genus" and does not contain any actual species info
    remaining_labels = taxanames(correct_genus_number & ~contains_forbidden_term & ~removed_sequence_length & ~length_split_remaining_labels);
    if toggle_coryne == 1
        % Deals with a genus naming issue in coryne
        remaining_labels = cellfun(@(x) replace(x,";Corynebacterium 1;", ";Corynebacterium;" ), remaining_labels, 'UniformOutput', false);
    end
    
    split_remaining_labels = cellfun(@(x) split(x, ' '), remaining_labels, 'UniformOutput', false);
    clean_remaining_labels = cellfun(@(x) [x{1} + " " + x{2}], split_remaining_labels);

    % Remove weird brackets 
    clean_remaining_labels2 = cellfun(@(x) erase(x, ["[","]"]), clean_remaining_labels, 'UniformOutput', false);

    % Cutibacterium specific cleaning
    if toggle_cuti_clean == 1
        clean_remaining_labels2 = replace(clean_remaining_labels2, ';Propionibacterium acnes', ';Cutibacterium acnes');
        clean_remaining_labels2 = replace(clean_remaining_labels2, ';Propionibacterium humerusii', ';Cutibacterium humerusii');
        clean_remaining_labels2 = replace(clean_remaining_labels2, ';Propionibacterium namnetense', ';Cutibacterium namnetense');
        clean_remaining_labels2 = replace(clean_remaining_labels2, ';Propionibacterium avidum', ';Cutibacterium avidum');
        clean_remaining_labels2 = replace(clean_remaining_labels2, ';Propionibacterium granulosum', ';Cutibacterium granulosum');
    end

    % Combine to create final data structures
    final_taxafile = fasta_names(correct_genus_number & ~contains_forbidden_term & ~removed_sequence_length & ~length_split_remaining_labels) + "	" + clean_remaining_labels2;
    final_taxafile = cellfun(@(x) string(erase(x, '>')), final_taxafile);   
    final_taxafile = replace(final_taxafile, "Neis8seriaceae", "Neisseriaceae");
    final_fasta_file = strings(length(final_taxafile).*2,1);
    final_fasta_file(1:2:end) = fasta_names(correct_genus_number & ~contains_forbidden_term & ~removed_sequence_length & ~length_split_remaining_labels);
    
    useable_squence_lines_temp = sequence_lines(correct_genus_number & ~contains_forbidden_term & ~removed_sequence_length & ~length_split_remaining_labels);
    useable_squence_lines = cellfun(@(x) x(1:180), useable_squence_lines_temp, 'UniformOutput', false);
    final_fasta_file(2:2:end) = useable_squence_lines;

    % Write to new file
    cleanedfilename = "SILVA_v132_" + genuses_to_clean_file_demon(ij) + "_CLEAN";

    fid = fopen(cleanedfilename + ".fasta", 'wt');
    fprintf(fid, '%s\n', final_fasta_file);
    fclose(fid);

    fid = fopen(cleanedfilename + ".tax", 'wt');
    fprintf(fid, '%s\n', final_taxafile);
    fclose(fid);
end

%% SPLIT FILE GET RID OF TAXA
%{
fid=fopen('initial_reads_SILVA132.fna','r');
alternate_textscan = textscan(fid,'%s','delimiter','\n');
alternate_textscan_decompressed = alternate_textscan{:};

header_lines = contains(alternate_textscan_decompressed,'>');
header_lines_text = alternate_textscan_decompressed(header_lines);

header_lines_split = cellfun(@(x) split(x, ' '), header_lines_text);
%}
%% OLD CODE 
%{
parfor ii=1:length(genuses_to_clean)
    writtenfilename = "SILVA_v132_" + genuses_to_clean_file_demon(ii) + "_unmodified.fasta";
    clean_filename = "SILVA_v132_" + genuses_to_clean_file_demon(ii) + "_CLEAN.fasta";
    tax_filename = "SILVA_v132_" + genuses_to_clean_file_demon(ii) + "_CLEAN.tax";
    
    fid = fopen(writtenfilename,'r');
    fasta_file_modified = strings(0);
    taxafile = strings(0);

    gross_mislabelling_cat = 0;
    gross_mislabels = strings(0);

    toggle_switch_skip_taxa = 0;

    while ~feof(fid)
        tline = fgetl(fid);
        if contains(tline, '>')
            toggle_switch_skip_taxa = 0;
            split_line = strsplit(tline,' B');
            fasta_name = split_line{1};
            fasta_name = fasta_name(2:end);
            split_space_label = split(split_line{2}, ' ');
            split_space_label_v2 = split(split_line{2}, [" ",";"]);
            number_Cut_Prop_mentions = sum(contains(lower(split_space_label_v2), ["cutibacterium", "propionibacterium","staphylococcus",'corynebacterium']));

            % Check for gross mislabelling
            is_uncultured_bacterium = sum(contains(split_line{2}, ["uncultured", "unidentified"]));
            if number_Cut_Prop_mentions<2 & ~is_uncultured_bacterium
                gross_mislabelling_cat = gross_mislabelling_cat + 1;
                gross_mislabels = [gross_mislabels, split_line{2}];
            end

            % Check if modifications are needed to the species name or if the
            % species needs to be removed 
            if length(split_space_label)>2 || sum(contains(split_space_label, ["sp.", "unidentified","uncultured"]))>0 || number_Cut_Prop_mentions<2
                if sum(contains(split_space_label, ["subsp","Type"]))>0 && number_Cut_Prop_mentions>=2
                    % Remove subspecies or type info, as it is redundant and
                    % screws up classification
                    taxaname = fasta_name + "	"+ 'B'+ split_space_label{1} + ' ' + split_space_label{2};
                elseif sum(contains(split_space_label, ["sp.", "uncultured", "unidentified"]))>0 || number_Cut_Prop_mentions<2
                    % Remove any entry that is unidentified or uncultured on a
                    % species level
                    % Note: this dataset is messy and has species names like
                    % "Aureobasidium melanogenum". Because of this the number
                    % of Cutibacterium mentions needs to be greater than 2 to
                    % check for these cases
                    toggle_switch_skip_taxa = 1;  
                else
                    taxaname = fasta_name + "	"+ 'B'+ split_space_label{1} + ' ' + split_space_label{2};
                end
            else
                taxaname = fasta_name + "	"+ 'B'+ split_line(2);
            end

            if ~toggle_switch_skip_taxa
                fasta_file_modified = [fasta_file_modified; string(split_line(1))];
                taxaname = erase(taxaname,["[","]"]);
                taxafile = [taxafile; taxaname];
            end
        else
            if ~toggle_switch_skip_taxa
                periods_removed = replace(string(tline),'.','?');
                fasta_file_modified = [fasta_file_modified; periods_removed];
            end
        end

    end
    fclose(fid);

    taxafile_modified = replace(taxafile, ';Propionibacterium acnes', ';Cutibacterium acnes');
    taxafile_modified = replace(taxafile_modified, ';Propionibacterium humerusii', ';Cutibacterium humerusii');
    taxafile_modified = replace(taxafile_modified, ';Propionibacterium namnetense', ';Cutibacterium namnetense');
    taxafile_modified = replace(taxafile_modified, ';Propionibacterium avidum', ';Cutibacterium avidum');
    taxafile_modified = replace(taxafile_modified, ';Propionibacterium granulosum', ';Cutibacterium granulosum');


    fid = fopen(clean_filename', 'wt');
    fprintf(fid, '%s\n', fasta_file_modified);
    fclose(fid);

    fid = fopen(tax_filename, 'wt');
    fprintf(fid, '%s\n', taxafile_modified);
    fclose(fid);
end
%}