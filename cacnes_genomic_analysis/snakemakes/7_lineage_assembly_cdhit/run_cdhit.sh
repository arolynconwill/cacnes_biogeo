mkdir -p output
cd-hit -i input/lineages_all.fasta -o output/lineages_all_db_95 -c 0.95 -n 5 -M 64000 -d 0