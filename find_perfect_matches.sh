#!/bin/bash


if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <query file> <subject file> <output file>"
    exit 1
fi

query_file=$1
subject_file=$2
output_file=$3

query_length=$(grep -v ">" "$query_file" | tr -d '\n' | wc -c)

blastn -query "$query_file" -subject "$subject_file" -task blastn-short \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' > temp_output.txt

awk -v len="$query_length" '$3 == 100 && $4 == len && $5 == 0 && $6 == 0' temp_output.txt > "$output_file"

perfect_matches_count=$(wc -l < "$output_file")
echo "Number of perfect matches: $perfect_matches_count"


rm temp_output.txt

