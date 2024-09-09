#!/bin/bash


if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <protein query file> <nucleotide database> <output file>"
    exit 1
fi


query_file=$1
subject_db=$2
output_file=$3


query_length=$(grep -v ">" "$query_file" | tr -d '\n' | wc -c)


tblastn -query "$query_file" -db "$subject_db" -out temp_output.txt \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'


awk -v qlen="$query_length" '$3 > 30 && $4 > 0.9*qlen' temp_output.txt > "$output_file"


matches_count=$(wc -l < "$output_file")
echo "Number of matches: $matches_count"


rm temp_output.txt


