#!/bin/bash


if [ "$#" -ne 3 ]; then
  echo "Usage: ./find_homologs.sh <query file> <subject file> <output file>"
  exit 1
fi


query_file=$1
subject_file=$2
output_file=$3


tblastn -query "$query_file" -subject "$subject_file" -out blast_results.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"


awk -v qlen="$qlen" '{if ($3 > 30 && $4 > (0.9 * qlen)) print $0}' blast_results.txt > "$output_file"

num_matches=$(wc -l < "$output_file")
echo "Number of matches found: $num_matches"


rm blast_results.txt
