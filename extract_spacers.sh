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

awk -v len="$query_length" '$3 == 100 && $4 == len && $5 == 0 && $6 == 0' temp_output.txt > perfect_matches.bed

awk '{
    if (NR > 1) {
        spacer_start = prev_end_pos + 1;
        spacer_end = $9 - 1;
        if (spacer_end > spacer_start) {
            print $2 "\t" spacer_start "\t" spacer_end;
        }
    }
    prev_end_pos = $10;
}' perfect_matches.bed > spacer_regions.bed


seqtk subseq "$subject_file" spacer_regions.bed > "$output_file"

num_spacers=$(wc -l < spacer_regions.bed)
echo "Number of spacers extracted: $num_spacers"


rm temp_output.txt perfect_matches.bed spacer_regions.bed


