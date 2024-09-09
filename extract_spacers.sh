#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <query file> <subject file> <output file>"
    exit 1
fi

query_file=$1  # CRISPR repeat sequence (query)
subject_file=$2  # Genome assembly (subject)
output_file=$3  # Output file for extracted spacers

# Get the length of the query sequence
query_length=$(grep -v ">" "$query_file" | tr -d '\n' | wc -c)

# Perform BLAST to find perfect matches (CRISPR repeats)
blastn -query "$query_file" -subject "$subject_file" -task blastn-short \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' > temp_output.txt

# Identify CRISPR repeat locations and prepare the BED file for spacers
awk -v len="$query_length" 'BEGIN { prev_end = 0 }
{
    # Identify perfect matches and print spacer coordinates to BED file
    if ($3 == 100 && $4 == len && $5 == 0 && $6 == 0) {
        if (prev_end > 0) {
            spacer_start = prev_end + 1;
            spacer_end = $7 - 1;
            if (spacer_start < spacer_end) {
                print $2 "\t" spacer_start "\t" spacer_end;
            }
        }
        prev_end = $8;
    }
}' temp_output.txt > spacers.bed

# Use seqtk to extract the spacer sequences from the genome assembly
seqtk subseq "$subject_file" spacers.bed > "$output_file"

# Count and print the number of spacers extracted
spacer_count=$(grep -c ">" "$output_file")
echo "Number of spacers extracted: $spacer_count"

# Clean up temporary files
rm temp_output.txt spacers.bed

