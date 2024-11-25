
# to identify what the bacterial contamination on chromosome 7 NCBI found is 

  # need to excise the region 9124443..12716725 from chromosome 7
  cd /projectnb/mullenl/alqassar/genome_bacterial_contam
  module load samtools
  samtools faidx /projectnb/mullenl/alqassar/synteny/satsuma_chromosemble_output/wcm_pseudochromosomes_renamed.fasta pseudochr_OX411251.1_7:9124443-12716725 > extracted_region_of_chr7.fasta


## split fasta into 300,000 bp chunks to be able to blast on the NCBI online portal

awk -v chunk_size=300000 '
BEGIN { header = ""; seq = ""; chunk_count = 0; }
{
    if ($0 ~ /^>/) {
        # Process the previous sequence if it exists
        if (seq != "") {
            while (length(seq) > 0) {
                chunk_count++;
                print header "_" chunk_count > (filename "_" chunk_count ".fasta");
                print substr(seq, 1, chunk_size) > (filename "_" chunk_count ".fasta");
                seq = substr(seq, chunk_size + 1);
            }
        }
        # Start a new header
        header = $0;
        chunk_count = 0;
        filename = substr(header, 2);  # Use header name for file naming
    } else {
        seq = seq $0;  # Concatenate sequence lines
    }
}
END {
    # Process the last sequence in the file
    if (seq != "") {
        while (length(seq) > 0) {
            chunk_count++;
            print header "_" chunk_count > (filename "_" chunk_count ".fasta");
            print substr(seq, 1, chunk_size) > (filename "_" chunk_count ".fasta");
            seq = substr(seq, chunk_size + 1);
        }
    }
}
' extracted_region_of_chr7.fasta



# to rename the names of chromosomes to NCBI_accession numbers for OSF repo in the synteny dataframes 

for i in *satsuma_summary.chained.out; do
#create backup of file
cp "$i" "${i}.bak"
# Use cut to isolate the first three columns and the fourth column separately
cut -f1-3 "$i" > "${i}.first_three_columns.txt"
    cut -f4 "$i" > "${i}.fourth_column.txt"
    cut -f5-8 "$i" > "${i}.last_three_column.txt"

while IFS=$'\t' read -r search_term replacement; do
  # Use sed to replace full word matches in input.txt
  sed -i "s/\b$search_term\b/$replacement/g" "${i}.fourth_column.txt"
done < replacement.txt

# Rejoin the columns: first three columns + modified fourth column
    paste -d'\t' "${i}.first_three_columns.txt" "${i}.fourth_column.txt" "${i}.last_three_column.txt" > "${i}.modified"

    #remove this weird ^M 
    sed -i 's/\r//g' "${i}.modified"

    # Rename the modified file
    mv "${i}.modified" "${i}_acc_renamed"

    # Clean up temporary files
    rm -f "${i}.first_three_columns.txt" "${i}.fourth_column.txt" "${i}.last_three_column.txt"
  
  echo "Replacement complete for $i. Output written to ${i}_acc_renamed."
done


# need to rename gene names in maker transcripts file to match new chromosome accession no
awk 'NR==FNR {map[$1] = $2; next} {for (key in map) gsub(key, map[key]); print}' replacement.txt wcm_rnd2_abinitio.all.maker.transcripts.fasta > wcm_rnd2_abinitio.all.maker.transcripts.renamed.fasta
