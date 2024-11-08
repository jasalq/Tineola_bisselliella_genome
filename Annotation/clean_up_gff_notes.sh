# code to clean up GFF 

# install GFFUtils 
# download https://github.com/fls-bioinformatics-core/GFFUtils/releases
tar xzf GFFUtils-0.10.3.tar.gz

virtualenv venv
source venv/bin/activate
pip install -r ./GFFUtils-0.10.3/requirements.txt
pip install ./GFFUtils-0.10.3/

#when I tried to run gff_cleaner it gave me an error so I had to modify:
 nano /projectnb/mullenl/alqassar/software/GFFUtils/venv/lib/python3.10/site-packages/GFFUtils/GFFFile.py
 	# from collections import Iterator 
 	# to:
 	# from collections.abc import Iterator

cd /projectnb/mullenl/alqassar/wcm_annotation/3_emapper_functional_annotations/
gff_cleaner --clean emapper_output.emapper.decorated.gff  -o emapper_output.emapper.decorated.clean.gff  

#there was this error:
gffcleaner 0.12.0
Input : emapper_output.emapper.decorated.gff
Output: emapper_output.emapper.decorated.clean.gff
Replacing 'Anc_*' and blanks with '0's in 'score' column
Traceback (most recent call last):
  File "/projectnb/mullenl/alqassar/software/GFFUtils/venv/bin/gff_cleaner", line 8, in <module>
    sys.exit(main())
  File "/projectnb/mullenl/alqassar/software/GFFUtils/venv/lib/python3.10/site-packages/GFFUtils/cli/gff_cleaner.py", line 267, in main
    score_unexpected_values = sorted(list(score_unexpected_values))
TypeError: '<' not supported between instances of 'str' and 'float'

# trying to fix this way:
	pip install pandas 

	#make python script check_gff.py
	
# now retry with cleaned gff 
gff_cleaner --clean emapper_output.emapper.decorated.pre_clean.gff  -o emapper_output.emapper.decorated.clean.gff  

# that did not work so now trying to install funannotate
micromamba create -n funannotate -c bioconda -c conda-forge funannotate
eval "$(micromamba shell hook --shell bash)" # needed this because of this error: 'micromamba' is running as a subprocess and can't modify the parent shell.
												# Thus you must initialize your shell before using activate and deactivate.

micromamba activate funannotate
funannotate -h #to check install



awk 'NR==FNR{map[$1]=$2; next} {for (i=1; i<=NF; i++) if ($i in map) $i=map[$i]} 1' replacement.txt  emapper_output.emapper.decorated.clean.gff > output.gff



awk 'NR==FNR {map[$1] = $2; next} 
{
    for (key in map) {
        gsub(key, map[key]);
    }
    print
}' replacement.txt emapper_output.emapper.decorated.clean.gff > output.gff
# I renamed this gff to be JBFQDO000000000_annotation.gff


# to get the number of genes annotated
awk '$3 == "gene" {split($9, a, ";"); print a[1]}' wcm_rnd2_abinitio_wcm_pseudochromosomes.all.maker.all.gff | cut -d'=' -f2 > old_gene_ids.txt
# to get the number f genes functionally annotated
grep -Ff old_gene_ids.txt emapper_output.emapper.annotations.tsv > annotated_genes.tsv


# I need to filter out any annotations in the hard masked regions due to bacterial contamination
#generated a file with all the regions that are hardmasked 
#now will use bedtools to see regions that do not intersect with these hardmasked regions 
# but first need to rename to NCBI accession no in bed file 
awk 'NR==FNR {map[$1] = $2; next} {for (key in map) gsub(key, map[key]); print}' replacement.txt hardmasked_regions.bed > hardmasked_regions_renamed.bed




bedtools intersect -v -a JBFQDO000000000_annotation.gff -b hardmasked_regions_renamed.bed > JBFQDO000000000_filtered_annotations.gff


# to get the number of genes annotated
awk '$3 == "gene" {split($9, a, ";"); print a[1]}' JBFQDO000000000_filtered_annotations.gff | cut -d'=' -f2 > filtered_gene_ids.txt
wc -l filtered_gene_ids.txt

# need to rename gene names in emapper table
awk 'NR==FNR {map[$1] = $2; next} {for (key in map) gsub(key, map[key]); print}' replacement.txt emapper_output.emapper.annotations.tsv > emapper_output_renamed.emapper.annotations.tsv

# now can get number of genes functionally annotated in filtered gff
grep -Ff filtered_gene_ids.txt emapper_output_renamed.emapper.annotations.tsv > filtered_annotated_genes.tsv
wc -l filtered_annotated_genes.tsv 
# 11,037 functionally annotated out of 11,259 total genes 

# need to assess if any of the annotation statistics changed from functional annotation plus removing stuff in regions that were hard masked due to bacterial contamination
module load perl/5.28.1
  module load bioperl/1.7.2
  module load agat/0.7.0
  agat_sp_statistics.pl --gff JBFQDO000000000_filtered_annotations.gff -o JBFQDO000000000_filtered_annotations_summary

#to rename chromosomes in protein file
awk 'NR==FNR {map[$1] = $2; next} {for (key in map) gsub(key, map[key]); print}' replacement.txt wcm_rnd2_abinitio.all.maker.proteins.fasta > wcm_rnd2_abinitio.all.maker.proteins.renamed.fasta


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





