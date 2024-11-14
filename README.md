<head> <h1> Tineola_bisselliella_genome </h1> </head>

<body>  

This repository contains all code associated with the manuscript "De Novo Genome Assembly and Annotation for the Webbing Clothes Moth (<em>Tineola bisselliella</em>): A Globally Distributed, Economically Important Pest", which has been submitted to <em>Genome Biology and Evolution</em>. The goal of this project was to assemble and annotate a PacBio sequenced <em> Tineola biselliella </em> genome. It used a modified version of the MAKER 3.01.04 annotation pipeline (Cantarel et al., 2008). We also generated an annotated transcriptome containing RNA sequence data from larvae and adults and performed synteny alignments between the <em>Tineola bisselliella </em> pseudochromosome genome assembly we generated and other lepidopteran genomes.

Contact info: Jasmine Alqassar (j.alqassar@gwu.edu)
</body>
<h2> Genome Sequencing and Assembly </h2>
<body>Two <em> Tineola bisselliella </em> larvae were obtained from an infested skunk fur; their DNA was extracted and sent for PacBio HiFi Sequencing at the University of Deleware. 11 Gb of HiFI reads were assembled using Improved Phased Assembler (IPA) (GitHub: https://github.com/PacificBiosciences/pbipa). This contig level assembly was then assembled into scaffolds and 30 pseudochromosomes using synteny assisted mapping to a closely related species with a chromosome-level assembly,<em> Tinea pellionella</em> (GenBank accession no. GCA_948150575.1), by Satsuma Chromosembler v.3.0 (Grabherr et al., 2010). BUSCO was performed using the Lepidoptera-odb10 lineage dataset (5,286 orthologs) to assess genome assembly quality. BUSCO was also performed using the Bacteria-odb9 lineage dataset (148 orthologs) and NCBI performed a contamination screen during the submission process to identify bacterial contamination. Small regions of bacterial contamination were masked from the final assmebly using BEDtools. </body>
  
<h2> Transcriptome Sequencing, Assembly, and Annotation </h2>
<body> Ten live <em>T. bisselliella</em> individuals (five early-stage larvae, three late-stage larvae, and two adults) were sampled from a colony that was established in the lab in June 2022. RNA was extracted from each individual using a Qiagen RNeasy Micro Kit. Samples were separated into three barcoded pools based on developmental stage (i.e., early-stage larvae, late-stage larvae, or adult) and were then prepped using the Illumina Stranded mRNA Prep method. 100 bp paired-end sequencing of all three libraries was performed on an Illumina NovaSeq Flowcell. The paired-end RNAseq data were assembled and annotated following the protocol presented in Rivera and Davies (2021). </body>

<h2> Genome Annotation </h2>
<body> The pseudochromosome genome assembly was annotated following the original MAKER pipeline detailed in Cantarel et al. (2008) using the most up-to-date version of the program, MAKER v3.01.04, with the modifications detailed in the manuscript. <body> 

<h2> Synteny Analysis </h2>
<body>Satsuma was used to perform the synteny analysis (Grabherr et al., 2010), and the R package RIdeogram was used to visualize the output of Satsuma (Hao et al., 2020). Additionally, the synteny alignment between the <em>Tineola bisselliella </em> pseudochromosome genome assembly and <em>Tinea pellionella </em> genome was visualized as a Circos plot using the R package Circlize (Gu et al., 2014). </body>

<h2> Master's Thesis Only Analyses: Differential Expression and Gene Ontology (GO) Enrichment Analyses </h2>
<body> As a part of my master's thesis preliminary differential expression and Gene Ontology (GO) enrichment analyses using RNAseq data from early-stage larvae, late-stage larvae, and adults.Trimmed RNAseq paired-end reads generated previously were mapped to the assembled transcriptome using Bowtie2 (--no-unal --local --score-min L,16,1 -L 16 -I 0 -X 1500) (Langmead & Salzberg, 2012). Using the tag-based RNAseq protocol developed in Meyer et al. (2011) read counts per gene were generated. The R package DESeq2 was used to identify differentially expressed genes (DEGs) across life stages (adult vs. larvae) (Love et al., 2014). A heat map of 5,066 DEGs across life stages was produced using the heatmap.2 R function. Rank-based gene ontology enrichment analysis was performed with the differentially expressed genes using the adaptive clustering of GO categories and Mann-Whitney U tests based on ranking of signed log p-values (GO-MWU) method from (Wright et al., 2015).
</body>

