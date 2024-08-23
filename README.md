fungal_diversity_in_soil
________________________

Workflow for genetic analysis of fungi in soil samples

%%%

DNA has been extracted from soil samples with an appropriate extraction kit, or soft extraction method. Here we used the PowerSoil kit from Qiagen.

Then DNA extracts were quantified with a fluorometry method (here Qbit) and their purity estimated with a micro-droplet spectrophotometer (here NanoVue).

Full 18S and ITS1-5.8S-ITS2 regions were amplified, indexed in libraries and sequenced on PacBio Revio (here by BMKGENE from BMK GmbH).

%%%

Fastq files in the directory "./1_fastq_files" are then analyzed, following this workflow:

    (dada2 in R) Run the script "./2_dada2_18S_script.R" and/or "./2_dada2_ITS_script.R". Trim the primers, draw read length distribution, filter reads on length and quality, dereplicate reads, learn errors, denoise, remove PCR chimeras. Read abundances by sample and the corresponding sequences are placed in "./3_dada2_files", respectively "st2_nochim_18S.txt" and "rep-seqs_18S.fna" for 18S, "st2_nochim_ITS.txt" and "rep-seqs_ITS.fna" for ITS.

    (qiime2 in Jupyter lab) By running the script "4_vsearch-blast_xxx_script.ipynb", make sequence consensus at 97% and 99% similarity, from "rep-seqs_xxx.fna" (replace xxx by 18S or ITS), in order to group ASVs (Amplicon Sequence Variants) into OTUs (Operational Taxonomic Units). Assign taxonomy to these sequences, based on SILVA and UNITE databases, respectively for 18S and ITS regions. OTU tables are placed in "./4_vsearch-analysis" as tsv files (feature-table-97_test202408_18S_PacBio.tsv and feature-table-97_test202408_ITS_PacBio.tsv), but taxonomic assignations are placed in "./5_blast-assignation/taxonomy.tsv".

    (table editor like LibreOffice Calc) Edit the files "./5_vsearch_files/feature-table-97_test202408_xxx_PacBio.tsv" and "./6_assignation_files/taxonomy.tsv" in order to produce OTU tables with taxonomy and read abundances ; save read abundances only in a csv file, with sample names at headers but no taxa names (OTU-table-97_test202408_18S_PacBio.csv and OTU-table-97_test202408_IT_PacBio.csv), and save taxonomy only in another csv file, with no headers but with one taxonomical rank per column, with strictly the same OTU ranking than read abundances (tax_97_18S_PacBio.csv and tax_97_ITS_PacBio.csv). A third csv file must give the sample names, with headers (samples_18S.csv and samples_ITS.csv). There are 3 new files per marker, placed in the directory "./8_phyloseq-files".

    (phyloseq in R) Run (and edit if needed) the script "./7_phyloseq_18S_script.R" or "./7_phyloseq_ITS_script.R" will make phyloseq objects with 18S and ITS outputs, delete singletons and Unassigned OTUs, possibility to make downstream analysis on composition, multivariate analysis etc.

%%%

Any inquiries, send me an email : tony.robinet@mnhn.fr

