FungiSol: fungal_diversity_in_soil
__________________________________

Dans le test 2024_Aug, nous avons essayé deux marqueurs nucléaires, un long (SSU : 18S, autour de 1800bp) et un court (les deux régions ITS, autour de 600bp). Les blancs d’extraction ont été amplifiés et séquencés, ils étaient parfaits (avec les petits contaminants habituels), donc aucune contamination notable des échantillons durant la manip labo. Ils sont dans les fichiers bruts mais enlevés des fichiers de résultats.

Les séquences uniques (ASV) ont été regroupées en OTU à 97% de similitude. On peut aussi effectuer un regroupement à 99% pour ITS (non présenté ici, mais facilement calculable), voire à 100% si on travaille sur les lignées.

Les taxons 18S ont été assignés par rapport à la base de références SILVA 138 SSU NR99. Les taxons ITS à la base UNITE 04.04.2024. L’assignation d'un taxon est donnée si elle représente au moins 50% de toutes les assignations connues de la séquence dans la base de références (min bootstrap=0.5), on peut aussi régler ce paramètre à 80% ou sur une autre valeur. En dessous de la valeur de bootstrap, l’assignation est notée NA pour le rang taxonomique en question.

Une fois les singletons (1 seul read par échantillon) et les non-assignés enlevés (NA au rang Kingdom), le marqueur 18S détecte 593 OTUs sur les 3 échantillons, dont 151 OTUs fongiques et assimilés fongiques (Ascomycota, Basidiomycota, Blastocladiomycota, Chytridiomycota, Cryptomycota, Hyphochytriomycetes, Labyrinthulomycetes, Mucoromycota, Myxogastria, Peronosporomycetes). Le reste étant des protistes et invertébrés du sol.

Le marqueur ITS détecte 491 OTUs fongiques, mais uniquement dans les champignons « vrais » (classés phylogénétiquement dans le règne des FUNGI (Ascomycota, Basidiomycota, Chytridiomycota, Fungi_Incertae_sedis, Glomeromycota, Mortierellomycota (=anciens Zygomycota), Mucoromycota, Neocallimastigomycota, Rozellomycota et Sanchytriomycota). Il y avait aussi 139 OTUs non fongiques, mais je ne les ai pas assignées.

Sur l'ITS, l'échantillon 1D est dominé par les Eurotiomycètes (Ascomycètes), le EP05 par les Dothideomycètes et les Sordariomycètes (Ascomycètes tous les deux), le EP68 par les Mortierellomycètes (Mucoromycètes ou Mortierellomycètes selon les auteurs) et les Sordariomycètes (Ascomycètes).


<img width="747" alt="tests_juin2024_diversite_ITS" src="https://github.com/user-attachments/assets/87145385-b9c2-43df-b863-d2867f677bdb)">


%%%

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

Any inquiries? Send an email at tony.robinet@mnhn.fr

