library(phyloseq)
library(microbiome)
library(metagMisc)
library(dplyr)
library(stringr)


#################### Creation of phyloseq objects
#################################################

samples <- read.table("~/sync/github/fungal_diversity_in_soil/2024_Aug/ITS/8_phyloseq_files/samples_ITS.csv", sep=",", header= T) %>% as.matrix()
rownames(samples) <- samples[,1]
samples %>% as_tibble()
sampledata = sample_data(data.frame(samples, row.names=rownames(samples), stringsAsFactors=FALSE))
sample_names(sampledata)


##########
###### ITS
abondITS <- read.csv2("~/sync/github/fungal_diversity_in_soil/2024_Aug/ITS/8_phyloseq_files/OTU_table-97_test202408_ITS_PacBio.csv", sep="\t", header= T) %>% as.data.frame()
rownames(abondITS) <- paste("i",(1:nrow(abondITS)))
abondITS <- mutate_all(abondITS, function(x) as.numeric(as.character(x)))
abondITS=as.matrix(abondITS)
colnames(abondITS) <- c("1D", "EP05", "EP68", "BL")
OTU_ITS = otu_table(abondITS, taxa_are_rows = TRUE)

taxITS <- read.csv2("~/sync/github/fungal_diversity_in_soil/2024_Aug/ITS/8_phyloseq_files/tax_97_ITS_PacBio.csv", sep=";", header= F) %>% as.data.frame()
dim(taxITS)
rownames(taxITS) <- rownames(abondITS)
taxITS <- taxITS[,-8]
colnames(taxITS) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxITS %>% head
taxITS <- as.matrix(taxITS)
TAX_ITS = tax_table(as.matrix(taxITS))

test202408_ITS = phyloseq(OTU_ITS, TAX_ITS, sampledata)

test202408_ITS_2reads <- core(test202408_ITS, detection = 2, prevalence = 0) # remove singletons

test202408_ITS_2reads <- subset_samples(test202408_ITS_2reads, sample!="BL")

test202408_ITS_2reads %>% sample_sums()

test202408_ITS_2reads@tax_table[,1] <- gsub(" d__", "", test202408_ITS_2reads@tax_table[,1], fixed = TRUE)
test202408_ITS_2reads@tax_table[,1] <- gsub("d__", "", test202408_ITS_2reads@tax_table[,1], fixed = TRUE)
test202408_ITS_2reads@tax_table[,1] <- gsub(" k__", "", test202408_ITS_2reads@tax_table[,1], fixed = TRUE)
test202408_ITS_2reads@tax_table[,1] <- gsub("k__", "", test202408_ITS_2reads@tax_table[,1], fixed = TRUE)
test202408_ITS_2reads@tax_table[,2] <- gsub(" p__", "", test202408_ITS_2reads@tax_table[,2], fixed = TRUE)
test202408_ITS_2reads@tax_table[,2] <- gsub("p__", "", test202408_ITS_2reads@tax_table[,2], fixed = TRUE)
test202408_ITS_2reads@tax_table[,3] <- gsub(" c__", "", test202408_ITS_2reads@tax_table[,3], fixed = TRUE)
test202408_ITS_2reads@tax_table[,3] <- gsub("c__", "", test202408_ITS_2reads@tax_table[,3], fixed = TRUE)
test202408_ITS_2reads@tax_table[,4] <- gsub(" o__", "", test202408_ITS_2reads@tax_table[,4], fixed = TRUE)
test202408_ITS_2reads@tax_table[,4] <- gsub("o__", "", test202408_ITS_2reads@tax_table[,4], fixed = TRUE)
test202408_ITS_2reads@tax_table[,5] <- gsub(" f__", "", test202408_ITS_2reads@tax_table[,5], fixed = TRUE)
test202408_ITS_2reads@tax_table[,5] <- gsub("f__", "", test202408_ITS_2reads@tax_table[,5], fixed = TRUE)
test202408_ITS_2reads@tax_table[,6] <- gsub(" g__", "", test202408_ITS_2reads@tax_table[,6], fixed = TRUE)
test202408_ITS_2reads@tax_table[,6] <- gsub("g__", "", test202408_ITS_2reads@tax_table[,6], fixed = TRUE)
test202408_ITS_2reads@tax_table[,7] <- gsub(" s__", "", test202408_ITS_2reads@tax_table[,7], fixed = TRUE)
test202408_ITS_2reads@tax_table[,7] <- gsub("s__", "", test202408_ITS_2reads@tax_table[,7], fixed = TRUE)


test202408_ITS_2reads <- test202408_ITS_2reads %>% subset_taxa(Kingdom!="Unassigned")

## Compute relative abundances
test202408_ITS_2reads_comp <- microbiome::transform(test202408_ITS_2reads, "compositional")


####################### Fungal composition
library(ggplot2)
library(RColorBrewer)
#install.packages("upstartr")
library(upstartr)
library(ggpubr)
theme_set(theme_classic())

test202408_ITS_2reads_comp_phy <- tax_glom(test202408_ITS_2reads_comp, taxrank = "Phylum") # 10 phyla
test202408_ITS_2reads_comp_cla <- tax_glom(test202408_ITS_2reads_comp, taxrank = "Class") # 29 classes
test202408_ITS_2reads_comp_fam <- tax_glom(test202408_ITS_2reads_comp, taxrank = "Family") # 142 familles

test202408_ITS_2reads_comp_fam_without_rares <- aggregate_rare(test202408_ITS_2reads_comp_fam, level = "Family", detection = 0.05, prevalence = 0.01) # 6 familles + others
test202408_ITS_2reads_comp_cla_without_rares <- aggregate_rare(test202408_ITS_2reads_comp_cla, level = "Class", detection = 0.05, prevalence = 0.01) # 5 classes + others
test202408_ITS_2reads_comp_phy_without_rares <- aggregate_rare(test202408_ITS_2reads_comp_phy, level = "Phylum", detection = 0.05, prevalence = 0.01) # 3 phyla + others

test202408_ITS_2reads_comp_cla_without_rares@tax_table


colourCount = test202408_ITS_2reads_comp_cla_without_rares@tax_table[,2] %>% as.factor() %>% levels %>% length()
getPalette = colorRampPalette(brewer.pal(colourCount[1], "Paired"))



plot_composition(test202408_ITS_2reads_comp_cla_without_rares) +
  scale_fill_manual(values = getPalette(colourCount[1])) + 
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_percent() +
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "ITS1-5.8S-ITS2 Fungi") + 
  #  theme_ipsum(grid="Y") +
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle=90, hjust=1, size=9),
        legend.text = element_text(face = "italic"), legend.position="bottom")


