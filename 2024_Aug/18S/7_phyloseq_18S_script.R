library(phyloseq)
library(microbiome)
library(metagMisc)
library(dplyr)
library(stringr)


#################### Creation of phyloseq objects
#################################################

samples <- read.table("~/sync/github/fungal_diversity_in_soil/2024_Aug/18S/8_phyloseq_files/samples_18S.csv", sep=",", header= T) %>% as.matrix()
rownames(samples) <- samples[,1]
samples %>% as_tibble()
sampledata = sample_data(data.frame(samples, row.names=rownames(samples), stringsAsFactors=FALSE))
sample_names(sampledata)


##########
###### 18S
abond18S <- read.csv2("~/sync/github/fungal_diversity_in_soil/2024_Aug/18S/8_phyloseq_files/OTU-table-97_test202408_18S_PacBio.csv", sep="\t", header= T) %>% as.data.frame()
rownames(abond18S) <- paste("e",(1:nrow(abond18S)))
abond18S <- mutate_all(abond18S, function(x) as.numeric(as.character(x)))
abond18S=as.matrix(abond18S)
colnames(abond18S) <- c("1D", "EP05", "EP68", "BL")
OTU_18S = otu_table(abond18S, taxa_are_rows = TRUE)

tax18S <- read.csv2("~/sync/github/fungal_diversity_in_soil/2024_Aug/18S/8_phyloseq_files/tax_97_18S_PacBio.csv", sep=";", header= F) %>% as.data.frame()
dim(tax18S)
rownames(tax18S) <- rownames(abond18S)
tax18S %>% head
colnames(tax18S) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
tax18S %>% head
tax18S <- as.matrix(tax18S)
TAX_18S = tax_table(as.matrix(tax18S))

test202408_18S = phyloseq(OTU_18S, TAX_18S, sampledata)

test202408_18S_2reads <- core(test202408_18S, detection = 2, prevalence = 0) # remove singletons

test202408_18S_2reads <- subset_samples(test202408_18S_2reads, sample!="BL")

test202408_18S_2reads %>% sample_sums()



test202408_18S_2reads@tax_table[,1] <- gsub(" d__", "", test202408_18S_2reads@tax_table[,1], fixed = TRUE)
test202408_18S_2reads@tax_table[,1] <- gsub("d__", "", test202408_18S_2reads@tax_table[,1], fixed = TRUE)
test202408_18S_2reads@tax_table[,1] <- gsub(" k__", "", test202408_18S_2reads@tax_table[,1], fixed = TRUE)
test202408_18S_2reads@tax_table[,1] <- gsub("k__", "", test202408_18S_2reads@tax_table[,1], fixed = TRUE)
test202408_18S_2reads@tax_table[,2] <- gsub(" p__", "", test202408_18S_2reads@tax_table[,2], fixed = TRUE)
test202408_18S_2reads@tax_table[,2] <- gsub("p__", "", test202408_18S_2reads@tax_table[,2], fixed = TRUE)
test202408_18S_2reads@tax_table[,3] <- gsub(" c__", "", test202408_18S_2reads@tax_table[,3], fixed = TRUE)
test202408_18S_2reads@tax_table[,3] <- gsub("c__", "", test202408_18S_2reads@tax_table[,3], fixed = TRUE)
test202408_18S_2reads@tax_table[,4] <- gsub(" o__", "", test202408_18S_2reads@tax_table[,4], fixed = TRUE)
test202408_18S_2reads@tax_table[,4] <- gsub("o__", "", test202408_18S_2reads@tax_table[,4], fixed = TRUE)
test202408_18S_2reads@tax_table[,5] <- gsub(" f__", "", test202408_18S_2reads@tax_table[,5], fixed = TRUE)
test202408_18S_2reads@tax_table[,5] <- gsub("f__", "", test202408_18S_2reads@tax_table[,5], fixed = TRUE)
test202408_18S_2reads@tax_table[,6] <- gsub(" g__", "", test202408_18S_2reads@tax_table[,6], fixed = TRUE)
test202408_18S_2reads@tax_table[,6] <- gsub("g__", "", test202408_18S_2reads@tax_table[,6], fixed = TRUE)


test202408_18S_2reads <- test202408_18S_2reads %>% subset_taxa(Kingdom!="Unassigned")
test202408_18S_2reads@tax_table[,2] %>% as.factor() %>% levels()
test202408_18S_2reads %>% subset_taxa(Phylum=="Ascomycota" | Phylum=="Basidiomycota" | Phylum=="Blastocladiomycota" | Phylum=="Chytridiomycota" |
                                        Phylum=="Cryptomycota" | Phylum=="Hyphochytriomycetes" | Phylum=="Hyphochytriomycetes" |
                                        Phylum=="Labyrinthulomycetes" | Phylum=="Mucoromycota" | Phylum=="Myxogastria" | Phylum=="Peronosporomycetes")

test202408_18S_2reads@tax_table[,2] %>% as.factor() %>% levels()

write.csv(test202408_18S_2reads@otu_table, "~/sync/github/fungal_diversity_in_soil/2024_Aug/18S/test202408_18S_2reads.csv")
write.csv(test202408_18S_2reads@tax_table, "~/sync/github/fungal_diversity_in_soil/2024_Aug/18S/taxa_test202408_18S_2reads.csv")



## Compute relative abundances
test202408_18S_2reads_comp <- microbiome::transform(test202408_18S_2reads, "compositional")

subset_species(test202408_18S_2reads)

####################### Fungal composition
library(ggplot2)
library(RColorBrewer)
#install.packages("upstartr")
library(upstartr)
library(ggpubr)
theme_set(theme_classic())

test202408_18S_2reads_comp_phy <- tax_glom(test202408_18S_2reads_comp, taxrank = "Phylum") # 43 phyla
test202408_18S_2reads_comp_cla <- tax_glom(test202408_18S_2reads_comp, taxrank = "Class") # 71 classes
test202408_18S_2reads_comp_fam <- tax_glom(test202408_18S_2reads_comp, taxrank = "Family") # 122 familles

test202408_18S_2reads_comp_fam_without_rares <- aggregate_rare(test202408_18S_2reads_comp_fam, level = "Family", detection = 0.05, prevalence = 0.01) # 11 familles + others
test202408_18S_2reads_comp_cla_without_rares <- aggregate_rare(test202408_18S_2reads_comp_cla, level = "Class", detection = 0.05, prevalence = 0.01) # 13 classes + others
test202408_18S_2reads_comp_phy_without_rares <- aggregate_rare(test202408_18S_2reads_comp_phy, level = "Phylum", detection = 0.05, prevalence = 0.01) # 9 phyla + others

test202408_18S_2reads_comp_cla_without_rares@tax_table


colourCount = test202408_18S_2reads_comp_cla_without_rares@tax_table[,2] %>% as.factor() %>% levels %>% length()
getPalette = colorRampPalette(brewer.pal(colourCount[1], "Paired"))



plot_composition(test202408_18S_2reads_comp_cla_without_rares) +
  scale_fill_manual(values = getPalette(colourCount[1])) + 
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_percent() +
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "18S complet Eukaryotes") + 
  #  theme_ipsum(grid="Y") +
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle=90, hjust=1, size=9),
        legend.text = element_text(face = "italic"), legend.position="bottom")


