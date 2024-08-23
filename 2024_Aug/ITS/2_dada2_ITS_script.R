# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("dada2")
library(dada2)
library(dplyr)
library(Biostrings)

## https://benjjneb.github.io/LRASManuscript/LRASms_fecal.html
## https://next-its.github.io/

## In order to install dada2 package from BiocManager, you need an access token to github repository
## first sign in (or sign up) to your github account
## then follow this : https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens
## after what you can set up your github parameters in R
## usethis::use_git_config(user.name = "xxx", user.email = "yyy")
## GitHub personal token (valid until nonv 2024):  zzz
## Sys.setenv(GITHUB_PAT = "zzz") # to update the GitHub token
## Sys.getenv("GITHUB_PAT") # to check what is the current token



##################### ITS
#######################################################@

path1 <-"~/sync/github/fungal_diversity_in_soil/2024_Aug/ITS/1_fastq_files" # change with your path
ITS1F <- "CTTGGTCATTTAGAGGAAGTAA"
ITS4 <- "TCCTCCGCTTATTGATATGC"
rc <- dada2:::rc

list.files(path1)
setwd(path1)
# il ne doit y avoir dans le path que les fichiers fastq (pas les .fastq.gz)

fnFs <- sort(list.files(path1, pattern="_clean.fastq.gz", full.names = TRUE))
fnFs <- sort(fnFs)

nops <- c("~/sync/metabarcoding/fonges_pathogenes_Nantes/test_240808/ITS/1_fastq/NOPS/Fungus_CC151-002P0001_clean.fastq_NOPS.gz",
          "~/sync/metabarcoding/fonges_pathogenes_Nantes/test_240808/ITS/1_fastq/NOPS/FFungus_CC151-002P0002_clean.fastq_NOPS.gz",
          "~/sync/metabarcoding/fonges_pathogenes_Nantes/test_240808/ITS/1_fastq/NOPS/FFungus_CC151-002P0003_clean.fastq_NOPS.gz",
          "~/sync/metabarcoding/fonges_pathogenes_Nantes/test_240808/ITS/1_fastq/NOPS/FFungus_CC151-003P0001_clean.fastq_NOPS.gz")


prim2 <- removePrimers(fnFs, nops, primer.fwd=ITS1F, primer.rev=dada2:::rc(ITS4), orient=TRUE)

lens.fn <- lapply(nops, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)



# filter
path2 <-"~/sync/github/fungal_diversity_in_soil/2024_Aug/ITS/1_fastq_files/NOPS/" # change with your path
filts2 <- file.path(path2, "noprimers", "filtered", basename(fnFs))
track2 <- filterAndTrim(nops, filts2, minQ=3, minLen=300, maxLen=900, maxN=0, rm.phix=FALSE, maxEE=2)
track2

# dereplicate
derepFs <- derepFastq(filts2, verbose=TRUE)

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "CC151_"), `[`, 1)
sample.names
names(derepFs) <- sample.names


# learn errors
err2 <- learnErrors(derepFs, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)

# denoise
dd2 <- dada(derepFs, err=err2, BAND_SIZE=32, multithread=TRUE)
cbind(ccs=prim2[,1], primers=prim2[,2], filtered=track2[,2], denoised=sapply(dd2, function(x) sum(x$denoised)))

# sequences table
st2 <- makeSequenceTable(dd2); dim(st2)
# Check chimeras
bim2 <- isBimeraDenovo(st2, minFoldParentOverAbundance=3.5, multithread=TRUE)
table(bim2)
sum(st2[,bim2])/sum(st2)
# Remove chimeras
st2_nochim <- removeBimeraDenovo(st2, method = "consensus", verbose = FALSE)
st2_nochim %>% dim()
t(st2_nochim) %>% head

# then we use qiime2 for using vsearch 0.97 / 0.99 and BLAST tools
write.table(t(st2_nochim), "~/sync/github/fungal_diversity_in_soil/2024_Aug/ITS/3_dada2_files/st2_nochim_ITS.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

uniquesToFasta(st2_nochim, "~/sync/github/fungal_diversity_in_soil/2024_Aug/ITS/3_dada2_files/rep-seqs_ITS.fna", 
               ids=colnames(st2_nochim))

