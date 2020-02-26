####DADA2 16S Amplicon Processing####
#When running this script it's best to run it interactively in the R environment to be sure each step works as expected
#there are also points within the protocol where you need to examine some output files to make decisions about quality
#filtering, trimming, etc. Once you've run through the script and are sure everything is up to snuff you can re-run
#from the terminal using: Rscript dada2.for_16s_dat.R whenever you need to (e.g., you noticed a mistake and you need to reprocess)

#may need to start R in the following way if you encounter vector memory full error -- can also add to bash profile: env R_MAX_VSIZE=700Gb Rscript dada2_processing.3.6.1.R

####CHANGE THESE####
FWPRI="GTGYCAGCMGCCGCGGTAA" #CHANGE ME to your forward primer
RVPRI="GGACTACNVGGGTWTCTAAT" #CHANGE ME to your reverse primer
FPAT="_R1_001.fastq.gz" #CHANGE ME to match the pattern in your forward reads
RVPAT="_R2_001.fastq.gz" #CHANGE ME to match the pattern in your reverse reads
BACTAX="01_data/ezbiocloud_id_taxonomy.qza" #change this to your bacterial reference db taxonomy file
BACREF="01_data/ezbiocloud_qiime_full.qza" #change this to your bacterial reference db fasta file
#####################

####Required Libraries####
library(dada2)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(qualpalr)
library(viridis)
library(ShortRead)
library(Biostrings)
library(seqinr)
library(phyloseq)

####Environment Setup####
theme_set(theme_bw())

####File Path Setup####
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
list.files("01_data")
fnFs <- sort(list.files("01_data", pattern=FPAT, full.names = TRUE))
fnRs <- sort(list.files("01_data", pattern=RVPAT, full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #CHANGE the delimiter in quotes and the number at the end of this command to decide how to split up the file name, and which element to extract for a unique sample name

####fastq Quality Plots####
# be sure to take a look at these, they should inform your read trimming decisions downstream
pdf("03_output/forward_quality_plot.pdf")
plotQualityProfile(fnFs[1:2]) #this plots the quality profiles for each sample, if you have a lot of samples, it's best to look at just a few of them, the plots take a minute or two to generate even only showing 10-20 samples.
dev.off()
pdf("03_output/reverse_quality_plot.pdf")
plotQualityProfile(fnRs[1:2])
dev.off()

####Primer Removal####
####identify primers####
FWD <- FWPRI
REV <- RVPRI
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path("01_data", "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path("01_data", "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[2]]), #add the index of the sample you'd like to use for this test (your first sample may be a blank/control and not have many sequences in it, be mindful of this)
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[2]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[2]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[2]]))

#### primer removal ####
system2(cutadapt, args = "--version")

path.cut <- file.path("01_data", "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

#Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#sanity check, should report zero for all orientations and read sets
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2", full.names = TRUE))

####filter and trim reads####
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

####trim & filter####
#filter and trim command. dada2 can canonically handle lots of errors, I am typically permissive in the maxEE parameter set here, in order to retain the maximum number of reads possible. error correction steps built into the dada2 pipeline have no trouble handling data with this many expected errors.
#check out the DADA2 tutorial to see what each of the options in the command below do
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, trimLeft=5, trimRight=25, minLen = c(150,120),
                     maxN=c(0,0), maxEE=c(2,2), truncQ=c(2,2), rm.phix=TRUE, matchIDs=TRUE,
                     compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
#check this table, you should retain most of your reads. If you lose a lot you might need to change the filterAndTrim quality thresholds
#View(retained)

####learn error rates####
#the next three sections (learn error rates, dereplication, sample inference) are the core of dada2's sequence processing pipeline.
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

pdf("03_output/error_plot.pdf")
plotErrors(errF, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered before moving on.
dev.off()

####dereplication####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names #this is just to ensure that all your R objects have the same sample names in them
names(derepFs) <- sample.names
names(derepRs) <- sample.names

####sample inference####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]

####OPTIONAL: remove low-sequence samples before merging####
#a "subscript out of bounds" error at the next step (merging) may indicate that you aren't merging any reads in one or more samples.
#samples_to_keep <- as.numeric(out[,"reads.out"]) > 500 #example of simple method used above after the filter and trim step. if you already did this but still got an error when merging, try the steps below
getN <- function(x) sum(getUniques(x)) #keeping track of read retention, number of unique sequences after ASV inference
track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN))
samples_to_keep <- track[,4] > 500 #your threshold. try different ones to get the lowest one that will work. #this method accounts for dereplication/ASVs left after inference
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)] #record names of samples you have the option of removing

####merge paired reads####
mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], derepRs[samples_to_keep], verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

####construct sequence table####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

####View Sequence Length Distribution Post-Merging####
#most useful with merged data. this plot will not show you much for forward reads only, which should have a uniform length distribution.
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
pdf("03_output/length_hist.pdf")
plot(x=length.histogram[,1], y=length.histogram[,2]) #view length distribution plot
dev.off()

####remove low-count singleton ASVs####
#create phyloseq otu_table
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

#some metrics from the sequence table
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV
length(otu_pres_abs_rowsums) #how many ASVs
length(which(otu_pres_abs_rowsums == 1)) #how many ASVs only present in one sample

#what are the counts of each ASV
otu_rowsums <- rowSums(otus) #raw counts per ASV
otu_singleton_rowsums <- as.data.frame(otu_rowsums[which(otu_pres_abs_rowsums == 1)]) #raw read counts in ASVs only presesnt in one sample
hist(otu_singleton_rowsums[,1], breaks=500, xlim = c(0,200), xlab="# Reads in ASV") #histogram plot of above
length(which(otu_singleton_rowsums <= 50)) #how many ASVs are there with N reads or fewer? (N=50 in example)

#IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's "transform_sample_counts" to create a relative abundance table, and then filter your ASVs by choosing a threshold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
dim(seqtab) # sanity check
dim(otus) # (this should be the same as last command, but the dimensions reversed)
otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x)) #create relative abundance table
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame
df[is.na(df)] <- 0 #if there are samples with no merged reads in them, and they passed the merge step (a possiblity, converting to a relative abundance table produes all NaNs for that sample. these need to be set to zero so we can do the calculations in the next steps.)
otus_rel_ab.rowsums <- rowSums(df) #compute row sums (sum of relative abundances per ASV. for those only present in one sample, this is a value we can use to filter them for relative abundance on a per-sample basis)
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1) #which ASVs are only present in one sample
b <- which(otus_rel_ab.rowsums <= 0.001) #here is where you set your relative abundance threshold #which ASVs pass our filter for relative abundance
length(intersect(a,b)) #how many of our singleton ASVs fail on this filter
rows_to_remove <- intersect(a,b) #A also in B (we remove singleton ASVs that have a lower relative abundance value than our threshold)
otus_filt <- otus[-rows_to_remove,] #filter OTU table we created earlier
dim(otus_filt) #how many ASVs did you retain?
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence table matrix to continue with dada2 pipeline

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="pooled", multithread=TRUE, verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset
dim(seqtab.nosingletons.nochim)
sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons) #proportion of nonchimeras #it should be relatively high after filtering out your singletons/low-count ASVs, even if you lose a lot of ASVs, the number of reads lost should be quite low

####track read retention through steps####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN), rowSums(seqtab.nosingletons), rowSums(seqtab.nosingletons.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(track, 100-track[,6]/track[,5]*100, 100-track[,7]/track[,6]*100)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nosingletons", "nochimeras", "percent_singletons", "percent_chimeras")
rownames(track) <- sample.names[samples_to_keep]

####save output from sequnce table construction steps####
write.table(data.frame("row_names"=rownames(track),track),"read_retention.16s.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"sequence_table.16s.merged.txt", row.names=FALSE, quote=F, sep="\t")
uniquesToFasta(seqtab.nosingletons.nochim, "rep_set.fa")
system("awk '/^>/{print \">ASV\" ++i; next}{print}' < rep_set.fa > rep_set_fix.fa")

################Load back into R

#if you must save your sequence table and load it back in before doing taxonomy assignments, here is how to reformat the object so that dada2 will accept it again
seqtab.nosingletons.nochim <- fread("sequence_table.16s.merged.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab.nosingletons.nochim) <- seqtab.nosingletons.nochim[,1] #set row names
seqtab.nosingletons.nochim <- seqtab.nosingletons.nochim[,-1] #remove column with the row names in it
seqtab.nosingletons.nochim <- as.matrix(seqtab.nosingletons.nochim) #cast the object as a matrix

# now replace the long ASV names (the actual sequences) with human-readable names, and save the new names and sequences as a .fasta file in your project working directory
my_otu_table <- t(as.data.frame(seqtab.nosingletons.nochim)) #transposed (OTUs are rows) data frame. unclassing the otu_table() output avoids type/class errors later on
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') #create new names
colnames(seqtab.nosingletons.nochim) <- ASV.num #rename your ASVs in the taxonomy table and sequence table objects

##assign taxonomy using VSEARCH and QIIME2
#first need to convert reference sequences to qiime2 format
#reference database should be in qiime format: qiime tools import --input-path ezbiocloud_qiime_full.fasta --output-path ezbiocloud_qiime_full.qza --type 'FeatureData[Sequence]'
#also taxonomy file: qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path ezbiocloud_id_taxonomy.txt --output-path ezbiocloud_id_taxonomy.qza
#pull the region of interest from full length 16S sequences (is better for classification): qiime feature-classifier extract-reads --i-sequences ezbiocloud_qiime_full.qza --p-f-primer GTGYCAGCMGCCGCGGTAA --p-r-primer GGACTACNVGGGTWTCTAAT --p-min-length 200 --p-max-length 300 --o-reads ezbiocloud_qiime_v4region.qza
system("source activate qiime2019.10 ; qiime tools import --input-path rep_set_fix.fa --output-path rep_set_fix.qza --type 'FeatureData[Sequence]'")
system("rm -r assigntax")
system(sprintf("source activate qiime2019.10 ; qiime feature-classifier classify-consensus-vsearch --i-query rep_set_fix.qza --i-reference-reads %s --i-reference-taxonomy %s --output-dir assigntax", BACREF, BACTAX))
system("unzip assigntax/classification.qza -d assigntax/")

#get file path for taxonomy file
tempfile <- subset(dir(path="assigntax"), !grepl("classification.qza", dir(path="assigntax/")))
newpath <- paste("assigntax/", tempfile, "/data/taxonomy.tsv", sep="")
#make 7 level taxonomy file
system(sprintf("awk '{print $2}' %s | sed '1d' > taxonomy_strings.txt", newpath))
system(sprintf("awk '{print $1}' %s | sed '1d' > asv_ids.txt", newpath))
system("python 02_code/fix_taxonomy_L7.py taxonomy_strings.txt > fix_string.txt") #the script fix_taxonomy_L7.py should be in your working directory
system("paste asv_ids.txt fix_string.txt > taxonomy_L7.txt")
system("rm taxonomy_strings.txt fix_string.txt asv_ids.txt")

####combine sequence and taxonomy tables into one####
#taxa will be the rows, columns will be samples, followed by each rank of taxonomy assignment, from rank1 (domain-level) to rank7/8 (species-level), followed by accession (if applicable)
#first check if the row names of the taxonomy table match the column headers of the sequence table
taxa <- read.table(newpath, header=T, sep="\t", row.names=1)
length(which(row.names(taxa) %in% colnames(seqtab.nosingletons.nochim)))
dim(taxa)
dim(seqtab.nosingletons.nochim)
#the number of taxa from the last three commands should match

#now ensure that the taxa in the tables are in the same order #this should be true if you haven't reordered one or the other of these matrices inadvertently
order.col <- row.names(taxa)
seqtab.nosingletons.nochim <- seqtab.nosingletons.nochim[,order.col]
row.names(taxa) == colnames(seqtab.nosingletons.nochim) #IMPORTANT: only proceed if this evaluation is true for every element. if it isn't you'll need to re-order your data. I'd suggest sorting both matrices by their rows after transposing the sequence table.

#as long as the ordering of taxa is the same, you can combine like this (note you need to transpose the sequence table so that the taxa are in the rows)
sequence_taxonomy_table <- cbind(t(seqtab.nosingletons.nochim), taxa)
colnames(sequence_taxonomy_table)[colnames(sequence_taxonomy_table)=="V2"] <- "taxonomy"
colnames(sequence_taxonomy_table)[colnames(sequence_taxonomy_table)=="V3"] <- "score"
colnames(sequence_taxonomy_table)[colnames(sequence_taxonomy_table)=="V4"] <- "count"

#now write to file
write.table(data.frame("row_names"=rownames(sequence_taxonomy_table),sequence_taxonomy_table),"sequence_taxonomy_table.16s.merged.txt", row.names=FALSE, quote=F, sep="\t")

#filter out unwanted taxonomic groups
system("grep -v 'Unassigned' sequence_taxonomy_table.16s.merged.txt | awk '{print $1}' | grep 'A' > wanted.ids")
wanted <- read.table("wanted.ids", header=F)
seqtab.filtered <- seqtab.nosingletons.nochim[, which(colnames(seqtab.nosingletons.nochim) %in% wanted$V1)]

#get representative seq tree
system("seqtk subseq rep_set_fix.fa wanted.ids > rep_set.filt.fa")
system("mafft --auto rep_set.filt.fa > rep_set.filt.align.fa")
system("fasttree -nt rep_set.filt.align.fa > rep_set.filt.tre")