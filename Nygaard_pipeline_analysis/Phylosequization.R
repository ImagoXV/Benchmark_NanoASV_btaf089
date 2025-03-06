#NanoASV phylosequisation
#Arthur Cousson - 2023
#Contact : arthur.cousson@ird.fr
#Workaround conflicting host vs local libraries

# .libPaths(c("/usr/local/lib/R/site-library","/usr/lib/R/site-library","/usr/lib/R/library"))
# 
# args <- commandArgs(trailingOnly = TRUE)
# DIR <- args[1]
# OUTPWD <- args[2]
# R_CLEANING <- args[3]
# TREE <- args[4]
# METADATA <- args[5]


# Nygaard pipeline output phyloseq generation

setwd("~/Documents/Thesis/NanoASV/Software_dev/Benchmark_Nygaard/Nygaard_pipeline_analysis/")

#Phylosequization -----
library(phyloseq)
library(ape) #To handle trees
metadata <- read.csv("metadata.csv", row.names = 1, header = TRUE, check.names = FALSE)
#metadata <- read.csv(paste0(METADATA,"/metadata.csv"), row.names = 1, header = TRUE, check.names = FALSE)
barcodes <- rownames(metadata)

##ASV tables ----
#Individuals ASV tables loading
temp_ASV = list.files(path = "./", pattern = "*.tsv")


temp_ASV <- list()
temp_ASV[[1]] <- read.csv("ASV_TABLE_filtered_choped_barcode01.fasta.blast.tab.tsv", check.names = F, sep = "\t", header = F, row.names = 1)
temp_ASV[[2]] <- read.csv("ASV_TABLE_filtered_choped_barcode02.fastq.gz.fasta.blast.tab.tsv", check.names = F, sep = "\t", header = F, row.names = 1)
temp_ASV[[3]] <- read.csv("ASV_TABLE_filtered_choped_barcode03.fastq.gz.fasta.blast.tab.tsv", check.names = F, sep = "\t", header = F, row.names = 1)
temp_ASV[[4]] <- read.csv("ASV_TABLE_filtered_choped_barcode04.fastq.gz.fasta.blast.tab.tsv", check.names = F, sep = "\t", header = F, row.names = 1)
temp_ASV[[5]] <- read.csv("ASV_TABLE_filtered_choped_barcode05.fastq.gz.fasta.blast.tab.tsv", check.names = F, sep = "\t", header = F, row.names = 1)
temp_ASV[[6]] <- read.csv("ASV_TABLE_filtered_choped_barcode06.fastq.gz.fasta.blast.tab.tsv", check.names = F, sep = "\t", header = F, row.names = 1)
temp_ASV[[7]] <- read.csv("ASV_TABLE_filtered_choped_barcode07.fastq.gz.fasta.blast.tab.tsv", check.names = F, sep = "\t", header = F, row.names = 1)
temp_ASV[[8]] <- read.csv("ASV_TABLE_filtered_choped_barcode08.fastq.gz.fasta.blast.tab.tsv", check.names = F, sep = "\t", header = F, row.names = 1)
temp_ASV[[9]] <- read.csv("ASV_TABLE_filtered_choped_barcode09.fastq.gz.fasta.blast.tab.tsv", check.names = F, sep = "\t", header = F, row.names = 1)
temp_ASV[[10]] <- read.csv("ASV_TABLE_filtered_choped_barcode10.fastq.gz.fasta.blast.tab.tsv", check.names = F, sep = "\t", header = F, row.names = 1)
temp_ASV[[11]] <- read.csv("ASV_TABLE_filtered_choped_barcode11.fastq.gz.fasta.blast.tab.tsv", check.names = F, sep = "\t", header = F, row.names = 1)
temp_ASV[[12]] <- read.csv("ASV_TABLE_filtered_choped_barcode12.fastq.gz.fasta.blast.tab.tsv", check.names = F, sep = "\t", header = F, row.names = 1)

names(temp_ASV) <- barcodes

colnames(temp_ASV[[1]])


#The following function will deal with empty barcodes (like blanks)
for(i in 1:length(temp_ASV)) {
  if (lapply(temp_ASV[i], function(df) nrow(df)) == 0) {
    temp_ASV[i] <- list(data.frame(V1 = rep(0, times = 1)))
    rownames(temp_ASV[[i]]) <- rownames(temp_ASV[[i-1]])[nrow(temp_ASV[[i-1]])]
  }
}

for (i in 1:length(temp_ASV)){
  if (rownames(temp_ASV[[i]])[1] == "*"){
    rn <- rownames(temp_ASV[[i]])[-1]
    temp_ASV[[i]] <- data.frame(temp_ASV[[i]][-1,], check.rows = F, check.names = F)
    rownames(temp_ASV[[i]]) <- rn
  }
}

for (i in 1:length(temp_ASV)) colnames(temp_ASV[[i]]) <- barcodes[i]

names(temp_ASV) <- barcodes[1:length(temp_ASV)]

##ASV Taxonomy ----
#Individual taxonomy tables loading
temp_TAX = list.files(path = paste0("."), pattern="*Taxonomy.txt")
for (i in 1:length(temp_TAX)) assign(temp_TAX[i], data.frame(read.csv2(file = paste("./",temp_TAX[i], sep = ""), sep = ";", header = F, check.names = F, fill = TRUE, 
                                                                       col.names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", 
                                                                                     "Species", "other1", "other2", "other3", "other4", "other5", "other6",
                                                                                     "others7", "others8", "others9", "others10", "others11", "others12", "others13", "others14", "others15", "others16", "others17", "others19"))))
#The reason we allow so many fields is that Chloroplasts and Mitochondria (mainly), have non-coherent taxonomy sizes compared to Bacterial full taxonomy.
#To avoid any "Out of field" errors or alike, we allow them to exist in a first time

temp_TAX <- mget(temp_TAX)

for(i in 1: length(temp_TAX)){
  if (lapply(temp_TAX[i], function(df) nrow(df)) != 0) {#This step allows to correct for Kingdom, which is currently merged with SILVA ID
    rownames(temp_TAX[[i]]) <- sapply(strsplit(temp_TAX[[i]][,1], split = " "), "[[", 1)
    temp_TAX[[i]][,1] <- sapply(strsplit(temp_TAX[[i]][,1], split = " "), "[", 2)
    # if(file.exists(paste0(OUTPWD,"/Results/Unknown_clusters/unknown_clusters.tsv"))){
    #   temp_TAX[[i]] <- rbind(temp_TAX[[i]], U_TAX)
    # }
  }
}

for(i in 1: length(temp_TAX)) {
  if (nrow(temp_TAX[[i]]) == 0) {
    temp_TAX[[i]][1,] <- temp_TAX[[i-1]][1,]
    rownames(temp_TAX[[i]]) <- rownames(temp_ASV[[i-1]])[length(rownames(temp_ASV[[i-1]]))]
  }
}

temp_phyloseq <- barcodes

##Phyloseq objects ----

physeq_list <- list()

for (i in 1:length(temp_ASV)) {
  # Create the phyloseq object
  physeq_object <- phyloseq::phyloseq(phyloseq::otu_table(temp_ASV[[i]], taxa_are_rows = TRUE), 
                                      phyloseq::tax_table(as.matrix(temp_TAX[[1]])), 
                                      phyloseq::sample_data(metadata))
  
  # Get the name from the barcodes vector
  barcode_name <- barcodes[i]
  # Assign the phyloseq object to the list with the dynamic name
  physeq_list[[barcode_name]] <- physeq_object
}


i<-1 #Reset the incrementation
#Initialize the phyloseq object
NanoASV <- phyloseq::merge_phyloseq(physeq_list[[i]], physeq_list[[i + 1]])

# If more than 2 samples, then, adding them all together
if (length(physeq_list) > 2) {
  for (i in 3:length(physeq_list)) {
    NanoASV <- phyloseq::merge_phyloseq(NanoASV, physeq_list[[i]])
  }
}

Nygaard <- NanoASV
##Dataset cleaning ----
#Delete bad entries such as Eukaryota, Cyanobacteria and Archea if any
Nygaard <- phyloseq::subset_taxa(Nygaard, Kingdom != "Eukaryota")
Nygaard <- phyloseq::subset_taxa(Nygaard, Family != "Mitochondria")
Nygaard <- phyloseq::subset_taxa(Nygaard, Order != "Chloroplast")
##Taxonomy cleaning ----
#After those functions, there is no more taxa with mixed up names so we can remove supp fields of taxa table
tax_table(Nygaard) <- phyloseq::tax_table(Nygaard)[,1:7]

#Phyloseq export ----
print("Exporting results")
save(Nygaard, file = paste0("./Nygaard.rdata"))

#CSV export -----
#Create the TAXOTU file, encompassing taxonomy and abundance
TAXOTU <- data.frame(Nygaard@tax_table, Nygaard@otu_table)
#Exporting the file for people needing this format as well
write.csv(TAXOTU, file = "./Taxonomy-Abundance_table.csv")

