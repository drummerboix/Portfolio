#QUESTION----
#Find species richness from data collected by country from BOLD and and thorough investigating, I want to check whether two genes yield the same phylogenetic relationships. I will be using ribosomal RNA genes COI-5P and ITS to reconstruct phylogenic trees for the phylum Platyhelminthes. I want to determine whether the same clusters or grouping of particular species can be observed in the dendrograms generated  for these genes. I will also try to determine which one of these genes is better for understanding phylogenetic relationships in the phylum. i would also determine the top 10 species in the database and check which species are endemic to particular regions.
#LOAD NECESSARY PACKAGES ----
library(tidyverse)
library(vegan)
library(dendextend)
library(muscle)
library(ggplot2)
library(Biostrings)
library(ape)
library(DECIPHER)
library(ggmap)
library(phytools)

#CODE######################################################################################################################################################
#Download data from BOLD
#Platyhelminthes <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Platyhelminthes&format=tsv")

#Write data set into a file
#write_tsv(Platyhelminthes, "Platyhelminthes.tsv")
Platyhelminthes <- read_tsv("Platyhelminthes.tsv")

#Perform exploration of the data to determine columns suitable for analysis
class(Platyhelminthes)
summary(Platyhelminthes)
names(Platyhelminthes)
sum(is.na(Platyhelminthes))

#Extract columns useful for phylogenetic studies
Platyhelminthes_filtered <- Platyhelminthes %>%
  select(processid, bin_uri, species_name, country, lat, lon, markercode, nucleotides)

#Checking number of unique species
length(unique(Platyhelminthes_filtered$species_name))

#Checking numbers of unique bin
length(unique(Platyhelminthes_filtered$bin_uri))

#Checking number of unique genus
length(unique(Platyhelminthes_filtered$genus_name))

#Checking number of countries that provided data 
length(unique(Platyhelminthes_filtered$country))

#Remove records with missing data for species name, genus name, country, latitude, longitude, marker codes and missing nucleotides
Platyhelminthes_clean <- Platyhelminthes_filtered %>%
  filter(!is.na(bin_uri)) %>%
  filter(!is.na(country)) %>%
  filter(!is.na(lat)) %>%
  filter(!is.na(lon)) %>%
  filter(!is.na(markercode)) %>%
  filter(!is.na(nucleotides)) 

#Attribute check for cleaned up dataframe
dim(Platyhelminthes_clean)
names(Platyhelminthes_clean)
unique(Platyhelminthes_clean$species_name)
unique(Platyhelminthes_clean$bin_uri)
unique(Platyhelminthes_clean$country)
sum(is.na(Platyhelminthes_clean$nucleotides))

#Count records by marker code
Platyhelminthes_clean.marker <- Platyhelminthes_clean %>%
  group_by(markercode) %>%
  summarise(n = length(processid)) %>%
  arrange(desc(n)) %>%
  print()

#Count number of records within each country and return counts in descending order
Platyhelminthes_clean.country <- Platyhelminthes_clean %>%
  group_by(country) %>%
  summarise(n = length(processid)) %>%
  arrange(desc(n)) %>%
  print()

#Checking for minimum, max, median and mean sequence lengths
min(str_length(Platyhelminthes_clean$nucleotides))
mean(str_length(Platyhelminthes_clean$nucleotides))
median(str_length(Platyhelminthes_clean$nucleotides))
max(str_length(Platyhelminthes_clean$nucleotides))

########################MAP BIN DISTRIBUTION#########################################################################Source of map is Stamen.com
#Select columns necessary for quickmap operation using the ggmap package
mapdata <- Platyhelminthes_clean %>%
  select(country, lon, lat, bin_uri)

#Visualize a map for the country column, just an extra to verify longitude and latitude column data
qmplot(lon, lat, data = mapdata, colour = country, zoom = 5, main = "Worldmap showing all countries with submitted data")

#plot map using the BIN numbers as distinctions or differentiators with 5x zoom for easy download
WorldMapBin <- qmplot(lon, lat, data = mapdata, colour = bin_uri, zoom = 5, xlim = c(-166.573, 163), ylim = c(-34.897, 69.195), main = "Worldmap showing worldwide distribution of BIN")
WorldMapBin #The legend is too large

#I implemented a function to reduce size of legend because it was too large
Small_legends <- function(myPlot, pointSize = 3, textSize = 5, spaceLegend = 1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
Small_legends(WorldMapBin) #I applied the function to the map for better visibility

#For finding out how much data each country contributed
table(mapdata$country) #Canada has the most data, followed by the United States, then China, I will be plotting these on individual maps to compare distribution

#I want to analyze the 3 countries with most contributions to this data set and plot the exact region where the data was collected
#Canada first
CanadaMapdata <- mapdata %>%
  filter(country == "Canada")

CanadaMap <- qmplot(lon, lat, data = CanadaMapdata, colour = bin_uri, zoom = 5, main = "Canada BIN distribution") #I am going to be viewing the data according to BIN code, a particular color represents each BIN
Small_legends(CanadaMap) #I am applying the filter to the map because of the large sample size leading to the map being obstructed by the legends

#Then the united States 
USAMapdata <- mapdata %>%
  filter(country == "United States") #Filter by United States

USAMap <- qmplot(lon, lat, data = USAMapdata, colour = bin_uri, zoom = 5, main = "United States BIN distribution")
Small_legends(USAMap)

#Then China
ChinaMapdata <- mapdata %>%
  filter(country == "China") #Filter by china

ChinaMap <- qmplot(lon, lat, data = ChinaMapdata, colour = bin_uri, zoom = 5, main = "China BIN ditribution")
Small_legends(ChinaMap)

#For the exact amount of data each country contributed
BIN.by.country <- mapdata %>%
  group_by(bin_uri) %>% #BIN code in order 
  select(bin_uri, country) #i am selecting country and BIN code as they are important for analysis

#Check which BIN occurs more worldwide
END <- BIN.by.country %>%
  group_by(bin_uri) %>%
  count(country) #Count per country
  
colnames(END) <- (c("bin_uri","country","Occurences")) #Rename columns for clarity
END

#CONCLUSION: The data is very prevalent to the North American Region, most coming from Canada with the highest species diversity, most species of the Platyhelminthes appears to be endemic to the North American region while the rest is spread out through other continents, China seems to have the largest sample after this but it is made up of 4 main BIN. BOLD:AAC3222 seems to be the most occuring species in the world.
##########################################################################################################################################################################################################################################

#Checking distribution of sequence lengths
hist(str_length(Platyhelminthes_clean$nucleotides), xlab = 'Sequence Length', ylab = 'Frequency', main = 'Figure 1. Distribution of Platyhelminthes Sequence Lengths Frequency')
#The histogram shows stragglers in term of sequence length, remoove sequences <400 and >800, this would help eliminate outliers to enable better alignment further down the line 

#Filtering out sequences longer than 800 bp
Platyhelminthes_clean.trim <- Platyhelminthes_clean %>%
  filter(!(str_length(Platyhelminthes_clean$nucleotides) > 800))

#Filtering out sequences shorter than 400 bp
Platyhelminthes_clean.trim <- Platyhelminthes_clean.trim %>%
  filter(!(str_length(Platyhelminthes_clean.trim$nucleotides) < 400))
length(unique(Platyhelminthes_clean.trim$species_name)) #109

#Checks to make sure sequences >800 bp and <400 bp are removed
min(str_length(Platyhelminthes_clean.trim$nucleotides)) #400 
mean(str_length(Platyhelminthes_clean.trim$nucleotides)) #532.9897
median(str_length(Platyhelminthes_clean.trim$nucleotides)) #463
max(str_length(Platyhelminthes_clean.trim$nucleotides)) #795

#Checking the number of unique species after trimming
length(unique(Platyhelminthes_clean.trim$species_name)) #109

#Checking the number of unique bin after trimming
length(unique(Platyhelminthes_clean.trim$bin_uri)) #161

#Checking distribution of sequence length using an histogram
hist(str_length(Platyhelminthes_clean.trim$nucleotides),xlab = 'Sequence Length', ylab = 'Frequency', main = 'Figure 2. Distribution of trimmed Platyhelminthes Sequence Lengths Frequency')


####################################################################################################################################################CoI Phylo#######################################################################
#Filter for COI-5P with nucleotide data
Platyhelminthes_clean.trimCOI <- Platyhelminthes_clean.trim %>%
  filter(markercode == "COI-5P") %>%
  filter(str_detect(nucleotides, "[ACGT]"))

#Reformat sequences into DNAStringSet
Platyhelminthes_clean.trimCOI$nucleotides <- DNAStringSet(Platyhelminthes_clean.trimCOI$nucleotides)

#Check if reformatting was done correctly
class(Platyhelminthes_clean.trimCOI$nucleotides)

#Assign BIN to their nucleotides as labels throughout downstream analysis, 
#I am going to be using BIN code instead of species name because most of the data don't have proper species name, most are place holders
names(Platyhelminthes_clean.trimCOI$nucleotides) <- Platyhelminthes_clean.trimCOI$bin_uri

###############################################################ITS Phylo############################################################################
#Filter for ITS with nucleotide data
Platyhelminthes_clean.trimITS <- Platyhelminthes_clean.trim %>%
  filter(markercode == "ITS") %>%
  filter(str_detect(nucleotides, "[ACGT]"))

#Reformat sequences into DNAStringSet
Platyhelminthes_clean.trimITS$nucleotides <- DNAStringSet(Platyhelminthes_clean.trimITS$nucleotides)

#Check if reformatting was done correctly
class(Platyhelminthes_clean.trimITS$nucleotides)

#Assign BIN to their nucleotides as labels throughoit downstream analysis
names(Platyhelminthes_clean.trimITS$nucleotides) <- Platyhelminthes_clean.trimITS$bin_uri

###############################################################################################

#Multiple Sequence Alignments
#MSA for COI-5P####################################################################################

#Perform MSA using the function muscle from the muscle package
PlatCOI.alignment <- DNAStringSet(muscle::muscle(Platyhelminthes_clean.trimCOI$nucleotides), use.names = T)

#writeXStringSet(PlatCOI.alignment, file = "PlatCOIalignment.fasta") 

PlatCOI.alignment2 <- DNAStringSet(muscle::muscle(Platyhelminthes_clean.trimCOI$nucleotides, gapopen = -10000), use.names = T) #Test with tight penalty

#Check alignment on a web browser
BrowseSeqs(PlatCOI.alignment)

#Calculate mean, max, min and histogram displaying distribution of gaps
min(unlist(lapply(PlatCOI.alignment, str_count, "-"))) #112
max(unlist(lapply(PlatCOI.alignment, str_count, "-"))) #578
mean(unlist(lapply(PlatCOI.alignment, str_count, "-"))) #384.4246
median(unlist(lapply(PlatCOI.alignment, str_count, "-"))) #429
hist(unlist(lapply(PlatCOI.alignment, str_count, "-")), xlab = 'Number of Gaps', ylab = 'Frequency', main = 'Figure 3. Distribution of Gaps in COI Sequence MSA')

#Convert data type into DNABin for more downstream analysis, I am using as.DNAbin from the ape package
PlatCOI.DNABin <- as.DNAbin(PlatCOI.alignment)

#Create a distance matrix using the "K80" model, I am using the K80 model because its been a preferred model for analysis and the fact that it produces similar distance matrices to other methods, I am using pairwise deletion because it helps preserve most of the data and only values that are missing in all variables will be deleted, this will also be applied to the ITS analysis
PlatCOI.DM <- dist.dna(PlatCOI.DNABin, model = "K80", as.matrix = T, pairwise.deletion = T)

#I am going to be clustering using the UPGMA (average) method with a 2% divergence threshold according to Herbert et al., 2003 who used a value close to 3% for lepidopterans, which is also an invertebrate like the platyhelminthes, I am using the UPGMA because it overcomes the limitation of single linkage that connects distant relatives closely and complete linkage that has compact clustering
PlatCOI.Cluster <- IdClusters(PlatCOI.DM,
                              method = "UPGMA",
                              cutoff = 0.02, #2% sequence difference
                              showPlot = T, #returns dendrogram
                              type = "cluster", #returns clusters output
                              verbose = T)
title("Figure 4. Phylogeny of Platyhelminthes based on COI-5P rRNA")

#MSA for ITS##########################################################################################

#Perform MSA using the function muscle from the muscle package
PlatITS.alignment <- DNAStringSet(muscle::muscle(Platyhelminthes_clean.trimITS$nucleotides), use.names = T)

#writeXStringSet(PlatITS.alignment, file = "PlatITSalignment.fasta") #Write to file for bootstrapping

#Display MSA in browser
BrowseSeqs(PlatITS.alignment)

#Calculate mean, max, min, median and creating histogram of distribution of gaps
min(unlist(lapply(PlatITS.alignment, str_count, "-"))) #552
max(unlist(lapply(PlatITS.alignment, str_count, "-"))) #928
mean(unlist(lapply(PlatITS.alignment, str_count, "-"))) #680.2821
median(unlist(lapply(PlatITS.alignment, str_count, "-"))) #665
hist(unlist(lapply(PlatITS.alignment, str_count, "-")),xlab = 'Number of Gaps', ylab = 'Frequency', main = 'Figure 5. Distribution of Gaps in ITS Sequence MSA')

#Convert data type into DNAbin for more downstream analysis
PlatITS.DNABin <- as.DNAbin(PlatITS.alignment)

#Create distance matrix using "K80" model
PlatITS.DM <- dist.dna(PlatITS.DNABin, model = "K80", as.matrix = T, pairwise.deletion = T)

#Clustering using the UPGMA method with 2% divergence threshold
PlatITS.Cluster <- IdClusters(PlatITS.DM,
                               method = "UPGMA",
                               cutoff = 0.02, #2% sequence difference as species boundary
                               showPlot = T, #Generates dendrogram
                               type = "clusters", #Returns clusters output
                               verbose = T)
title("Figure 6. Phylogeny of Platyhelminthes based on ITS rRNA")

#The MSA for ITS produced much higher number of gaps compared to the COI data which means a lot of gaps have to be applied in order to create the best alignmment, this could mean that the ITS sequences contain lots of differences in gene sequences caused by deletions.

#The dendrogram generated is just too large and crowded and very hard to visualize properly, so I would be selecting one random BIN number to visualize and eventually create a tanglegram to give a more conclusive answer on the phylogenetic relationships of the Platyhelminthes.

######Extracting one record per bin code from COI#######################################

#Random selction of one representative per BIN
PlatCOI.Sample <- Platyhelminthes_clean.trimCOI %>%
  group_by(bin_uri) %>%
  sample_n(1)

#Checking number of unique BIN
length(unique(PlatCOI.Sample$bin_uri))

#Converting Sequences into DNAStringSet data type
PlatCOI.Sample$nucleotides <- DNAStringSet(PlatCOI.Sample$nucleotides)

#Setting BIN as identifier for each sequence
names(PlatCOI.Sample$nucleotides) <- PlatCOI.Sample$bin_uri

#Performing MSA on the sample
PlatCOI.Sample.Align <- DNAStringSet(muscle::muscle(PlatCOI.Sample$nucleotides), use.names = T)

#Display MSA on browser
BrowseSeqs(PlatCOI.Sample.Align)

#Convert alignment into DNAbin
PlatCOI.Sample.DNABin <- as.DNAbin(PlatCOI.Sample.Align)

#Create distance matrix using "K80" model
PlatCOI.Sample.DM <- dist.dna(PlatCOI.Sample.DNABin, model = "K80", as.matrix = T, pairwise.deletion = T)

#####Extracting one record per BIN from ITS###########################################################

#Random selection of one representative
PlatITS.Sample <- Platyhelminthes_clean.trimITS %>%
  group_by(bin_uri) %>%
  sample_n(1)

#Checking number of unique bin
length(unique(PlatITS.Sample$bin_uri))

#Convert sequences into DNAStringset
PlatITS.Sample$nucleotides <- DNAStringSet(PlatITS.Sample$nucleotides)

#Setting BIN as identifier for each sequence
names(PlatITS.Sample$nucleotides) <- PlatITS.Sample$bin_uri

#Perform MSA on each representative sequence
PlatITS.Sample.Align <- DNAStringSet(muscle::muscle(PlatITS.Sample$nucleotides), use.names = T)

#Display MSA on browser
BrowseSeqs(PlatITS.Sample.Align)

#Convert Alignment data into DNAbin
PlatITS.Sample.DNABin <- as.DNAbin(PlatITS.Sample.Align)

#Create distance matrix using "K80" model
PlatITS.Sample.DM <- dist.dna(PlatITS.Sample.DNABin, model = "K80", as.matrix = T, pairwise.deletion = T)

############################################################################################################################################################

#Visulalizing the phylogenetic trees for COI and ITS individually########################################

#Creating dendrogram for COI using the hclust function and using the set(den) function from the dendextend package to set parameters
PlatCOI.treehang <- PlatCOI.Sample.DNABin %>%
  dist %>%
  hclust(method = "average") %>% #clustering based on average linkage methodology
  as.dendrogram %>% #Convert cluster to a dendrogram
  set("branches_k_color",k=161) %>% #Coloring based on clusters for branches, i am coloring by unique BIN
  set("labels_cex", c(.6,.6)) %>% #resizing the label size
  set("nodes_pch", 19) %>% #Choosing shape for nodes
  set("nodes_col", c("black")) %>% #color of nodes
  set("hang_leaves") #creates hanging dendrogram
plot(PlatCOI.treehang)

#The dendrogram is too cluttered and packed, so i decided to create a circle dendrogram with the circlize_dendrogram function from the dendextend package.

PlatCOI.tree <- PlatCOI.Sample.DNABin %>%
  dist %>%
  hclust(method = "average") %>%
  as.dendrogram %>%
  set("branches_k_color",k=161) %>%
  set("labels_colors") %>%
  set("labels_cex", c(.6,.6)) %>%
  set("nodes_pch", 19) %>%
  set("nodes_col", c("black"))
circlize_dendrogram(PlatCOI.tree)
title("Figure 7. UPGMA tree for COI")

#I removed tip labels for cleaner visualization
circlize_dendrogram(PlatCOI.tree, labels = F)
title("Figure 7. UPGMA tree for COI")

########################################################################################################

#Creating dendrogram for ITS
PlatITS.tree <- PlatITS.Sample.DNABin %>%
  dist %>%
  hclust(method = "average") %>%
  as.dendrogram %>%
  set("branches_k_color",k=24) %>%
  set("labels_colors") %>%
  set("labels_cex", c(.6,.6)) %>%
  set("nodes_pch", 19) %>%
  set("nodes_col", c("black"))
circlize_dendrogram(PlatITS.tree)
title("Figure 8. UPGMA tree for ITS")
circlize_dendrogram(PlatITS.tree, labels = F)
title("Figure 8. UPGMA tree for ITS")

#Comparing dendrograms generated by COI and ITS rRNA genes##################################################################################################

#I want to compare the two trees using the tanglegram function from the dendextend package, this function creates a mirror image and prunes the branches that the COI and the ITS doo not share and generates a comparison of the available variables

#I will be creating a dendrogram using the COI MSA distance matrix
PlatCOI.dend <- PlatCOI.Sample.DNABin %>%
  dist %>%
  hclust(method = "average") %>% #Making clusters based on the average-linkage method
  as.dendrogram %>% #Converts the clusters to dendrogram
  set("branches_k_color", k=161) #Setting the colors of the branches based on clusters, i am using unique BIN as color

#Phylogenetic signal
COIphylo <- as.phylo(PlatCOI.dend) #Convert to class phylo
x <- fastBM(COIphylo)
phylosig(COIphylo,x,method = "lambda", test = T)

PlatITS.dend <- PlatITS.Sample.DNABin %>%
  dist %>%
  hclust(method = "average") %>% #Making clusters based on the average-linkage method
  as.dendrogram %>% #Converts the clusters to dendrogram
  set("branches_k_color", k=24) #Setting the colors of the branches based on clusters, i am using unique BIN as color

ITSphylo <- as.phylo(PlatITS.dend) #Convert to class phylo
y <- fastBM(ITSphylo)
phylosig(ITSphylo,x,method = "lambda", test = T)

#Compare the dendrograms from the COI and ITS rRNA using tanglegram from the dendextend package
tanglegram(PlatCOI.dend, PlatITS.dend, main = "Figure 9. Phylogeny Reconstruction", main_left = "Platyhelminthes ITS rRNA", main_right = "Platyhelminthes COI-5P rRNA", common_subtrees_color_branches = T, columns_width = c(5, 2, 5), margin_inner = 11, cex_main = 1, cex_main_left = 1, cex_main_right = 1, sort = F, rank_branches = T)

########CONCLUSION##########################################################################################################################################
#The COI and the ITS rRNA do not seem to produce the same phylogenetic hypothesis. An observation can be made that the BIN numbers are placed in different clusters in Figure 9. In Figures 4 and 6, we can see that the distance scales are very different, the scale in COI is more than four times that of ITS dendrogram. The difference might be due to the differences in the available data considering there was more data for the COI compared to the ITS, the clustering and the DNA distance methods might not have been the best fit for the data analysis. Some adjustments could be made to improve results, we can use RAxML for bootstrapping the alignments and performing cophoenetic correlation tests on all the dendrograms to check their validity, trying all the distance models and tests could help in deciding which methods could have been the best choice. Collecting data from other sources will probably help to improve the study, a CADM test could also be useful.
