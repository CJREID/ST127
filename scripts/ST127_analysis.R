####-----PACKAGES-----####
# Load required packages
library(tidyverse)
library(ggpubr)
library(reshape2)
library(pheatmap)
library(ggtree)
library(RColorBrewer)
library(abricateR)
library(gridExtra)
library(paletteer)

####-----SETUP-----####
# SET WORKING DIRECTORY
wrkdir <- "/Volumes/126451/WORK/Colleagues/Paarthiphan/ST127/revision/ST127_project/"
setwd(wrkdir)

# Define 'not in' function for subsetting
'%notin%' <- Negate('%in%')

# Make output directory
if (dir.exists("outputs")){
} else {
  dir.create("outputs")
  dir.create("outputs/data")
  dir.create("outputs/figures")
}

####-----IMPORT DATA-----####
#Get input filenames
infiles <-c(list.files("ST127_pipelord2_results/summaries/", pattern = "\\.txt"))

#Read in and name data.frames with whatever comes before .txt
for (file in infiles){
  indir <- c("ST127_pipelord2_results/summaries")
  f <- read_delim(paste(indir, file, sep = "/"), delim = "\t", col_names = TRUE, trim_ws = TRUE)
  assign(paste(substr(file, 1, nchar(file)-4), sep = ""), f)
}

# Read in raw metadata without accessions
metadata <- read_delim("meta/ST127_Meta_edited.csv", delim = ",") %>% 
  select(everything(), -`Short Read Accession`, -BioProject, -BioSample) %>%
  mutate(`Collection Date` = as.character(`Collection Date`))

# Read in raw metadata with accessions (used later to make Table S1)
all_meta <- read_delim("meta/ST127_Meta_edited.csv", delim = ",") %>% 
  mutate(`Collection Date` = as.character(`Collection Date`))

####-----PROCESS ABRICATE DATA-----####
#Load paths to files needed for abricateR
abricate_path <- "ST127_pipelord2_results/summaries/genotype.txt"
pointfinder_path <- "ST127_pipelord2_results/summaries/pointfinder.txt"
pMLST_data <- "ST127_pipelord2_results/summaries/pMLST.txt"

#Provide output names
output_name <- "ST127"
output_dir <- "processed_R"

#Run abricateR
abricateR::abricateR(
  abricate_in = abricate_path,
  output = output_name,
  identity = 90,
  length = 95,
  output_directory = output_dir,
  writecsv = FALSE,
  pointfinder_data = pointfinder_path,
  pMLST_data = pMLST_data
)

####-----PROCESS pUTI89 ALIGNMENT-----####
# Get tree path and abricate alignment data
path_to_tree <- "trees/CGA_full.prokka.og.rooted.removed.contree"
path_to_abricate <- "ST127_pipelord2_results/abricate/plasmids/pUTI89/ST127.pUTI89.abricate.tab"
plasrefname <- "pUTI89"

# Minimum hit thresholds
# Minimum hit length (as a percentage [i.e. 0.5 = 0.5%])
min_hit_length <- 0.5
# Minimum nucleotide ID (also as a percentage [i.e. 90 = 90%])
min_hit_id <- 90

# Read in the abricate genotype data sheet
#(small number of rows for colname reassignment)
#This is to reduce memory requirements
abricate_hits <-
  read_delim(
    path_to_abricate,
    "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    n_max = 10
  )

# Rename some columns for clarity
colnames(abricate_hits)[c(1, 10:11)] <-
  c("name", "perc_coverage", "perc_identity")

# Extract column names for later reassignment
abricate_hits_colnames <- colnames(abricate_hits)

# Read in full abricate data sheet
abricate_hits <-
  read_delim(
    path_to_abricate,
    "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    col_names = FALSE,
    skip = 1
  )

# Remove cases where there are multiple headers from concatenation of abricate reports
abricate_hits <- abricate_hits %>% filter(X2 != "SEQUENCE")

# Assign column names
colnames(abricate_hits) <- abricate_hits_colnames

# Convert percent coverage and identity to numeric type to allow filtering
abricate_hits$perc_coverage <- as.numeric(abricate_hits$perc_coverage)
abricate_hits$perc_identity <- as.numeric(abricate_hits$perc_identity)

#Filter to perc_identity > 90%
#abricate_hits <-
abricate_hits <- abricate_hits %>% filter(perc_identity > min_hit_id)
abricate_hits <- abricate_hits %>% filter(perc_coverage > min_hit_length)

# Trim excess characters in the assembly names and reassign this to row names
abricate_hits$name <- gsub("\\..*", "", abricate_hits$name)

#Read in the tree file
tree <-
  read.tree(file = path_to_tree)

# Trim the names of the assemblies in the tree tip labels
tree$tip.label <- gsub("\\..*", "", tree$tip.label)

# Subset the hits to strains within the tree
abricate_hits <- abricate_hits %>% filter(name %in% tree$tip.label)

# Extract from coverage column the coverage range and the length of the reference
abricate_hits$coverage_range <- gsub("\\/.*","", abricate_hits$COVERAGE)
abricate_hits$ref_length <- gsub(".*\\/","", abricate_hits$COVERAGE)

# Save length of plasmid reference as a variable for later use
ref_length <- as.numeric(unique(abricate_hits$ref_length))

# Replace the '-' with a ':' in the coverage
abricate_hits$coverage_range <- gsub("-",":", abricate_hits$coverage_range)

# Create a column for start and end coordinates of hits
abricate_hits$end <- gsub("[0-9]+:","", abricate_hits$coverage_range)
abricate_hits$start <- gsub(":[0-9]+","", abricate_hits$coverage_range)

# Select columns of interest
abricate_hits <- abricate_hits %>%
  select(name, gene = GENE, ref_length, start, end, percentage = perc_coverage)

# Convert start and end coordinates to numeric
abricate_hits$start <- as.numeric(abricate_hits$start)
abricate_hits$end <- as.numeric(abricate_hits$end)

# Create an empty matrix equal to length of ref plasmid
empty_plasrow <- rep(0, times = unique(abricate_hits$ref_length))

# Create an empty matrix with n rows (n = sample size) with ncol == length(ref plasmid)
empty_plasmatrix <- matrix(rep(empty_plasrow,
                               times = length(unique(abricate_hits$name))),
                           nrow = length(unique(abricate_hits$name)))

# Create a list of levels for sample names, a list of start coords and a list of end coords
# and bind these in a list of lists
start_ends <- list(as.list(as.integer(as.factor(abricate_hits$name))),
                   as.list(as.integer(abricate_hits$start)),
                   as.list(as.integer(abricate_hits$end)))

# Create a counter
counter <- 0

# Map the BLAST hits to our matrix of bp coordinates
for (i in 1:nrow(abricate_hits)){
  sample <- start_ends[[1]][[i]]
  start_coord <- start_ends[[2]][[i]]
  end_coord <- start_ends[[3]][[i]]
  empty_plasmatrix[sample, start_coord:end_coord] <- 1
  counter <- counter + 1
  if(counter %% 1000 == 0){
    message(paste(counter, "out of ", nrow(abricate_hits), "hits processed"))}
}

# Rename matrix
base_matrix <- empty_plasmatrix

# Remove old matrix
rm(empty_plasmatrix)

# Convert matrix to a dataframe
base_matrix <- as.data.frame(base_matrix, stringsAsFactors = FALSE)

# Assign sample names to rows
rownames(base_matrix) <- unique(abricate_hits$name)

# Generate indices that cover blocks of 100 columns in base_matrix
bin_ranges <- c(seq(from = 1, to = ref_length, by = 100))
bin_ranges2 <- c(seq(from = 100, to = ref_length, by = 100), ref_length)

# Split the ranges into two lists of [1] start and [2] end indices 
bin_splits <- list(bin_ranges, bin_ranges2)

# Initialise empty vector for loop below
binned_hits <- vector()

#Create a counter
counter <- 0

# Binning loop - this loop creates a matrix of percentage hits for 100 bp chunks of the reference plasmid for each strain
for (i in 1:length(bin_ranges2)){
  # If the last bin is only 1 base long then the as.matrix line won't work,
  #so we have to include the if statement below:
  if(i == length(bin_ranges) & bin_ranges[length(bin_ranges)] == bin_ranges2[length(bin_ranges2)]){
    row_sum <- 1
  }else{
    # Generate row sums (i.e. Number of matching bases) for 100 column chunks of the base_matrix
    row_sum <- as.matrix(rowSums(base_matrix[,bin_splits[[1]][i]:bin_splits[[2]][i]]))
  }
  # Bind them together in a new vector
  binned_hits <- cbind(binned_hits, row_sum)
  counter <- counter + 1
  if(counter %% 50 == 0){
    message(paste(counter, " samples out of ", length(bin_splits[[1]]), " processed" ))}
}

# Save row names
nems <- rownames(binned_hits)

# Convert binned_hits to a data frame
binned_hits <- as.data.frame(binned_hits)

# Assign row names
rownames(binned_hits) <- nems

# Generate dataframe wiht a summary statistic for how many bins were covered at >90%ID
hits_df <- as.data.frame(rowSums(binned_hits))

# Add names column
hits_df$working_names <- rownames(hits_df)

# Rename columns
colnames(hits_df) <- c("pUTI89_ID","working_name")

# Convert value to a percentage
hits_df$pUTI89_ID <- round((hits_df$pUTI89_ID/ref_length) * 100)

####-----PROCESS METADATA-----####
#Select ColV data to add to meta
colv <- ST127_simple_summary_N90L95 %>% select(Name = name, ColV, IncF_RST, IncI1_MLST)

# Compile metadata
metadata <- left_join(metadata, mlst, by = c("Name" = "name")) %>% mutate(ST = str_extract(ST, "^[0-9]{3}")) %>% select(-scheme)
metadata <- left_join(metadata, colv) %>% select(-IncI1_MLST)

# Categorise ColV
metadata <- metadata %>% mutate(ColV = case_when(ColV == "0" ~ "No", ColV == "1" ~ "Yes"))

# Categorise seqeunces that are F29:A-:B10 carriers or >90% coverage of pUTI89 as pUTI89+
metadata <- left_join(metadata, hits_df, by = c("Name"="working_name")) %>%
  mutate(pUTI89 = case_when(pUTI89_ID >= 78 | IncF_RST == "F29:A-:B10" ~ "Yes", TRUE ~ "No")) %>%
  select(-pUTI89_ID)

# Extract working names for collection
working_names <- as.vector(metadata$Name)

# Create 'geno_meta' df with all metadata and gene screening data
geno_meta <- left_join(metadata, ST127_simple_summary_N90L95 %>% select(-ColV, -IncF_RST, -IncI1_MLST), by = c("Name" = "name"))

# Define gene columns as those that are integers
gene_cols <- names(geno_meta %>% select(where(is.integer)))

# Recode multiple hits as a single hit
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(2, 1, .x)))
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(3, 1, .x)))
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(4, 1, .x)))

# Convert back to integer (gsub makes everything character class)
geno_meta <- geno_meta %>% mutate(across(all_of(gene_cols), as.integer))

# Extract metadata column names for later use with the geno_meta dataframe
meta_cols <- names(metadata)

####-----BAPS-----####
# Read in BAP cluster data
bap_gene <- read_csv("fastbaps/coregene_MVC107.fastbaps_clusters.csv") %>% select(Name = Isolates, everything())

# Extract names
strain_names <- bap_gene$Name

# Re-assign names to rownames
rownames(bap_gene) <- strain_names

# Generate data frame for plotting panels for BAP clustering of each alignment
plot.df.gene <- data.frame(id = rownames(bap_gene), 
                           fastbaps1 = bap_gene$`Level 1`, 
                           stringsAsFactors = FALSE)

# Create a character column of the values to make plotting easier
plot.df.gene$gene_bap <- as.character(plot.df.gene$fastbaps1)

# Add BAP prefix to the cluster numbers
plot.df.gene$gene_bap <- paste0("BAP", plot.df.gene$gene_bap)

# Process data into casted format for Scoary analysis prior to excluding the outgroup strain
baps.scoary <- plot.df.gene

# Create value column for casting
baps.scoary$value <- rep(1, nrow(plot.df.gene))

# Select relevant columns
baps.scoary <- baps.scoary %>% select(Name = id, Cluster = gene_bap, value)

# Cast into wide format
baps.scoary <- baps.scoary %>% reshape2::dcast(Name ~ Cluster)

# Recode NAs to zeroes
baps.scoary[is.na(baps.scoary)] <- 0

# Write file for Scoary
# write_csv(baps.scoary, "ST127.BAPS.scoary.csv")

# Remove outgroup strain from BAP data
plot.df.gene <- plot.df.gene %>% filter(id != "MVC107")

# Create ggtrtee object from ST127 phylogenetic tree
prgo.bap.tree <- ggtree(tree, branch.length = "none") 

# Create facet plots showing fastbaps grouping for each phylogeny
prgo.facet_gene <- facet_plot(prgo.bap.tree, 
                         panel = "fastbaps",
                         data = plot.df.gene,
                         geom = geom_raster, 
                         aes(x = fastbaps1, 
                             fill = gene_bap))

####-----TREE METADATA-----####
# Add BAP cluster data to metadata
metadata <- left_join(metadata, plot.df.gene, by = c("Name" = "id")) %>% select(everything(), -fastbaps1, Cluster = gene_bap)

# Extract just what we need for tree visualisation
tree_meta <- metadata %>% 
  select(Name,
         Cluster,
         Source, 
         Date = `Collection Date`, 
         Continent, 
         Country,
         ColV,
         IncF_RST,
         pUTI89)

# Get names
tree_names <- tree_meta$Name

# Refine tree metadata
tree_meta <- tree_meta %>% select(Cluster, Source, Continent, pUTI89)

# Assign names as row names (for gheatmap function)
rownames(tree_meta) <- tree_names

####-----DEFINE COLOURS-----####
# pUTI89 colours
puti_vars <- unique(tree_meta$pUTI89)
puti_clrs <- c("#38c9b1", "white")
names(puti_clrs) <- rev(puti_vars)

# Colours for BAP clusters
bap_vars <- unique(tree_meta$Cluster)
bap_clrs <- paletteer_dynamic("cartography::multi.pal", 13)
names(bap_clrs) <- sort(bap_vars)

# Colours for isolate source
source_vars <- unique(tree_meta$Source)
source_clrs <- colorRampPalette(brewer.pal(7, "Set1"))(7)
names(source_clrs) <- sort(source_vars)
edit <- c("#f0f035", "#FF7F00")
names(edit) <- c("Human", "Livestock")
source_clrs[names(edit)] <- edit

# Colours for continents
continent_vars <- unique(tree_meta$Continent)
continent_clrs <- colorRampPalette(brewer.pal(5, "Set2"))(5)
names(continent_clrs) <- sort(continent_vars)

# Combine colours for gheatmap
tree_vars <- c(source_clrs, continent_clrs, bap_clrs, puti_clrs)

####-----FIGURE 1 SOURCE SUMMARY-----####
# Plot isolate continents of origin stratified by source
figure1 <- ggplot(metadata, aes(Continent)) +
  geom_bar(aes(fill = Source))+
  scale_fill_manual(name = "Source", values = source_clrs) +
  scale_x_discrete(name = "Continent", labels = function(x) str_wrap(x, width = 10))+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 200), n.breaks = 8) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 12, vjust = , face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 16, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size =12),
        legend.text = element_text(size = 10)) 

# Save as Figure 1
ggsave("Figure1_metadata.png",
       figure1, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 150)

####-----FIGURE 2 TREE-----####
# Create ggtree object from phylogenetic tree
prgo.tree.heat <- ggtree(tree, layout = "fan", open.angle = 7, size = .2)

# Draw main tree wiht metadata
figure2 <- gheatmap(prgo.tree.heat, 
         tree_meta, 
         width = .2,
         font.size = 1.7,
         colnames_offset_x = ,
         colnames_offset_y = 2.4,
         colnames_position = "top",
         colnames_angle = 88,
         hjust = 0.5,
         color = rgb(0, 0, 0, alpha = .2)) +
  scale_fill_manual(name = "Data", values = tree_vars) + 
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 10))

# Save as Figure 2
ggsave("Figure2_tree.png",
       figure2, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

####-----FIGURE 3 SNP HEATMAP-----####
# The purpose of this section is to visualise strain pairs from different sources that differ by less than 30 SNPs 
# across the 3,266,764 bp core gene alignment with a focus on the dominant human and companion animal isolates

# Read in pairwise SNP-dists file
core_snps <- read_csv(file = "snps/ST127_coreSNPdists.csv") %>% rename(Name = ...1)

# Extract names
rows <- core_snps$Name

# Melt into long format
melt.core <- reshape2::melt(core_snps)

# Join metadata to first isolate in each pair
melt.core.a <- left_join(melt.core, metadata, by=c("Name")) %>% 
  select(Name, variable, value, Source, Country, `Collection Year` = `Collection Date`, Cluster) %>% 
  rename(Name1 = Name,
         Name = variable,
         Source1 = Source,
         Country1 = Country,
         Clustera = Cluster,
         `Collection Year1` = `Collection Year`)

# Join metadata to second isolate in each pair, then filter to equal to or less than 30 SNPs. Exclude self comparisons by filtering out 0 SNPS
# Filter out identical source pairs
# Filter Source 1 to Human and Companion animal and Source two to anything but Companion animals
# This results in 57 unique pairs of human or companion animal vs all comparisons
melt.core.unique <- left_join(melt.core.a, metadata, by=c("Name")) %>%
  select(Name1, Name, Source1, Source, Country1, Country, `Collection Year1`, `Collection Year`= `Collection Date`, Clustera, Cluster, value)%>% 
  rename(Name2 = Name, Source2 = Source, Country2 = Country, `Collection Year2` = `Collection Year`, Clusterb = Cluster, `Core SNPs` = value) %>% 
  filter(`Source1` != `Source2`, `Core SNPs` <=30, `Core SNPs` != 0, Source1 == "Human" | Source1 == "Companion Animal", Source2 != "Companion Animal")

# Extract the countries for visualisation later
country1 <-  melt.core.unique %>% filter(`Core SNPs` <= 30) %>% pull(Country1) %>% unique()
country2 <-  melt.core.unique %>% filter(`Core SNPs` <= 30) %>% pull(Country2) %>% unique()
snp_countries <- unique(c(country1, country2))

# Extract annotation row data for pheatmap
anno_rows.unique <- melt.core.unique %>% 
  filter(`Core SNPs` <= 30) %>% 
  select(Name1, 
         Source = Source1,
         Country = Country1,
         Cluster = Clustera) %>%
  distinct()

# Convert Name column to rownames
rownem.unique <- anno_rows.unique$Name1
anno_rows.unique <- anno_rows.unique %>% select(-Name1)
rownames(anno_rows.unique) <- rownem.unique

# As above except for columns
anno_cols_snp.unique <- melt.core.unique %>% 
  filter(`Core SNPs` <= 30) %>%
  select(Name2, 
         Source = Source2,
         Country = Country2,
         Cluster = Clusterb) %>% 
  distinct()

# Convert name to rownames
colnem.unique <-anno_cols_snp.unique$Name2
anno_cols_snp.unique <- anno_cols_snp.unique %>% select(-Name2)
rownames(anno_cols_snp.unique) <- colnem.unique

# Create DF of paired names and SNP counts
pairs.unique <- melt.core.unique %>% select(Name1, `Core SNPs`, Name2) %>% filter(`Core SNPs` <= 30)

# Cast to SNP matrix
paircast.unique <- reshape2::dcast(pairs.unique, Name1 ~ Name2, value.var = "Core SNPs")

# Move Name column to rownames and convert to matrices
SNP_rows.unique <- paircast.unique$Name1
paircast_mat.unique <- as.matrix(paircast.unique %>% select(-Name1))
rownames(paircast_mat.unique) <- SNP_rows.unique

# Set NAs outside the SNP range so they can be ignored
paircast_mat.unique[is.na(paircast_mat.unique)] <- 100

# Make quick heatmap with clustering to group the most closely related strains
out.unique <- pheatmap(mat = paircast_mat.unique,
                  cluster_cols = TRUE,
                  cluster_rows = TRUE)

# Extract row and column orders from the above heatmap
roword.unique <- rownames(paircast_mat.unique[out.unique$tree_row[["order"]],])
colord.unique <- colnames(paircast_mat.unique[,out.unique$tree_col[["order"]]])

# Reorder matrix with the order extracted above
paircast_mat.unique <- paircast_mat.unique[roword.unique,]
paircast_mat.unique <- paircast_mat.unique[,colord.unique]

# Reset 100 values as NAs
paircast_mat.unique[paircast_mat.unique == 100] <- NA

# Define country colours for heatmap annotation
country_clrs <- colorRampPalette(brewer.pal(8, "Set1"))(9)
names(country_clrs) <- snp_countries

# Combine all colours for heatmap annotation
anno_colors <- list(Source = source_clrs,
                    Country = country_clrs, Cluster = bap_clrs)

# Create Figure 3
figure3 <- pheatmap(mat = paircast_mat.unique,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_col = anno_cols_snp.unique,
         annotation_row = anno_rows.unique,
         color = colorRampPalette(brewer.pal(9, "BuPu"))(14),
         breaks = c(seq(from = 0, to = 30, by = 2)),
         fontsize_col = 8,
         fontsize_row =8,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_colors = anno_colors,
         annotation_names_col = TRUE,
         annotation_names_row = TRUE,
         legend = TRUE,
         annotation_legend = FALSE,
         border_color = "black",
         na.cols = "grey")

# Save figure 3
ggsave("Figure3_SNPheatmap_all.png", 
       figure3, 
       path = "outputs/figures/", 
       device = "png", 
       width = 297,
       height = 297,
       unit ="mm", 
       dpi = 150)

####-----FIGURE 4 pUTI89 HEATMAP-----####
# Due to the density of 'columns' arising from the hit matrix in the PROCESS pUTI89 ALIGNMENT section, this figure had to be visualised in two parts and combined
# The first part is the tree with metadata and the second is the map of hits to the pUTI89 sequence - please see the publication figure.

# Create ggtree object for pUTI89 alignment visualisation
pUTI_tree <- ggtree(tree, branch.length = "none") %<+% metadata

# Refine meta data for use in this figure
meta_heat <- metadata %>% select(Name, Cluster, Source, pUTI89)
metaheatnames <- meta_heat$Name
meta_heat <- meta_heat %>% select(Cluster, Source, pUTI89)
rownames(meta_heat) <- metaheatnames

# Define colours
alt_ptui89_clrs <- c("Yes" = "#df03fc","No" = "white")

# Visualise the tree with metadata
figure4a <- gheatmap(pUTI_tree, 
                      meta_heat,
                      offset = ,
                      width = .1,
                      font.size = 2,
                      colnames= TRUE,
                      colnames_position = "top",
                      colnames_angle = 70,
                      colnames_offset_y = 5,
                      hjust = .5,
                      color = rgb(0, 0, 0, alpha = 0)) +
  scale_fill_manual(name = "Data", values = c(alt_ptui89_clrs, bap_clrs, source_clrs))+
  ylim(c(0,305))+
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

# Visualise the tree with hits to pUTI89
figure4b <- gheatmap(pUTI_tree, 
                      binned_hits,
                      offset = 35,
                      width = 30,
                      font.size = 2,
                      colnames= FALSE,
                      hjust = 0,
                      color = rgb(0, 0, 0, alpha = 0)) +
  scale_fill_gradient(low = "white", high = "#8dd3c7", na.value = "white")+
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

# Save Figure 4a
ggsave("Figure4_pUTI89_heatmapA.png",
       figure4a, 
       path = "outputs/figures/", 
       device = "png", 
       width= 130, 
       height = 210, 
       unit ="mm", 
       dpi = 300)

# Save Figure 4b
ggsave("Figure4_pUTI89_heatmapB.png",
       figure4b, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 300)

####-----HEATMAP DATA PROCESSING-----####
# Split ABRicate gene hits into their functional groups
# Get all the hits from CARD database and intI1 and intI2. Fix all the messy names. Filter out genes present in >90% of isolates. These are housekeeping
# genes that sometimes mutate to confer AMR phenotypes but we are only concerned wiht acquired resistance genes
args <- geno_meta %>%
  select(all_of(meta_cols), starts_with("card"), contains("intI")) %>%
  rename_with(~ gsub("card_", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("Escherichia_coli_", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("_beta-lactamase", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("(", "_", .x, fixed = TRUE)) %>%
  rename_with(~ gsub(")", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("PC1__", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("EC_custom_intI1.*", "intI1", .x)) %>%
  rename_with(~ gsub("EC_custom_intI2.*", "intI2", .x)) %>%
  rename_with(~ gsub("Shigella_flexneri_chloramphenicol_acetyltransferase", "catA1", .x, fixed = TRUE)) %>% 
  select(where(is.character), where( ~ is.integer(.x) && sum(.x) <= 270)) %>%
  select(sort(names(.))) %>%
  relocate(all_of(meta_cols), contains("intI"))

# Calculate total carriage of each gene in the collection
arg_totals <- t(args %>% summarise(across(where(is.integer), sum)))
arg_totals <- as_tibble(arg_totals, rownames = "Gene") %>% rename(Total = V1) %>% mutate(Percentage = round(Total/299*100, 2))

# Get hits from VFDB
vags <- geno_meta %>% 
  select(all_of(meta_cols), starts_with("vfdb")) %>%
  rename_with(~ gsub("vfdb_", "", .x))

# Selett additional virulence genes from our custom database
custom_vags <- geno_meta %>% 
    select(Name, contains("EC_custom")) %>%
  rename_with(~ gsub("EC_custom_", "", .x, fixed = TRUE)) %>%
  select(Name,
         starts_with(c("cba", "cbi", "cjr", 
                       "cva", "cvi","eit",
                       "fecA", "hek", "hyl",
                       "iha","iss", "merA",
                       "ompT", "silA",
                       "terA", "traT", "usp"))) %>%
  rename_with(~ gsub("_[A-Z]{1,2}.*", "", .x)) %>%
  rename_with(~ gsub("_pAPEC-O1-ColBM", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_pUTI89", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_type3", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_|VFG1539", "", .x, fixed =TRUE))

# Join VFDB and custom hits
vags <- left_join(vags, custom_vags) %>% select(sort(everything())) %>%
  relocate(all_of(meta_cols))

# Calculate total carriage of each gene in the collection
vag_totals <- t(vags %>% summarise(across(where(is.integer), sum)))
vag_totals <- as_tibble(vag_totals, rownames = "Gene") %>% rename(Total = V1) %>% mutate(Percentage = round(Total/299*100, 2))

# Extract mobile genetic elements from the ISFinder database. NB: this data is not covered in the paper
mges <- geno_meta %>% 
  select(all_of(meta_cols), starts_with("ISfinder_")) %>%
  rename_with(~ gsub("ISfinder_Feb_2020_", "", .x)) %>%
  rename_with(~ gsub(":.*", "", .x))

# Extract plasmid related genes from the plasmidfinder database and tidy up messy names
plas <- geno_meta %>% 
  select(Name, starts_with("plasmidfinder")) %>%
  rename_with(~ gsub("plasmidfinder_", "", .x)) %>%
  mutate(IncBOKZ = rowSums(select(.,starts_with("IncB")))) %>%
  select(-starts_with("IncB/"), -!matches("Inc|Col"), Name) %>%
  rename_with(~ gsub("_.*$", "", .x)) %>%
  rename_with(~ gsub("(", "_", .x, fixed = TRUE)) %>%
  rename_with(~ gsub(")", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("/", "", .x, fixed = TRUE)) %>% 
  mutate(IncBOKZ = as.integer(IncBOKZ))

# Calculate total carriage of each gene in the collection
plas_totals <- t(plas %>% summarise(across(where(is.integer), sum)))
plas_totals <- as_tibble(plas_totals, rownames = "Gene") %>% 
                        rename(Total = V1) %>% 
                        mutate(Percentage = round(Total/299*100, 2))

# Summary of IncF RSTs in the collection
F_RST_totals <- geno_meta %>% 
  filter(IncF_RST != "F-:A-:B-") %>% 
  group_by(IncF_RST,ColV, pUTI89) %>%
  summarise(Total = n()) %>% 
  mutate(Percentage = round(Total/299*100, 2))

# Put dataframes into a list so they can be converted to matrices for heatmap visualisation
gene_list <- list(args = args, vags = vags, mges = mges, plas = plas)

# Loop that converts the gene screening dataframes into binary heatmaps named `geneprefix_heat` e.g args heatmap is called args_heat
for (f in 1:length(gene_list)){
  heat <- gene_list[[f]]
  heat <- heat %>% select(where(is_integer)) %>% as.matrix()
  rownames(heat) <- gene_list[[f]]$Name
  heat[is.na(heat)] <- 0
  heat[heat >= 1] <- "Present"
  heat[heat == 0] <- "Absent"
  heat <- as.data.frame(cbind(heat, pUTI89 = metadata$pUTI89, Cluster = metadata$Cluster, Source = metadata$Source))
  heat <- heat %>% select(Cluster, Source, pUTI89, everything())
  assign(paste(names(gene_list[f]), "heat", sep="_"), heat)
}

# Additional colours for visualisation
heat_clrs <- c("Present" = "#8dd3c7", "Absent" = "#ededed", "Yes" = "#df03fc","No" = "white")

####-----FIGURE 5 PUTI89 vs ARGs and VAGs-----####
# Add a column counting total resistance genes wihtout counting integrase genes
args2 <- args %>% select(-intI1, -intI2) %>% mutate(`Total ARGs` = rowSums(across(where(is.integer))))

# Convert 1s and 0s to Yes/No 
args2 <- args2 %>% mutate(across(`AAC_3-IId`:tet_B, ~ gsub(1, "Yes", .x)), across(`AAC_3-IId`:tet_B, ~ gsub(0, "No", .x)))

# Plot Total ARGs by pUTI89 carriage
pUTI89_argplot <-  ggboxplot(args2, x = "pUTI89", y = "Total ARGs",
                           color = "pUTI89", palette = c("#48e0c7", "#f52020"),
                           order = c("Yes", "No"),
                           ylab = "ARGs", xlab = "pUTI89") +
  scale_x_discrete(label = c("Yes", "No"))

# Kruskal test to determine if there is a difference between means for the categories
kruskal.test(`Total ARGs` ~ pUTI89, data = args2)

# Wilcoxon test to quantify the difference with p-values
wilcox.test(x = as.numeric(unlist(args2 %>% filter(pUTI89 == "Yes") %>% select(`Total ARGs`))), 
            y = as.numeric(unlist(args2 %>% filter(pUTI89 == "No") %>% select(`Total ARGs`))),
            conf.int = TRUE)

# Visualise plot with p-values
pUTI89_argplot_p <- pUTI89_argplot + stat_compare_means(method = "wilcox.test", label.x = 1.25)    

# IntI1 vs ARGs
# Subset intI1 hits
int <- args %>% select(Name, intI1)

#Join back to arg df
args3 <- left_join(args2, int)

# Convert 1/2s and 0s to Yes/No 
args3 <- args3 %>% mutate(intI1 = gsub(1, "Yes", intI1), intI1 = gsub(0, "No", intI1), intI1 = gsub(2, "Yes", intI1))

# Plot ARG counts by intI1 presence/absence
intI1_argplot <-  ggboxplot(args3, x = "intI1", y = "Total ARGs",
                            color = "intI1", palette = c("#48e0c7", "#f52020"),
                            order = c("Yes", "No"),
                            ylab = "ARGs", xlab = "intI1") +
  scale_x_discrete(label = c("Yes", "No"))

# Check if there is a difference in ARGs by intI1 status
kruskal.test(`Total ARGs` ~ intI1, data = args3)

# Wilcoxon test to quantify the difference with p-values
wilcox.test(x = as.numeric(unlist(args3 %>% filter(intI1 == "Yes") %>% select(`Total ARGs`))), 
            y = as.numeric(unlist(args3 %>% filter(intI1 == "No") %>% select(`Total ARGs`))),
            conf.int = TRUE)

# Plot with p value
intI1_argplot_p <- intI1_argplot + stat_compare_means(method = "wilcox.test", label.x = 1.25) 

# Combine the two plots for figure 5
figure5 <- ggarrange(pUTI89_argplot_p, intI1_argplot_p, 
                     ncol = 2, nrow = 1, labels =c("a)","b)"), legend = "none")

# Save Figure 5
ggsave("Figure5_pUTI89vsARGs.png",
       figure5,
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 110, 
       unit ="mm", 
       dpi = 300)

####-----FIGURE S1 - SOURCES, CLUSTERS, pUTI89-----####
# Some additional figures showing how sources, clusters and pUTI89 carriage relate to one another
# Plot BAP clusters by source
cluster_source <- ggplot(metadata, aes(factor(Cluster, levels = c(paste0(rep("BAP",13), seq(1:13))), 
                                              labels = c(paste0(rep("BAP",13), seq(1:13)))))) + 
  geom_bar(aes(fill = Source)) +
  scale_fill_manual(values = source_clrs) +
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 175), n.breaks = 8) +
  scale_x_discrete(name = "Cluster")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 8))

# Plot BAP clusters by pUTI89 carriage
cluster_pUTI89 <- ggplot(metadata, aes(factor(Cluster, levels = c(paste0(rep("BAP",13), seq(1:13))), 
                                              labels = c(paste0(rep("BAP",13), seq(1:13)))))) + 
  geom_bar(aes(fill = pUTI89)) +
  scale_fill_manual(name = "pUTI89", values = c("Yes" = "#48e0c7", "No" = "#f52020")) +
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 175), n.breaks = 8) +
  scale_x_discrete(name = "Cluster")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = 1, face = "plain", angle = 45),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 8))

# Plot pUTI89 carriage by source
pUTI89_source <- ggplot(metadata, aes(fct_rev(fct_infreq(factor(pUTI89))))) +
  geom_bar(aes(fill = `Source`)) +
  scale_fill_manual(values = source_clrs) +
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 200)) +
  scale_x_discrete(name = "pUTI89", labels = c("Yes", "No"))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 8)) 

# Plot sources by pUTI89 carriage
source_pUTI89 <- ggplot(metadata, aes(Source)) +
  geom_bar(aes(fill = pUTI89)) +
  scale_fill_manual(name = "pUTI89", values = c("Yes" = "#48e0c7","No" = "#f52020")) +
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 175)) +
  scale_x_discrete(name = "Source", labels = function(x) str_wrap(x, width = 10))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 8, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 8)) 

# Plot them all together
figureS1 <- ggarrange(cluster_source, cluster_pUTI89, pUTI89_source, source_pUTI89,
                     ncol = 2, nrow = 2, 
                     labels =c("a)","b)", "c)", "d)"), 
                     hjust = c(-.05, -.05, -0.5, -0.5),
                     legend = "right",
                     align = "h")

# Save Figure 3
ggsave("FigureS1_clusters.png",
       figureS1, 
       path = "outputs/figures/", 
       device = "png", 
       width= 310, 
       height = 210, 
       unit ="mm", 
       dpi = 300)

####-----FIGURE S2-4 GENE HEATMAPS-----####
# Alignment of gene presence/absence with tree and metadata
# Ggtree object for heatmap visualisation
prgo.tree.heatmap <- ggtree(tree, branch.length = "none") %<+% metadata

## AMR HEATMAP ##
figureS2 <- gheatmap(prgo.tree.heatmap, 
         args_heat, 
         width = 30,
         font.size = 2,
         colnames_offset_x = 0.00001,
         colnames_offset_y = 2.2,
         colnames_position = "top",
         colnames_angle = 45,
         hjust = 0,
         color = rgb(0, 0, 0, alpha = .1)) +
  scale_fill_manual(name = "Data", values = c(heat_clrs, bap_clrs, source_clrs)) +
  ylim(c(0,310)) +
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

## VIRULENCE HEATMAP ##
figureS3 <- gheatmap(prgo.tree.heatmap, 
         vags_heat, 
         width = 30,
         font.size = 1,
         colnames_offset_x = 0.00001,
         colnames_offset_y = 2.2,
         colnames_position = "top",
         colnames_angle = 45,
         hjust = 0,
         color = rgb(0, 0, 0, alpha = .2)) +
  scale_fill_manual(name = "Data", values = c(heat_clrs, bap_clrs, source_clrs)) +
  ylim(c(0,310)) +
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

## PLASMID REPLICON HEATMAP ##
figureS4 <- gheatmap(prgo.tree.heatmap, 
         plas_heat, 
         width = 30,
         font.size = 2,
         colnames_offset_x = 0.00001,
         colnames_offset_y = 2.2,
         colnames_position = "top",
         colnames_angle = 45,
         hjust = 0,
         color = rgb(0, 0, 0, alpha = .2)) +
  scale_fill_manual(name = "Data", values = c(heat_clrs, bap_clrs, source_clrs)) +
  ylim(c(0,310)) +
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

# Save Figure S2
ggsave("FigureS2_ARG_heatmap.png",
       figureS2, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 300)

# Save Figure S3
ggsave("FigureS3_VAG_heatmap.png",
       figureS3, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 300)

# Save Figure S4
ggsave("FigureS4_plas_heatmap.png",
       figureS4, 
       path = "outputs/figures/", 
       device = "png", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 300)

####-----SCOARY ANALYSIS-----####
# Generate file list of Scoary results for BAP groups
filelist <- list.files(path="scoary/prokka", pattern="^BAP_.*csv", full.names = TRUE)

# Read in thefiles that actually have data - used file size as a proxy for this
for (f in 1:length(filelist)){
  if(file.size(filelist[f]) > 229){
    assign(paste0(str_extract(filelist[f], "BAP_\\d"), "_scoary"), read_csv(filelist[f]))
  }
}

# Combine into list
bap_scoary_list <- mget(ls(pattern = "BAP_.*_scoary"))

# Processing loop for each file
for (f in 1:length(bap_scoary_list)){
  data <- bap_scoary_list[[f]]
  name <- str_extract(names(bap_scoary_list[f]), "BAP_\\d")
  # Filter hypotheticals and combine gene and non-unique gene name to create a unique name that splits paralogs
  data <- data %>% 
    mutate(BAP = rep(name, nrow(data)), 
           Gene_Unique = paste(`Non-unique Gene name`, Gene, sep = "_")) %>% 
    filter(!grepl("hypothetical", Annotation), Benjamini_H_p < 1E-30) %>% 
    select(Gene_Unique, everything())
  
  # Get rid of NAs
  data$Gene_Unique <- gsub("NA_", "", data$Gene_Unique)
  
  # Split into over and under-represented ColV clade genes based on the Odds Ratio
  over_rep <- data %>% filter(Odds_ratio > 1 | Odds_ratio == "inf") %>%
    # slice(0:20) %>%
    dplyr::select(-Sensitivity,-Specificity,-Bonferroni_p) %>%
    dplyr::rename(
      Pos_present = Number_pos_present_in,
      Neg_present = Number_neg_present_in,
      Pos_absent = Number_pos_not_present_in,
      Neg_absent = Number_neg_not_present_in
    ) 
  
  under_rep <- data %>% filter(Odds_ratio < 1) %>%
    # slice(0:20) %>%
    dplyr::select(-Sensitivity,-Specificity,-Bonferroni_p) %>%
    dplyr::rename(
      Pos_present = Number_pos_present_in,
      Neg_present = Number_neg_present_in,
      Pos_absent = Number_pos_not_present_in,
      Neg_absent = Number_neg_not_present_in
    ) 
  # Assign names to the outputs
  assign(paste(str_extract(names(bap_scoary_list[f]), "BAP_\\d"), "over", sep = "_"), over_rep)
  assign(paste(str_extract(names(bap_scoary_list[f]), "BAP_\\d"), "under", sep = "_"), under_rep)
  
  # Remove unnecessary objects
  rm(data)
  rm(over_rep)
  rm(under_rep)
  rm(name)
}

# Get a list of the data frames from the above loop
bap_output_list <- mget(ls(pattern = "BAP_.*_[over|under]"))

# Create empty list for all significant BAP genes
all_bap_list <- list()

# Process all significant BAP genes into one dataframe
for (f in 1:length(bap_output_list)){
  if(nrow(bap_output_list[[f]]) != 0){
    all_bap_list <- append(all_bap_list, bap_output_list[f])
    all_bap <- bind_rows(all_bap_list)
  }
}

# Create empty list for overrepresented BAP genes
over_bap_list <- list()

# Process overrepresented BAP genes into one dataframe
for (f in 1:length(bap_output_list)){
  if(grepl("over", names(bap_output_list[f])) == TRUE && nrow(bap_output_list[[f]]) != 0){
    over_bap_list <- append(over_bap_list, bap_output_list[f])
    over_bap <- bind_rows(over_bap_list)
  }
}

# Create empty list for underrepresented BAP genes
under_bap_list <- list()

# Process underrepresented BAP genes into one dataframe
for (f in 1:length(bap_output_list)){
  if(grepl("under", names(bap_output_list[f])) == TRUE && nrow(bap_output_list[[f]]) != 0){
    under_bap_list <- append(under_bap_list, bap_output_list[f])
    under_bap <- bind_rows(under_bap_list)
  }
}

# Write out the BAP4 results - this is the only group with appreciably over-represented genes worth talking about
write_csv(BAP_4_over, "outputs/data/ST127.BAP4.scoary.over.csv")

# Bring in the Roary data to confirm that over-represented genes are mostly alternative alleles of genes that are carried by all ST127
roary_raw <- read_delim("../../data/roary_out_MVC107/gene_presence_absence.Rtab") %>% mutate(Gene = gsub("_.*", "", Gene))

# Simplify gene names
roary_raw <- roary_raw %>% mutate(Gene = gsub("_.*", "", Gene))

# Select the BAP4 over-represented genes
roary_raw <-  roary_raw %>% filter(grepl("lpxD|flhB|cheZ|galE|kpsT|kpsM|glcA|glcB|ybhI|dctD|sucB|dctD|sucC",Gene))

# Summary of counts for each gene in the whole collection indicates most of them (but not all) are specific alleles of core genes
BAP4_gene_totals <- roary_raw %>%
  rowwise() %>%
  mutate(
    total = sum(c_across(-Gene))
  ) %>% select(Gene, total)

BAP4_gene_totals %>% group_by(Gene) %>% summarise(sum(total))

####-----TABLES S1, S2-----####
## Process gene screening dataframes to generate Table S1
S1_args <- args %>% select(Name, where(is.integer))
S1_vags <- vags %>% select(Name, where(is.integer))
S1_plas <- plas %>% select(Name, where(is.integer))
S1_mges <- mges %>% select(Name, where(is.integer))

# Join them together
tableS1 <- left_join(S1_args, S1_vags)
tableS1 <- left_join(tableS1, S1_plas)
tableS1 <- left_join(tableS1, S1_mges)
tableS1 <- left_join(all_meta, tableS1)

# Write to file
write_csv(tableS1, "outputs/data/TableS1_ST127.csv")

## Process gene co-occurence dataframe to show contigs with intI1 and ARGs for Table S2
tableS2 <- ST127_co_occurence_N90L95 %>% 
  filter(grepl("intI1", same_scaff), name %in% metadata$Name) %>% 
  mutate(same_scaff = gsub("card_", "", same_scaff), 
         same_scaff = gsub("EC_custom_", "", same_scaff),
         same_scaff = gsub("Shigella_flexneri_chloramphenicol_acetyltransferase", "catA1", same_scaff),
         same_scaff = gsub("ISfinder_Feb_2020_", "", same_scaff),
         same_scaff = gsub("_HQ730118.1", "", same_scaff),
         same_scaff = gsub("_CP022086.1","", same_scaff),
         same_scaff = gsub("plasmidfinder_", "", same_scaff),
         same_scaff = gsub(":[^[:space:]]*","", same_scaff)) %>%
  rename(Name = name, `Assembly Contig` = SEQUENCE, Genes = same_scaff)

# Write to file
write_csv(tableS2, "outputs/data/TableS2_ST127.csv")