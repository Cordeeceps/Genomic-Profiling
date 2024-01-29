###########################################################################
# Group Members: Meghan Bothoff-Shanahan, Justin Vance, 
#                and Julia Christensen
###########################################################################
# Part 1

# Children's Brain Tumor Tissue Consortium (CBTTC)

# About the Project: We are identifying differentially expressed genes between 
#                    Caucasian individuals and African American individuals in 
#                    child brain tumor tissue samples. 

###########################################################################

library(dplyr)
library(ggplot2)
library(UCSCXenaTools) # needed to retrieve data
library(edgeR) # needed for processing, such as TMM
library(limma) # needed to find DE probes
library(class) # for knn

###########################################################################

# Step 2: 

# Retrieve the expression and phenotype data

###########################################################################

data(XenaData)

# limit to desired cohort
cbttc <- XenaData %>% filter(XenaCohorts == 'Pediatric Brain Tumor Atlas: CBTTC')

# Get the clinical data
cli_query = cbttc %>%
  filter(Label == "Participants information") %>%  # select clinical dataset
  XenaGenerate() %>%  # generate a XenaHub object
  XenaQuery() %>%     # generate the query
  XenaDownload()      # download the data


# prepare (load) the data into R
cbttc_part = XenaPrepare(cli_query)

# Get the RSEM expected_count
# instead of HTSeq - Counts in our case
cli_query <- cbttc %>% filter(Label == 'RSEM expected_count') %>%
  XenaGenerate() %>%  # generate a XenaHub object
  XenaQuery() %>%
  XenaDownload(download_probeMap = TRUE)

# prepare (load) the data into R
cbttc_counts <- XenaPrepare(cli_query)


###########################################################################

# Step 3: 

# Pre-process the data as appropriate (remove unwanted samples, make sure the
# order of samples in the clinical data matches the order of samples in the
# expression data).

###########################################################################

X <- data.frame(cbttc_counts$CBTTC_rsem.genes_expected_count.txt.gz)
rownames(X) <- X$gene_id
X <- X[,-1] # this removes the probes names

# probeMap = probe names
probeMap <- cbttc_counts$cavaticaKidsFirst_gencode.v27.primary_assembly.annotation.gene.probemap

# participants info
Y <- cbttc_part

# remove dots and replace with dashes for formatting consistancy
colnames(X) <- gsub('\\.', '-', colnames(X))

# filter out unknown races
Y <- filter(Y, Race != 'Reported Unknown')
Y <- filter(Y, Race != "More Than One Race" )

# match expression data to clinical data
common_samples <- intersect(colnames(X), Y$xena_sample)
mx <- match(common_samples, colnames(X))
my <- match(common_samples, Y$xena_sample)

X <- X[,mx]
Y <- Y[my,]

# Make sure that the samples match -- if they don't, this will produce an error
stopifnot(all(colnames(X) == Y$xena_sample))


###########################################################################

# Step 4:

# Process the expression data by removing probes with low expression and using
# TMM normalization. If your data is on the scale x = log2 (count + 1), then you will
# need to get the count data back by calculating 2x – 1.

###########################################################################

X <- round(2**X - 1)

# Removing genes with lower counts
dge <- DGEList(counts = X)
keep <- filterByExpr(dge,  min.prop = .10)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# Applying normalization
dge <- calcNormFactors(dge, method = 'TMM')

# Calculating logCPM values
logCPM <- cpm(dge, log=TRUE, prior.count = 3)

###########################################################################

# Step 5: 

# Generate a boxplot of the first 10 samples of your data, using
# boxplot(logCPM[,1:10]), to show that your data have been normalized. Your
# boxplot should have a meaningful title (using the main argument) and y-axis
# label (using the ylab argument).

###########################################################################

boxplot(logCPM[,1:10], main='CBTTC Brain Tumor', 
        ylab = "Log 2 Expression or Gene Expression")


###########################################################################

# Step 6: 

# How many samples were profiled, and how many probes are there in the
# dataset?

###########################################################################

# the total number of probes 
nrow(X)

# the total number of samples
ncol(X)


###########################################################################

# Step 7: Number of samples in each group

# Extract the column from the phenotype data that contains the categories that
# you would like to compare (e.g., the gender column), and use R to output the
# number of samples in each group. Note: in some cases, you may have to process
# the data first, for example if the values were Male1, Male2, Female1, Male3,
# Female2, etc, then the number would need to be removed. 

###########################################################################

table(Y$Race)

###########################################################################

# Step 8: 

# Using limma, find the probes that are differentially expressed across the 
# groups you are comparing, using a false discovery rate (FDR) of 10%. Output 
# the number of probes identified, using the nrow function. Note: it is possible
# for this number to be 0. For some datasets, particularly when the sample size
# is small, you may end up with an FDR of 100% for all probes! If the number of
# probes found using an FDR of 10% is less than 30, then find the top 30 
# probes. (Note: only do this last part if the number of probes is less than 30).

###########################################################################

# construct design matrix 
race <- Y$Race
table(race)

# change column names
design <- model.matrix(~-1 + race)
colnames(design) <- c("Alaskan_Native", "Asian", "African_American", 
                      "Hawaiian_Native", "Caucasian")
head(design)

# fit the linear model to each row of the expression matrix
fit <- lmFit(logCPM, design)
con.matrix <- makeContrasts(Caucasian - African_American, levels = design)
fit <- contrasts.fit(fit, con.matrix)

# calculate moderated t-statistics
fit.de <- eBayes(fit, trend = TRUE)
tt.10 <- topTable(fit.de, sort.by = "p", p.value = 0.10, number = Inf)

#output the # of differentially expressed probes across all groups
nrow(tt.10)

###########################################################################

# Step 9: 

# For the top probe (with the lowest adjusted p-value), construct a boxplot 
# (using ggplot) to compare the expression of that probe across groups. The 
# title of the boxplot should consist of the fold change (FC) and the FDR for 
# the probe, and the boxplot must be constructed using ggplot. Note that the 
# title includes the FC and not the log FC. If you were unable to complete the 
# previous question, you can use the first probe (1st row of logCPM) for this 
# step.

###########################################################################

probe <- rownames(tt.10)[1]
m <- match(probe, rownames(logCPM))
df <- data.frame(expr=logCPM[m,], race = race)

logFC <- tt.10$logFC[1]
FC <- paste0("FC = ", round(2**logFC, 2))
main <- paste0("Expression of Top Probe ", probe, ", ", FC)

df <- filter(df, race == c('Black or African American', 'White'))

ggplot(df, aes(x = race, y = expr, fill = race)) + geom_boxplot() +
  ylab("log2 expression") + ggtitle(main) + 
  scale_fill_manual(values = c("yellow", "darkblue")) +
  theme_classic() + theme(legend.position = "none")


###########################################################################

# Step 10: 

# Find the gene names corresponding to all probes from question #8 
# (or 9 if this was completed) and create a data frame with the following 
# columns: gene names, probe names, logFC, adjusted p-values. Your data 
# frame should not contain any other columns. Output the first 5 rows to 
# display the top 5 probes.

###########################################################################

# creates data frame containing probe name, logFC, and adjusted p value
df2 <- data.frame(ProbeName=(rownames(tt.10)), logFC = tt.10$logFC,
                  adj_pValue = tt.10$adj.P.Val)

#matches probes from df2 with probes in probeMap
m.3 <- match(df2$ProbeName, probeMap$id)

#adds gene names to df2
df2 <- mutate(df2, geneName = probeMap$gene[m.3])

# displays the first 5 rows
df2[1:5,]

###########################################################################

# Step 11: 

# Generate a heatmap showing the expression of the top 30 genes across samples.
# Use the gene names to label the rows of the heatmap, rather than the probe
# names. Color code the samples based on the two groups that you have analyzed.

###########################################################################

gene <- logCPM[1:30,]
probes <- rownames(gene)
m3 <- match(probes, rownames(logCPM))
expr <- logCPM[m3,]
m3 <- match(probes, probeMap$id)
genes <- paste0(probeMap$gene[m3], '(', probeMap$chrom[m3], ')')
rownames(expr) <- genes

col.heat <- colorRampPalette(c('orange', 'darkblue'))(200)

col.race <- as.integer(as.factor(race))
col.race <- c('lightgreen', 'magenta')[col.race]

heatmap(expr, ColSideColors = col.race, col = col.heat, scale = 'none')

###########################################################################

# Step 12: 

# Using DAVID, identify Gene Ontology (GO) terms and KEGG pathways that are
# associated with the differentially expressed genes identified in (11). 
# A file containing the gene names and a screenshot of these results should 
# be submitted with your R notebook.

###########################################################################

m4 <- match(rownames(gene), probeMap$id)
genes <- probeMap$gene[m4]
genes <- unique(genes)

write.table(genes, row.names = FALSE, quote = FALSE)

###########################################################################

# Step 13: 

# Summarize your results based on your analysis. Your summary should include
# the name of the dataset you’ve analyzed, with a link, a description of the 
# samples or individuals, including the number of samples and number of probes 
# that were profiled, the number of samples in each group, the number of 
# differentially expressed probes with an FDR of 10% or less, the names of the 
# top 3 genes, and the top GO terms or pathways associated with the phenotype 
# that you analyzed.

###########################################################################

# Dataset Name: 
# Children's Brain Tumor Tissue Consortium (CBTTC)

# Link: 
# https://xenabrowser.net/datapages/?cohort=Pediatric%20Brain%20Tumor%20Atlas%3A%20CBTTC&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

# Description: 
# This data set includes tumor tissue samples from children with pediatric brain cancer.
# We have filtered out duplicate data as well as data from races that are not Caucasian 
# or African American. 

# Total Number of Probes (after filtering): 
# 58,347

# Total Number of Samples (after filtering): 
# 796

# Total Number of Caucasian Samples: 
# 682

# Total Number of African American Samples: 
# 75

# Number of differentially expressed probes (FDR of <= 10%): 
# 2172 

# Names of the Top 3 Genes: 
# AL008721.2, TBC1D3L, PWP2

# Top GO Terms or Pathways associated with our data set: 
# cholesterol biosynthetic process
# positive regulation of insulin secretion involved in ceullar response to glucose stimulus 
# negative regulation of endothelian cell proliferation
# cell redox homeostasis 
# protein-DNA complex
# extracellular exosome
# protein binding 
