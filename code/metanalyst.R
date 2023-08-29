library(MetaboAnalystR)

download.file("https://www.metaboanalyst.ca/MetaboAnalyst/resources/data/data1.csv", "data1.csv")
download.file("https://www.metaboanalyst.ca/MetaboAnalyst/resources/data/data2.csv", "data2.csv")
download.file("https://www.metaboanalyst.ca/MetaboAnalyst/resources/data/data3.csv", "data3.csv")
download.file("https://www.metaboanalyst.ca/MetaboAnalyst/resources/data/data4.csv", "data4.csv")
# Create objects for storing processed data from meta-analysis
mSet <- InitDataObjects("conc", "metadata", FALSE)
# Read in example data: adenocarcinoma data2
mSet <- ReadIndData(mSet, "data1.csv", "colu");
# Sanity check data to ensure it is ready for analysis
mSet <- SanityCheckIndData(mSet, "data1.csv")
## to view any messages created during the sanity check
mSet$dataSet$check.msg
# Perform log-transformation
mSet <- PerformIndNormalization(mSet, "data1.csv", "log", 1);
#Perform differential expression analysis to identify DE features
mSet <- PerformLimmaDE(mSet, "data1.csv", 0.05, 0.0);
# Repeat steps for example data3
mSet <- ReadIndData(mSet, "data3.csv", "colu");
mSet <- SanityCheckIndData(mSet, "data3.csv")
mSet <- PerformIndNormalization(mSet, "data3.csv", "log", 1);
mSet <- PerformLimmaDE(mSet, "data3.csv", 0.05, 0.0);
# Repeat steps for example data4
mSet <- ReadIndData(mSet, "data4.csv", "colu");
mSet <- SanityCheckIndData(mSet, "data4.csv")
mSet <- PerformIndNormalization(mSet, "data4.csv", "log", 1);
mSet <- PerformLimmaDE(mSet, "data4.csv", 0.05, 0.0);
# Check if meta-data between all uploaded datasets are consistent
mSet <- CheckMetaDataConsistency(mSet, F);
###*** Choose one of 3 methods to perform meta-analysis *** ###
###*** OPTION 1 - COMBINE P-VALUES *** ###
mSet <- PerformPvalCombination(mSet, "fisher", 0.05)
###*** OPTION 2 - PERFORM VOTE COUNTING *** ###
mSet <- PerformVoteCounting(mSet, 0.05, 2.0)
###*** OPTION 3 - MERGE INTO MEGA-DATASET *** ###
mSet <- PerformMetaMerge(mSet, 0.05)
# Create results table
mSet <- GetMetaResultMatrix(mSet, "fc")
mSet <- PlotSelectedFeature(mSet, "pyrophosphate")
# Prepare data for the Venn Diagram, which will create a Integrated Venn diagram in your working directory (two overlapping circles, highlighting novel biomarker features from the meta-analysis, biomarkers that were consistently identified using meta-analysis and individual DE expression, and biomarkers that were only identified using individual DE expression.)
mSet <- PrepareVennData(mSet);
# Explore the Venn Diagram in the "vennData" object created
# Get names of features overlapping between selected datasets from "vennData"
mSet <- GetVennGeneNames(mSet, "data1.csvdata3.csvmeta_dat");
# Enter the object below to obtain the names of the features that overlap between all of the studies
mSet$dataSet$venn_overlap



#enrichment analysis
## When input is a list
# Create vector consisting of compounds for enrichment analysis
tmp.vec <- c("Acetoacetic acid", "Beta-Alanine", "Creatine", "Dimethylglycine", "Fumaric acid", "Glycine", "Homocysteine", "L-Cysteine", "L-Isolucine", "L-Phenylalanine", "L-Serine", "L-Threonine", "L-Tyrosine", "L-Valine", "Phenylpyruvic acid", "Propionic acid", "Pyruvic acid", "Sarcosine")
# Create mSetObj
mSet<-InitDataObjects("conc", "msetora", FALSE)
#Set up mSetObj with the list of compounds
mSet<-Setup.MapData(mSet, tmp.vec);
# Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
mSet<-CrossReferencing(mSet, "name")
mSet$name.map
