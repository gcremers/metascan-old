# metascan
Metabolic scanning and annotation of Metagenomes

Metascan is an Metagenomic Scanning and Annotation tool, with an emphasis on metabolic genes.
The heart of Metascan is a set of metabolic core genes, that are used to paint a picture of the metabolic capacity of the sample.
Furthermore, it utilizes the Kegg pathways for a complete metabolic overview of each sample.

Samples can be analyzed as eiter binned or unbinned metagenome.

Metascan consists of a perl script, a few auxillary textfiles and a large set of HMM profiles, created by clustering TrEmbl proteins, based on Kegg K-numbers.
Since it is a Prokka adaptation, it can therefor be run on any system which can run Prokka, just by downloading the script and the databases, without the use of inception-software. The only modification that needs to be done is to direct the script to the right location of the database.

