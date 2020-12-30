# Using biomaRt in R

BiomaRt is a great tool for locating and visualizing information on chromosomes and genes, and biomaRt allows us to employ R in our exploration and analysis.

The project uses GenomeGraphs package to visualize the HERC2 gene. 
Here I'm modifying the base code from Jonathan Pevsner visualising the HBB gene.

We will access Ensembl's Biomart RESTful service for information on the HERC2 gene, then render a web page of the information we get back from Biomart as visualized using the R package GenomeGraphs.


References: 
1. herc2 gene details (Ensemble id, chromosome and location of HERC2): https://uswest.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000128731;r=15:28111040-28322172
2. Jonathan Pevsner, Bioinformatics and Functional Genomics  (August, 2015), Chapter 8: DNA: The Eukaryotic Chromosome)



