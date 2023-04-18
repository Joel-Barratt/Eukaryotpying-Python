# Eukaryotpying-Python
Python implementation of Plucinski's Bayesian algorithm and Barratt's heuristic definition of genetic distance

The two methods described here compute genetic distances for subsequent heirarcical clustering. The method requires MLST data as input represented specifically in the form of a Haplotype Data Sheet or HDS. An example of this HDS format is provided in this repository in the "haplotype_sheets" directory. 

These two genetic distance computation methods (Plucinski's Bayesian method and Barratt's heuristic definition of genetic distance) were developed to address four main issues relating to the analysis of massive and complex MLST datasets:

(1) The common occurrence of missing data in genotyping datasets (e.g. such as a situation where 2 or 3 out of 6 multi-locus sequence typing loci fail to amplify for a subset of your specimens).

(2) The use of MLST methods in the context of sexually reproducing populations, where even closely related individuals may not possess the same genotype (and may be heterozygous) due to chromosomal crossover and random reassortment of chromosomes as occurs during meiosis.

(3) The issue of analyzing specimens which may be extremely complex, potentially representing mixed populations of individuals. Ever try to construct a phylogeny or generate a cluster dendrogram in a situation where for one MLST marker you detect one haplotype, at another you detect three, at another you detect four and another you detect two - in the same specimen? This is essentially what we deal with when we attempt to genotype Cyclospora cayetanensis directly from human stool. It gets extremely complicated.

(4) The absence of distance statistics that appropriately consider all aspects of genetic data (e.g. allele frequency, entropy of loci, nuclear versus extranuclear inheritance). Simpler metrics such as Bray-Curtis dissimilarty and Jaccard distances fail to consider these aspects of genetic data.
