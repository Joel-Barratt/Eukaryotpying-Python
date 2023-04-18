# Eukaryotpying-Python
This repository contains Python implementations of Plucinski's Bayesian algorithm and Barratt's heuristic definition of genetic distance. These two methods compute genetic distances for subsequent heirarcical clustering. The methods requires MLST data as input represented specifically in the form of a Haplotype Data Sheet or HDS. An example of this HDS format is provided in this repository in the "haplotype_sheets" directory. The code will then compute four genetic distances matrices and print these in the "ensemble_matrices" directory.

These two genetic distance computation methods (Plucinski's Bayesian method and Barratt's heuristic definition of genetic distance - nicknamed the "Eukaryotyping" method) were developed to address four main issues relating to the analysis of massive and complex MLST datasets:

 1. The common occurrence of missing data in genotyping datasets (e.g. such as a situation where 2 or 3 out of 6 multi-locus sequence typing loci fail to amplify for a subset of your specimens).

 2. The use of MLST methods in the context of sexually reproducing populations, where even closely related individuals may not possess the same genotype (and may be heterozygous) due to chromosomal crossover and random reassortment of chromosomes as occurs during meiosis.

 3. The issue of analyzing specimens which may be extremely complex, potentially representing mixed populations of individuals. Ever try to construct a phylogeny or generate a cluster dendrogram in a situation where for one MLST marker you detect one haplotype, at another you detect three, at another you detect four and another you detect two - in the same specimen? This is essentially what we deal with when we attempt to genotype *Cyclospora cayetanensis* directly from human stool. It gets extremely complicated.
 
4. The absence of distance statistics that appropriately consider all aspects of genetic data (e.g. allele frequency, entropy of loci, nuclear versus extranuclear inheritance). Simpler metrics such as Bray-Curtis dissimilarty and Jaccard distances fail to consider these aspects of genetic data.

_Please cite the following manuscripts:_

```
1. Barratt, JLN, S Park, FS Nascimento, J Hofstetter, M Plucinski, S Casillas, RS Bradbury, MJ Arrowood, Y Qvarnstrom, E Talundzic (2019) Genotyping genetically heterogeneous Cyclospora cayetanensis infections to complement epidemiological case linkage. Parasitology:1–9 doi:10.1017/S0031182019000581


2. Nascimento, FS, JLN Barratt, K Houghton, M Plucinski, J Kelley, S Casillas, C Bennett, C Snider, R Tuladhar, J Zhang, B Clemons, S Madison-Antenucci, A Russell, E Cebelinski, J Haan, T Robinson, MJ Arrowood, E Talundzic, RS Bradbury, and Y Qvarnstrom (2020) Evaluation of an ensemble-based distance statistic for clustering MLST datasets using epidemiologically defined clusters of cyclosporiasis. Epidemiology & Infection: 148, e172, 1–10. https://doi.org/10.1017/
S0950268820001697

3. Jacobson, D., Y Zheng, MM Plucinski, Y Qvarnstrom, JLN Barratt (2022) Evaluation of various distance computation methods for construction of haplotype-based phylogenies from large MLST datasets. Molecular Phylogenetics and Evolution: 177, 107608. https://doi.org/10.1016/j.ympev.2022.107608
```


### Running this code

>Update the "euktypDir" in the "run_BayHeur.sh" script. Next, within the "/DISTCOMP/Pycode_distcomp/" directory update the "markerList.txt" file to indicate your marker requirements. Then, while in the folder with all the files from the cloned Eukaryotyping-Python github run:

```bash
bash run_BayHeur.sh
```
> This will analyze the haplotype data sheet (HDS) within the haplotype_sheets directory and produce four pairwise genetic distance matrices that can be used for downstream analysis.
> TIP: Open the "run_BayHeur.sh" script and take a look at the comments in the script. This will provide information on how to modify certain flags to fine-tune the minimum data availability requirements. Dependencies are fairly standard (numpy, pandas etc) and you will likely already have these installed. If you obtain a ModuleNotFoundError when trying to run the code then just install the missing module.



## Acknowledgments

* Sincere thanks to Dr Yueli Zhang who wrote the majority of this code.


## License
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This soruce code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.


## Contributing
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.
