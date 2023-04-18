# Eukaryotpying-Python
Python implementation of Plucinski's Bayesian algorithm and Barratt's heuristic definition of genetic distance

The two methods described here compute genetic distances for subsequent heirarcical clustering. The method requires MLST data as input represented specifically in the form of a Haplotype Data Sheet or HDS. An example of this HDS format is provided in this repository in the "haplotype_sheets" directory. 

These two genetic distance computation methods (Plucinski's Bayesian method and Barratt's heuristic definition of genetic distance) were developed to address four main issues relating to the analysis of massive and complex MLST datasets:

(1) The common occurrence of missing data in genotyping datasets (e.g. such as a situation where 2 or 3 out of 6 multi-locus sequence typing loci fail to amplify for a subset of your specimens).

(2) The use of MLST methods in the context of sexually reproducing populations, where even closely related individuals may not possess the same genotype (and may be heterozygous) due to chromosomal crossover and random reassortment of chromosomes as occurs during meiosis.

(3) The issue of analyzing specimens which may be extremely complex, potentially representing mixed populations of individuals. Ever try to construct a phylogeny or generate a cluster dendrogram in a situation where for one MLST marker you detect one haplotype, at another you detect three, at another you detect four and another you detect two - in the same specimen? This is essentially what we deal with when we attempt to genotype Cyclospora cayetanensis directly from human stool. It gets extremely complicated.

(4) The absence of distance statistics that appropriately consider all aspects of genetic data (e.g. allele frequency, entropy of loci, nuclear versus extranuclear inheritance). Simpler metrics such as Bray-Curtis dissimilarty and Jaccard distances fail to consider these aspects of genetic data.


### Running this code

>While in the folder with all the files from the cloned Eukaryotyping github run:

```bash
bash run_BayHeur.sh
```
> This will analyze 99 samples and produce a pairwise matrix of distances that can be used for downstream analysis.  


## Additional Information

For additional detailed information on how these algorithms work please refer to our project [background](background.md).


## Deployment

<!-- need to update once on SciComp and CDCgov github -->

This section will be updated in the future.


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
