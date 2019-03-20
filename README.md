# rhto-GEM: Genome-scale metabolic model of _Rhodotorula toruloides_

[![DOI](https://zenodo.org/badge/133515828.svg)](https://zenodo.org/badge/latestdoi/133515828) [![GitHub version](https://badge.fury.io/gh/sysbiochalmers%2Frhto-gem.svg)](https://badge.fury.io/gh/sysbiochalmers%2Frhto-gem) 

- Brief model description

This repository contains the current genome-scale metabolic model of _Rhodotorula toruloides_, named **rhto-GEM**. For the latest updated release see [here](https://github.com/SysBioChalmers/rhto-GEM/releases).

- Abstract

_Rhodosporidium toruloides_ (syn. _Rhodotorula toruloides_) is a basidiomycetous yeast belonging to the subphylum Pucciniomycotina and occurs naturally in a wide range of habitats including surfaces of leaves, soil and sea water. The broad substrate range of _R. toruloides_ and its ability to accumulate lipids exceeding half of its cell dry weight has made this yeast a popular system for production of biological oils from inedible substrates (e.g. pentose sugars, crude glycerol). The _R. toruloides_ lipid fraction contains omega-3 linolenic acid and heptadecenoic acid, which makes this yeast a promising organism for production of pharma- and nutraceuticals. _R. toruloides_ also produces a number of carotenoid pigments including torularhodin, torulene, γ-carotene and β-carotene. Other applications of _R. toruloides_ include the production of the enzymes L-phenylalanine ammonia-lyase and D-amino acid oxidase.

- Model keywords

**GEM category:** Species; **Utilisation:** experimental data reconstruction; **Field:** metabolic-network reconstruction; **Type of model:** reconstruction, curated; **Model source:** [yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM) & [iYali](https://github.com/SysBioChalmers/Yarrowia_lipolytica_W29-GEM); **Omic source:** genomics; **Taxonomy:** _Rhodotorula toruloides_; **Metabolic system:** general metabolism; **Bioreactor**; **Strain:** NP11; **Condition:** minimal medium;

- Reference:  
> Tiukova IA _et al_. (2019) "Genome-scale model of _Rhodotorula toruloides_ metabolism" bioRxiv doi:[10.1101/528489](https://doi.org/10.1101/528489)

- Last update: 2019-03-20

- Main model descriptors:

| Taxonomy | Template Model | Reactions | Metabolites | Genes |
| ------------- |:-------------:|:-------------:|:-------------:|:-----:|
| _Rhodotorula toruloides_|	[yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM) & [iYali](https://github.com/SysBioChalmers/Yarrowia_lipolytica_W29-GEM) | 4930 | 3374 | 926 |

A [Memote](https://memote.readthedocs.io/en/latest/) snapshot report of the most recent release is available [here](https://SysBioChalmers.github.io/rhto-GEM).

This repository is administered by [@edkerk](https://github.com/edkerk/), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.

## Citation

* If you use rhto-GEM in your research, whether for simulations, data analysis, model reconstruction of other purposes, we ask you to cite the rhto-GEM paper on [bioRxiv](https://doi.org/10.1101/528489).
* In addition, it is good practice to cite the specific version of rhto-GEM that you used, to improve reproducibility. All rhto-GEM releases are archived in [Zenodo](https://zenodo.org/badge/latestdoi/133515828). Find the corresponding DOI for each release [here](https://zenodo.org/search?page=1&size=20&q=conceptrecid:2547988&sort=-publication_date&all_versions=True).

## Installation

### Required software

  * This model is recommended to be used with the [**RAVEN toolbox for MATLAB**](https://github.com/SysBioChalmers/RAVEN) (version 2.0).
  * Alternatively, the model can also be directly used with the [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox), [cobrapy](https://github.com/opencobra/cobrapy), or other software or toolboxes that support genome-scale models in SBML L3V1 FBCv2 format.

### Dependencies
* Please see the [RAVEN toolbox](https://github.com/SysBioChalmers/RAVEN) repository for dependencies regarding RAVEN.
* For contribution to development: a [git wrapper](https://github.com/manur/MATLAB-git) added to the search path.

### Installation instructions
* Just want to use the model? Clone it from [`master`](https://github.com/SysBioChalmers/rhto-GEM) in the Github repo, or just download [the latest release](https://github.com/SysBioChalmers/rhto-GEM/releases).
* Wish to also contribute? Fork it to your Github account, and create a new branch from [`devel`](https://github.com/SysBioChalmers/rhto-GEM/tree/devel).

## Model files

The model is available in `.xml`, `.txt`, `.yml`, `.mat` and `.xlsx` (the last 2 extensions only in `master`). Additionally, versions of toolboxes & SBML used for saving the model are tracked in `dependencies.txt`.

### Complementary scripts

* `newCommit.m`: prepares files from a modified model for a new GitHub commit in a development branch.
* `newCommit.m`: prepares files from a modified model for a new GitHub release in the master branch.
* `analysis`: folder with scripts performing analyses on rhto-GEM.
* `curation`: folder with additional curation scripts, after the initial reconstruction
* `experimental`: folder with scripts that modify _rhto-GEM_ to incorporate experimental data.
* `reconstruction`: folder with scripts that detail the model reconstruction and curation process.
* `validation`: folder with scripts validating performance of rhto-GEM.

### Complementary data

* `data`: experimentally measured data.
* `genome`: protein fasta files of _R. toruloides_, _S. cerevisiae_ and _Y. lipolytica_, as used for identifying orthologs.
* `meneco`: input and output files for meneco, as run for `r4_gapfilling`.
* `reconstruction`: miscellaneous files used during model reconstruction.
* `validation`: experimental data used to validate _rhto-GEM_.

## Contributors
* [Ievgeniia Tiukova](https://www.chalmers.se/en/staff/Pages/tiukova.aspx) ([@itaova](https://github.com/itaova)), Chalmers University of Technology, Gothenburg Sweden
* [Eduard J. Kerkhoven](https://www.chalmers.se/en/staff/Pages/Eduard-Kerkhoven.aspx) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Gothenburg Sweden
