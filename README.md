# rhto-GEM: Genome-scale metabolic model of _Rhodosporidium toruloides_

- Brief model description

This repository contains the current genome-scale metabolic model of _Rhodosporidium toruloides_, named **rhto-GEM**. For the latest updated release see [here](https://github.com/SysBioChalmers/rhto-GEM/releases).

- Abstract

_Rhodosporidium toruloides_ (syn. _Rhodotorula toruloides_) is a basidiomycetous yeast belonging to the subphylum Pucciniomycotina and occurs naturally in a wide range of habitats including surfaces of leaves, soil and sea water. The broad substrate range of _R. toruloides_ and its ability to accumulate lipids exceeding half of its cell dry weight has made this yeast a popular system for production of biological oils from inedible substrates (e.g. pentose sugars, crude glycerol). The _R. toruloides_ lipid fraction contains omega-3 linolenic acid and heptadecenoic acid, which makes this yeast a promising organism for production of pharma- and nutraceuticals. _R. toruloides_ also produces a number of carotenoid pigments including torularhodin, torulene, γ-carotene and β-carotene. Other applications of _R. toruloides_ include the production of the enzymes L-phenylalanine ammonia-lyase and D-amino acid oxidase.

- Model keywords

**GEM category:** Species; **Utilisation:** experimental data reconstruction; **Field:** metabolic-network reconstruction; **Type of model:** reconstruction, curated; **Model source:** [yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM); **Omic source:** genomics; **Taxonomy:** _Rhodosporidium toruloides_; **Metabolic system:** general metabolism; **Bioreactor**; **Strain:** NP11; **Condition:** minimal medium;

- Reference:  
>bioRxiv

- Last update: 2019-01-22

- Main model descriptors:

| Taxonomy | Template Model | Reactions | Metabolites | Genes |
| ------------- |:-------------:|:-------------:|:-------------:|:-----:|
| _Rhodosporidium toruloides_|	[yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM) | 4869 | 3334 | 896 |


This repository is administered by [@edkerk](https://github.com/edkerk/), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.


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

* `newCommit.m`: RAVEN function that prepares files from a modified model for a new GitHub commit in a development branch.
* `newCommit.m`: RAVEN function that prepares files from a modified model for a new GitHub release in the master branch.
* `experimental`: folder with scripts that modify _rhto-GEM_ to incorporate experimental data
  * `adjustRhtoBiomass.m`: RAVEN function that adjusts the biomass equation to match the specified lipid data.
  * `readTiukovaData.m`: RAVEN function that reads experimental lipid class, lipid chain length and exchange flux data, as provided in ComplementaryData.
  * `scaleAbundancesRhto.m`: RAVEN function as part of the SLIMEr approach, which scales lipid backbones or chains to reconcile experimental measurements. Called by `adjustRhtoBiomass.m`.
  * `simulateRhtoGrowth.m`: RAVEN function that simulates growth according to the provided flux data, called by `scaleAbundancesRhto.m`.
* `reconstruction`: folder with scripts that detail the model reconstruction and curation process.
  * `makeRxns.m`: RAVEN function that takes lipid templates reactions to produce SLIMEr reactions, called by `r5_lipids_curation.m`.
  * `r1_draft_from_homology.m`: RAVEN script that generates a first draft GEM from homology to _S. cerevisiae_.
  * `r2_additional_yeastGEM_reactions.m`: RAVEN script that performs manual curation based on yeast-GEM reactions.
  * `r3_add_rhto_specific_reactions.m`: RAVEN script that adds _R. toruloides_ specific reactions through manual curation.
  * `r4_gapfilling.m`: RAVEN script that details the gapfilling performed to reach a functional model.
  * `r5_lipids_curation.m`: RAVEN scripts that curates lipid metabolism by e.g. adding relevant SLIMEr reactions.
  * `r6_grRules_curation.m`: RAVEN script that corrects gene associations.
  * `menecoGapFill.sh`: BASH script that specifies meneco commands run to identify reactions for gapfilling, as used in `r4_gapfilling.m`.
  * `lipidTemplates.txt` and `lipidTransport.txt`: template reactions for lipid metabolism, as called by `r5_lipids_curation.m`.

### Complementary data

* `data`: experimentally measured data.
* `genome`: protein fasta files of _R. toruloides_ and _S. cerevisiae_, as used for identifying orthologs.
* `reconstruction`: miscellaneous files used during model reconstruction. 


## Contributors
* [Ievgeniia Tiukova](https://www.chalmers.se/en/staff/Pages/tiukova.aspx) ([@itaova](https://github.com/itaova)), Chalmers University of Technology, Gothenburg Sweden
* [Eduard J. Kerkhoven](https://www.chalmers.se/en/staff/Pages/Eduard-Kerkhoven.aspx) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Gothenburg Sweden
