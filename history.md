# History

### rhto v1.2.1:
* Features:
  * Curate grRules and existence of several reactions based on [doi:10.1186/s12934-015-0217-5](https://doi.org/10.1186/s12934-015-0217-5) (PR #34).
  * Curate grRules of xylose and glucose transporters based on differential expression in proteomics study [doi:10.1186/s13068-019-1478-8](https://doi.org/10.1186/s13068-019-1478-8) (PR #34).
* Fix:
  * Bug in `sumBioMass` resulted in incorrect fitted carbohydrate fraction of biomass equation, also resulting in incorrect NGAM (PR #33).
  * Remove four duplicate reactions (`r_4235`, `y200008`, `r_1000`, `r_4262`) (PR #34).
* Refactor:
  * Standardize two complex gene associations (`r_0831`, `r_0832`) (PR #34).

### rhto v1.2.0:
* Features
  * Set alternative carbon source exchange reactions as reversible (glycerol, xylose) (PR #28).
  * `init_rhto-GEM.m` script that prepares MATLAB environment with path to repo's folders (PR #28).
  * `makeLipidRxns` now takes list of acyl-chain configurations, instead of defining these by itself from a list of allowed acyl-chains.  (PR #28).
  * Correct gene associations of mitochondrial beta-oxidation reaction (PR #29).
  * Curate subsystems and EC numbers (PR #29).
  * Update DNA, RNA and amino acid pseudoreactions, based on frequency of nucleotides and amino acids in genome, transcript and protein FASTA files (PR #30).
  * Scale carbohydrate fraction to ensure that biomass pseudoreaction generates 1 gDCW (PR #30).
  * Refit energetics (PR #30).
* Fix
  * Curate acyl-chains: simplify number of possible chain combinations, based on measured chain distribution (e.g. 16:1 level is around 1% and therefore not included in the biomass), and distribution of acyl-chains on TAGs ([10.1007/s00253-017-8126-7](http://doi.org/10.1007/s00253-017-8126-7)):  sn-1 is 16:0 or 18:0, sn-2 is 18:1, 18:2 or 18:3; sn-3 is 16:0, 18:0 or 18:1. (PR #28).
  * Use template SLIMER reactions (PR #28).
  * Various curations of lipid metabolism (PR #28, PR #30).

### rhto v1.1.1:
* Features
  * Make Memote snapshot report via Travis with every push to `master`, report available at https://sysbiochalmers.github.io/rhto-GEM/ (PR #23).
  * Sort reactions and metabolites when saving model, for easier `diff` (PR #23).
  * Include phosphoketolase and phosphate transacetalyase reactions (PR #23).
  * Increase annotation coverage for metabolites and reactions by mapping to MetaNetX (PR #23).
  * Script for gene essentiality validation (PR #23).
* Fix
  * Remove duplicated reactions, introduced at various steps of reconstruction (PR #23).
  * Set fatty-CoA-ligases as irreversible, as they could otherwise form ATP regenerating loops that are unlikely to appear _in vivo_ (PR #23).
  * Refit energetics with new `NGAM` and `GAM` values (PR #23).
* Refactor
  * Change style of graph with growth rate validation (PR #23).
  * Refactor various simulation scripts (PR #23).
  * All new (non `yeast-GEM`) reactions are assigned identifiers with prefix `t_`  (PR #23).

### rhto v1.1.0:
* Features
  * Include model of _Yarrowia lipolytica_ [iYali4.1.1](https://github.com/SysBioChalmers/Yarrowia_lipolytica_W29-GEM/releases/tag/4.1.1) as template during reconstruction (PR #18).
  * Standardize identifiers: all new metabolites have ids `m_xxxx`, while new reactions are `t_xxxx` or `tr_xxxx` (for transport reactions) (PR #18).
  * Reduced the amount of hard-coded ids in scripts (PR #18).
* Fix
  * Include missing reactions such as complex I (PR #18).
  * Re-fit energetics (NGAM and GAM) (PR #18).
* Reorder
  * Reorganize reconstruction process: use TSV files for lists of curations (PR #18).
  * Include legend in validation graph (PR #18).

### rhto v1.0.0:
* Documentation:
  * Update readme.md to include correct links to Zenodo and bioRxiv (PR #15).

### rhto v0.0.2:
* Features
  * Include script `ComplementaryScripts\validation\growthPrediction.m` to valide model by comparing growth predictions (PR #13).
  * Include script `ComplementaryScripts\analysis\fseof_target_prediction.m` for FSEOF target prediction (PR #13).
* Fixes
  * Energetics script `r7_energetics.m` correctly reads file with flux data (PR #13).
  * Correct unlogical and duplicate grRules (PR #13).
  * Correct subSystems naming, removing orphan indices from template model (PR #13).
* Documentation
  * `Readme.md` has detailed text on how to cite the model (PR #13).
  * `Readme.md` contains version badge and Zenodo doi badge (PR #13).

### rhto v0.0.1:
Initial release of rhto-GEM.