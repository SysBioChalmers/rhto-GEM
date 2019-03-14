# History

### rhto v1.1.0:
* Features
  * Include model of _Yarrowia lipolytica_ [iYali4.1.1](https://github.com/SysBioChalmers/Yarrowia_lipolytica_W29-GEM/releases/tag/4.1.1) as template during reconstruction (PR #18).
  * Standardize identifiers: all new metabolites have ids `m_xxxx`, while new reactions are `t_xxxx` or `tr_xxxx` (for transport reactions) (PR #18).
  * Reduced the amount of hard-coded ids in scripts (PR #18).
* Fix
  * Include missing reactions such as complex I (PR #18).
  * Re-fit energetics (NGAM and GAM) (PR #18).
* Reorder:
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