# Changelog
All notable changes to this project will be documented in this file.

## [1.1.2] 2026-05-30
Support .cif files as input (auto-detected and converted).
Formally support AF3 predicted structures (including multimers).

## [1.1.2] - 2026-05-24
Adapt to AF3 per-atom pLDDT: use CA atom pLDDT and `resindices` for
cross-chain safe residue filtering in domain parsing.

## [1.1.1] - 2026-05-22

Reduce memory usage of result presentation.
Add `--pae-only` option to fetch script for downloading PAE files separately.

## [1.1.0] - 2026-05-14
Add symmetry transformation. Users can specify custom rotation/translation
operations via JSON files; DomainSeeker generates symmetry-related density
copies and extends crosslink data accordingly.

### [1.0.1] - 2026-04-13
More detailed status bar in ChimeraX while calculating prior probabilities.

### [1.0.1] - 2026-04-09
Adpat Mac OS.

### [1.0.1] - 2025-12-16
Add the *max_domain_size* to the domain parsing module.

### [1.0.1] - 2025-12-05
**release**  
Show running progress in ChimeraX status bar. 

## [1.0.0] - 2025-11-24
**Initial release.**