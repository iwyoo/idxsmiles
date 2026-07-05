# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - Unreleased

Initial release.

### Added

- `mol_to_idx_ordered_smiles(mol)`: a DFS-based SMILES writer for RDKit
  `Mol` objects that visits atoms in strict ascending atom-index order at
  every branch point, instead of RDKit's canonical rank or the ring-biased
  order that `Chem.MolToSmiles(mol, canonical=False)` actually produces.
- Ring closures, including multi-digit `%(NN)` closures.
- Multi-fragment SMILES (`.`-separated) output.
- Aromaticity, charges, isotopes, and atom maps.
- Standard tetrahedral chirality (`@`/`@@`), including 2-real-neighbor,
  lone-pair-bearing stereocentres (e.g. chiral P, S).
- Double-bond cis/trans stereo (`/`, `\`).
- Dative bonds (`->`, `<-`).
- Test suite covering round-trip structural correctness under random atom
  renumbering, a regression test for RDKit's ring-bond traversal bias, and
  an independent check that atoms are actually emitted in ascending-index
  order.

### Limitations

- Non-tetrahedral chirality (square planar, trigonal bipyramidal,
  octahedral) is not supported and raises `NotImplementedError`.
- Cumulated/shared stereo double bonds (one single bond flanking two
  stereo double bonds at once, e.g. conjugated diazo systems) are not
  supported and raise `NotImplementedError`.

[0.1.0]: https://github.com/iwyoo/idxsmiles/releases/tag/v0.1.0
