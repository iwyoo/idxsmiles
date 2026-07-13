# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-07-13

### Added

- `idxsmiles.reaction.align_reaction`: given an atom-mapped reaction SMILES
  (`"reactants>agents>products"`), reorders one side's atoms to follow the
  other side's atom order -- not just from a shared root atom, but across
  as much of the full atom order as connectivity allows -- then
  re-serializes with `mol_to_smiles` and strips atom maps.
  - `align_to="reactant"` (default): fixes the reactants and reorders the
    product to match them (forward reaction prediction).
  - `align_to="product"`: fixes the product and reorders the reactants
    instead (retrosynthesis, the R-SMILES/UAlign convention).
  - `random=True` (with an optional `seed`) fully shuffles the fixed
    side's atom indices before aligning, for R-SMILES-style
    root-augmentation.
  - Warns (rather than silently degrading to an unaligned passthrough)
    when the atom-map numbers can't drive an alignment: no atom-map
    numbers shared between reactants and product, or a duplicate
    atom-map number on one side.
  - Handles multi-fragment reactants/products, leaving groups, agents,
    and dative bonds (`'>'` inside `'->'` is not mistaken for a
    reactants/agents/products separator).
- `atom_visit_order(mol)`: returns each atom's 0-based rank in the
  ascending atom-index DFS traversal `mol_to_smiles` writes atoms in,
  without serializing anything.

### Changed

- Moved the test suite into `tests/` (was two files at the repository
  root); CI updated accordingly.

## [0.1.1] - 2026-07-06

### Added

- Full type hints on all public and internal functions, plus a `py.typed`
  marker (PEP 561) so type checkers pick up the package's types.

### Changed

- Restructured the package from a flat `idxsmiles.py` module into an
  `idxsmiles/` package (`idxsmiles/__init__.py`) so that `py.typed` can be
  shipped per PEP 561. The public API (`from idxsmiles import
  mol_to_smiles`) is unchanged.
- Declared the MIT license in `pyproject.toml` using the classic
  `license = {text = "MIT"}` table form and added a matching
  `License :: OSI Approved :: MIT License` classifier, so PyPI displays the
  license correctly (it previously showed as unset).

## [0.1.0] - 2026-07-05

Initial release.

### Added

- `mol_to_smiles(mol)`: a DFS-based SMILES writer for RDKit `Mol` objects
  that visits atoms in strict ascending atom-index order at every branch
  point, instead of RDKit's canonical rank or the ring-biased order that
  `Chem.MolToSmiles(mol, canonical=False)` actually produces.
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

[0.2.0]: https://github.com/iwyoo/idxsmiles/releases/tag/v0.2.0
[0.1.1]: https://github.com/iwyoo/idxsmiles/releases/tag/v0.1.1
[0.1.0]: https://github.com/iwyoo/idxsmiles/releases/tag/v0.1.0
