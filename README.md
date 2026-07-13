# idxsmiles

[![PyPI](https://img.shields.io/pypi/v/idxsmiles.svg)](https://pypi.org/project/idxsmiles/)
[![CI](https://github.com/iwyoo/idxsmiles/actions/workflows/ci.yml/badge.svg)](https://github.com/iwyoo/idxsmiles/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

An RDKit SMILES writer that orders atoms by their **original atom index**
instead of RDKit's canonical rank — deterministic, human-aligned output,
with stereochemistry preserved.

## Why

`Chem.MolToSmiles(mol)` (canonical=True) picks the next atom to visit using
canonical invariant ranks (`Canon::rankMolAtoms`), which have no relationship
to the atom's index in the `Mol` object. The resulting SMILES atom order can
look arbitrary relative to how the molecule was built or numbered.

`Chem.MolToSmiles(mol, canonical=False)` looks like the fix — internally its
`ranks` array is just `[0, 1, ..., N-1]` (the atom indices themselves). But
`Canon::dfsBuildStack` (RDKit's C++ traversal code) doesn't sort candidates
by that rank alone: if a candidate bond is part of a ring, a large
bond-type-dependent term gets added to its sort key, so a **lower-index ring
neighbor can still be visited after a higher-index non-ring neighbor**. This
is reproducible: renumber `C1=CC(C)C1` so the ring-closure atom gets index 3
and the exocyclic methyl gets index 4, and `canonical=False` still emits
atom 4 before atom 3.

Because that bias is baked into RDKit's compiled writer and can't be
disabled from Python, this library implements its own DFS-based SMILES
writer: atoms are always visited in strict ascending atom-index order, with
no exceptions for ring bonds — while still producing valid, stereo-correct
SMILES that RDKit can re-parse into the identical structure.

## Installation

```bash
pip install idxsmiles
```

## Usage

```python
from rdkit import Chem
from idxsmiles import mol_to_smiles

mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
smiles = mol_to_smiles(mol)
```

`mol_to_smiles(mol: Chem.Mol) -> str` takes an RDKit `Mol` and
returns a SMILES string whose atom order follows the input `Mol`'s atom
indices as closely as the molecular graph allows (ring/branch structure
still constrains what's possible — the result is not necessarily perfectly
monotonic, but every branch point prefers the lowest available index).
Re-parsing the output reproduces the same molecule, including
stereochemistry.

`atom_visit_order(mol: Chem.Mol) -> list[int]` returns, for each atom
index, its 0-based rank in that same traversal — i.e. the position it
would appear at in `mol_to_smiles`'s output — without serializing anything.

## Reaction SMILES alignment

`idxsmiles.reaction.align_reaction` takes an atom-mapped reaction SMILES
(`"reactants>agents>products"`, agents may be empty) and returns the
unmapped reaction with one side reordered so that atoms it shares with the
other side are visited in the same relative order the other side visits
them in — not just from a shared root atom, but across as much of the
full atom order as connectivity allows.

```python
from idxsmiles.reaction import align_reaction

mapped = ("[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][CH2:6][OH:7]"
          ">>[CH3:1][C:2](=[O:3])[O:7][CH2:6][CH3:5]")
align_reaction(mapped)
# 'CC(=O)O.CCO>>CC(=O)OCC'
align_reaction(mapped, align_to="product")
# 'CC(=O)O.OCC>>CC(=O)OCC'
```

`align_to` picks which side stays fixed and which is reordered to match
it:
- `"reactant"` (default): the reactants are left as-is and the **product**
  is reordered to follow them — for forward reaction prediction (reactants
  are the fixed model input, the product is the generation target).
- `"product"`: the product is left as-is and the **reactants** are
  reordered to follow it — the R-SMILES/UAlign convention for
  retrosynthesis (product fixed, reactants are the generation target).
  Note how the reactants' ethanol fragment (`OCC`) comes out as a literal
  substring match of the product's tail (`OCC`) in the second example
  above — that's the alignment.

`random=True` fully shuffles the fixed side's atom indices (seeded by
`seed`, if given) before aligning, so it's written from a random root with
random branch-order tie-breaks each time and the other side follows *that*
— useful for generating multiple distinct-but-equivalent training examples
per reaction, as in R-SMILES's root-augmentation:

```python
align_reaction(mapped, random=True, seed=0)
align_reaction(mapped, random=True, seed=1)  # a different, still-valid alignment
```

Other behavior:
- Atoms on the moving side with no counterpart on the fixed side (leaving
  groups, unmapped atoms) are ordered after all mapped atoms, keeping
  their original relative order.
- Agents, if present, are re-serialized in their original order with
  atom-map numbers stripped, but are not reordered.
- Atom-map numbers are stripped from the returned SMILES.
- A `'>'` inside a `'->'` dative bond is not mistaken for a
  reactants/agents/products separator.
- Warns (rather than raising) instead of silently producing an unaligned
  result when the atom-map numbers can't drive an alignment: no atom-map
  numbers shared between reactants and product, or a duplicate atom-map
  number on one side.

Order agreement between the two sides is only guaranteed *within* each
moving-side connected component: a fragment that only bonds to the rest
of the molecule on the other side (e.g. a second reagent joining
mid-chain) still has to be written as one contiguous block, even where the
fixed side interleaves it via a branch — and more generally, a DFS can
never visit a child before its own parent, even where the fixed side's
order would ask for that. Both are structural limitations of SMILES
itself, not bugs (see `test_disconnected_fragment_limitation` and
`test_dfs_parent_before_child_limitation` in `tests/test_reaction.py`).

## What's handled

- Ring closures (including multi-digit `%(NN)` closures)
- Multi-fragment SMILES (`.`-separated)
- Aromaticity, charges, isotopes, atom maps
- Standard tetrahedral chirality (`@`/`@@`), including 2-real-neighbor,
  lone-pair-bearing stereocentres (e.g. chiral P, S)
- Double-bond cis/trans stereo (`/`, `\`)
- Dative bonds (`->`, `<-`)

## Out of scope

Non-tetrahedral chirality (square planar, trigonal bipyramidal, octahedral)
and cumulated/shared stereo double bonds (one single bond flanking two
stereo double bonds at once, e.g. conjugated diazo systems) are not
supported. These raise `NotImplementedError` rather than silently emitting
incorrect stereochemistry.

## Testing

```bash
pip install rdkit pytest
pytest
```

`tests/test_idxsmiles.py` checks round-trip structural correctness under
random atom renumbering across a curated set of molecules (rings,
chirality, cis/trans, dative bonds, multi-fragment), a dedicated
regression test for the ring-bond bias described above, and an independent
check that atoms are *actually* written in ascending-index order (not just
that the round trip happens to produce an equivalent structure) — verified
by tagging atoms with their index as an atom-map number and confirming the
re-parsed output reads them back in the same order as a separately
implemented reference traversal.

`tests/test_reaction.py` checks `align_reaction` the same way, in both
`align_to` directions: structural round-trip correctness across a curated
set of reactions, and — since `align_reaction` strips atom-map numbers
from its output — an isotope-tag trick to independently confirm the
*ordering* itself (that shared atoms really do come out in the fixed
side's relative order within each moving-side fragment), plus
`random`/`seed`, the atom-map warnings, the dative-bond-separator, and
malformed-input edge cases.

## How it works

RDKit is still used for the low-level pieces:
- `Atom.GetSmarts()` / `Bond.GetSmarts()` for atom/bond text (these happen
  to call straight into RDKit's own SMILES atom/bond writer,
  `SmilesWrite::GetAtomSmiles` / `GetBondSmiles`)
- The chirality- and double-bond-stereo bookkeeping is a Python port of the
  logic in `Canon::canonicalizeFragment` / `Canon::canonicalizeDoubleBond`,
  adapted to run against this library's own traversal order instead of
  RDKit's canonical one.

Everything else — the DFS traversal order, ring-closure digit assignment,
and branch/parenthesis structure — is implemented independently in
`idxsmiles.py`.

## License

MIT — see [LICENSE](LICENSE).
