# idxsmiles

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
pytest test_idxsmiles.py
```

The test suite checks round-trip structural correctness under random atom
renumbering across a curated set of molecules (rings, chirality, cis/trans,
dative bonds, multi-fragment), a dedicated regression test for the
ring-bond bias described above, and an independent check that atoms are
*actually* written in ascending-index order (not just that the round trip
happens to produce an equivalent structure) — verified by tagging atoms
with their index as an atom-map number and confirming the re-parsed output
reads them back in the same order as a separately implemented reference
traversal.

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
