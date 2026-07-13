"""Align a mapped reaction's reactant/product atom order to each other.

Background
----------
Given an atom-mapped reaction SMILES, `align_reaction` renumbers one side
(the "moving" side) so that atoms it shares with the other side (the
"fixed" side) are visited by `mol_to_smiles` in the same relative order
the fixed side visits them in -- across the whole molecule, not just from
a chosen root atom. Concretely, with `align_to="reactant"` (the default):

  1. `atom_visit_order(reactant_mol)` gives each reactant atom's rank in
     the order `mol_to_smiles` would write it.
  2. Each product atom that shares an atom-map number with a reactant atom
     is assigned that reactant atom's rank.
  3. Product atoms with no counterpart in the reactants (rare with a
     proper atom mapping, but handled the same way) are ranked after all
     mapped atoms, keeping their original relative order.
  4. Product atoms are renumbered by that rank and re-serialized with
     `mol_to_smiles`, so its ascending-index DFS traversal follows the
     reactants' atom order wherever the product's own connectivity allows.

This suits forward reaction prediction (reactants are the fixed model
input, the product is the generation target). `align_to="product"` runs
the same procedure with the roles swapped: the product is left as-is and
the reactants are reordered to follow it instead -- the R-SMILES/UAlign
convention for retrosynthesis (product fixed, reactants are the
generation target).

The fixed side is left in its input atom order (i.e. whatever root/order
the caller's SMILES already implies) unless `random=True` (see
`align_reaction`). Atom-map numbers are stripped from the output; agents
(if present) are re-serialized in their original order with atom maps
stripped, but are not reordered.

Order agreement between the two sides is only guaranteed *within* each
moving-side connected component (molecule): a fragment gets renumbered so
its own DFS follows the fixed side's order as closely as its internal
connectivity allows, but a fragment that is still disconnected from the
rest of the moving side (e.g. a second reagent that only gets bonded to
the main substrate on the other side) must still be written as one
contiguous block -- it cannot be interleaved mid-traversal with another
fragment the way it may be interleaved (via a branch) on the fixed side.
"""

from __future__ import annotations

import random as _random  # aliased: align_reaction has a `random` bool param
import re
import warnings
from collections import Counter

from rdkit import Chem

from idxsmiles import atom_visit_order, mol_to_smiles

# A reaction SMILES separator '>' is any '>' not immediately preceded by
# '-' -- the exclusion is needed because idxsmiles (and this module) also
# supports '->' dative bonds, whose '>' must not be mistaken for one of the
# two reactants/agents/products separators.
_REACTION_SEPARATOR = re.compile(r"(?<!-)>")


def _split_reaction_smiles(rxn_smiles: str) -> tuple[str, str, str]:
    """Split a reaction SMILES into (reactants, agents, products) on its two
    top-level '>' separators, ignoring '>' that is part of a '->' dative
    bond anywhere in the string."""
    seps = [m.start() for m in _REACTION_SEPARATOR.finditer(rxn_smiles)]
    if len(seps) != 2:
        raise ValueError(
            "expected a reaction SMILES in 'reactants>agents>products' "
            "form with exactly two non-dative-bond '>' separators (found "
            f"{len(seps)}): {rxn_smiles!r}")
    i, j = seps
    return rxn_smiles[:i], rxn_smiles[i + 1:j], rxn_smiles[j + 1:]


def _aligned_atom_order(moving_mol: Chem.Mol,
                        fixed_rank_by_map: dict[int, int]) -> list[int]:
    """Return a `Chem.RenumberAtoms`-style newOrder for `moving_mol`:
    newOrder[i] is the (pre-renumbering) atom index that should become
    atom i, so that atoms shared with the fixed side (by atom-map number)
    come out in the fixed side's rank order, and unshared atoms are pushed
    after them, keeping their original relative order."""
    n = moving_mol.GetNumAtoms()
    max_fixed_rank = max(fixed_rank_by_map.values(), default=-1)

    rank: list[int | None] = [None] * n
    for atom in moving_mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num and map_num in fixed_rank_by_map:
            rank[atom.GetIdx()] = fixed_rank_by_map[map_num]

    unranked = [i for i in range(n) if rank[i] is None]
    for offset, idx in enumerate(unranked):
        rank[idx] = max_fixed_rank + 1 + offset

    return sorted(range(n), key=lambda i: rank[i])


def _strip_atom_maps(mol: Chem.Mol) -> None:
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)


def _warn_about_unusable_atom_maps(
        reactant_mol: Chem.Mol, product_mol: Chem.Mol) -> None:
    """Warn (rather than fail) when the atom-map numbers can't actually
    drive an alignment -- align_reaction still runs in these cases, but
    silently degrades to leaving the moving side in its input order, which
    is easy to miss without a warning."""
    reactant_maps = [a.GetAtomMapNum() for a in reactant_mol.GetAtoms()
                     if a.GetAtomMapNum()]
    product_maps = [a.GetAtomMapNum() for a in product_mol.GetAtoms()
                    if a.GetAtomMapNum()]

    for side, maps in (("reactants", reactant_maps), ("product", product_maps)):
        dups = sorted(m for m, c in Counter(maps).items() if c > 1)
        if dups:
            warnings.warn(
                f"duplicate atom-map number(s) in {side}: {dups} -- "
                "alignment for those atoms is undefined", stacklevel=3)

    if not (set(reactant_maps) & set(product_maps)):
        warnings.warn(
            "no atom-map numbers shared between reactants and product -- "
            "neither side will be reordered", stacklevel=3)


def align_reaction(mapped_rxn_smiles: str, align_to: str = "reactant",
                   random: bool = False, seed: int | None = None) -> str:
    """Given a reaction SMILES with atom-mapped reactants and product(s)
    (`"reactants>agents>products"`, agents may be empty), return the
    equivalent unmapped reaction SMILES with one side reordered to follow
    the other's atom order (see module docstring).

    `align_to` picks which side stays fixed and which is reordered to
    match it: `"reactant"` (default) fixes the reactants and reorders the
    product (forward reaction prediction-style); `"product"` fixes the
    product and reorders the reactants instead (retrosynthesis-style, as
    in R-SMILES/UAlign).

    If `random` is true, the fixed side's atom indices are fully shuffled
    (seeded by `seed`, if given) before alignment, so it is written from a
    random root with random branch-order tie-breaks and the other side is
    aligned to *that* -- useful for generating multiple
    distinct-but-equivalent training examples per reaction (as in
    R-SMILES's root-augmentation). With `random=False` (the default) the
    fixed side keeps the atom order implied by the input SMILES.

    Warns (rather than raising) if the atom-map numbers can't drive an
    alignment: no atom-map numbers shared between reactants and product,
    or a duplicate atom-map number on one side. In both cases the affected
    atoms are left in their input relative order instead of being aligned.
    """
    if align_to not in ("product", "reactant"):
        raise ValueError(
            f"align_to must be 'product' or 'reactant', got {align_to!r}")

    reactants_smi, agents_smi, products_smi = _split_reaction_smiles(
        mapped_rxn_smiles)

    product_mol = Chem.MolFromSmiles(products_smi)
    if product_mol is None:
        raise ValueError(f"could not parse product SMILES: {products_smi!r}")
    reactant_mol = Chem.MolFromSmiles(reactants_smi)
    if reactant_mol is None:
        raise ValueError(
            f"could not parse reactant SMILES: {reactants_smi!r}")

    _warn_about_unusable_atom_maps(reactant_mol, product_mol)

    if align_to == "product":
        fixed_mol, moving_mol = product_mol, reactant_mol
    else:
        fixed_mol, moving_mol = reactant_mol, product_mol

    if random:
        order = list(range(fixed_mol.GetNumAtoms()))
        _random.Random(seed).shuffle(order)
        fixed_mol = Chem.RenumberAtoms(fixed_mol, order)

    fixed_rank = atom_visit_order(fixed_mol)
    fixed_rank_by_map = {
        atom.GetAtomMapNum(): fixed_rank[atom.GetIdx()]
        for atom in fixed_mol.GetAtoms() if atom.GetAtomMapNum()
    }
    new_order = _aligned_atom_order(moving_mol, fixed_rank_by_map)
    moving_mol = Chem.RenumberAtoms(moving_mol, new_order)

    if align_to == "product":
        product_mol, reactant_mol = fixed_mol, moving_mol
    else:
        reactant_mol, product_mol = fixed_mol, moving_mol

    _strip_atom_maps(reactant_mol)
    _strip_atom_maps(product_mol)

    aligned_agents = ""
    if agents_smi:
        agents_mol = Chem.MolFromSmiles(agents_smi)
        if agents_mol is None:
            raise ValueError(f"could not parse agents SMILES: {agents_smi!r}")
        _strip_atom_maps(agents_mol)
        aligned_agents = mol_to_smiles(agents_mol)

    return (f"{mol_to_smiles(reactant_mol)}>{aligned_agents}>"
            f"{mol_to_smiles(product_mol)}")
