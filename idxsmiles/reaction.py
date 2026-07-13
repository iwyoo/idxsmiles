"""Align a mapped reaction's reactant atom order to its product's atom order.

Background
----------
Given an atom-mapped reaction SMILES, `align_reaction` renumbers reactant
atoms so that atoms shared with the product are visited by `mol_to_smiles`
in the same relative order the product atoms are -- across the whole
molecule, not just from a chosen root atom. Concretely:

  1. `atom_visit_order(product_mol)` gives each product atom's rank in the
     order `mol_to_smiles` would write it.
  2. Each reactant atom that shares an atom-map number with a product atom
     is assigned that product atom's rank.
  3. Reactant atoms with no counterpart in the product (leaving groups,
     unmapped atoms) are ranked after all mapped atoms, keeping their
     original relative order.
  4. Reactant atoms are renumbered by that rank and re-serialized with
     `mol_to_smiles`, so its ascending-index DFS traversal follows the
     product's atom order wherever the reactant's own connectivity allows.

The product side is left in its input atom order (i.e. whatever root/order
the caller's SMILES already implies) -- this function only aligns reactants
to it. Atom-map numbers are stripped from the output; agents (if present)
are re-serialized in their original order with atom maps stripped, but are
not reordered.

Order agreement with the product is only guaranteed *within* each reactant
connected component (molecule): a reactant fragment gets renumbered so its
own DFS follows the product's order as closely as its internal connectivity
allows, but a fragment that is still disconnected from the rest of the
reactants (e.g. a second reagent that only gets bonded to the main
substrate in the product) must still be written as one contiguous block --
it cannot be interleaved mid-traversal with another fragment the way it may
be interleaved (via a branch) in the product.
"""

from __future__ import annotations

import re

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


def _reactant_atom_order(reactant_mol: Chem.Mol,
                         product_rank_by_map: dict[int, int]) -> list[int]:
    """Return a `Chem.RenumberAtoms`-style newOrder: newOrder[i] is the
    (pre-renumbering) atom index that should become atom i."""
    n = reactant_mol.GetNumAtoms()
    max_product_rank = max(product_rank_by_map.values(), default=-1)

    rank: list[int | None] = [None] * n
    for atom in reactant_mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num and map_num in product_rank_by_map:
            rank[atom.GetIdx()] = product_rank_by_map[map_num]

    unranked = [i for i in range(n) if rank[i] is None]
    for offset, idx in enumerate(unranked):
        rank[idx] = max_product_rank + 1 + offset

    return sorted(range(n), key=lambda i: rank[i])


def _strip_atom_maps(mol: Chem.Mol) -> None:
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)


def align_reaction(mapped_rxn_smiles: str) -> str:
    """Given a reaction SMILES with atom-mapped reactants and product(s)
    (`"reactants>agents>products"`, agents may be empty), return the
    equivalent unmapped reaction SMILES with reactant atoms reordered to
    follow the product's atom order (see module docstring)."""
    reactants_smi, agents_smi, products_smi = _split_reaction_smiles(
        mapped_rxn_smiles)

    product_mol = Chem.MolFromSmiles(products_smi)
    if product_mol is None:
        raise ValueError(f"could not parse product SMILES: {products_smi!r}")
    reactant_mol = Chem.MolFromSmiles(reactants_smi)
    if reactant_mol is None:
        raise ValueError(
            f"could not parse reactant SMILES: {reactants_smi!r}")

    product_rank = atom_visit_order(product_mol)
    product_rank_by_map = {
        atom.GetAtomMapNum(): product_rank[atom.GetIdx()]
        for atom in product_mol.GetAtoms() if atom.GetAtomMapNum()
    }

    new_order = _reactant_atom_order(reactant_mol, product_rank_by_map)
    reactant_mol = Chem.RenumberAtoms(reactant_mol, new_order)

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
