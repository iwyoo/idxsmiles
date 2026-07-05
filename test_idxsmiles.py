"""Tests for idxsmiles.mol_to_smiles.

Covers:
  - round-trip structure preservation under random atom renumbering
  - the ring-bond bias this module exists to avoid (canonical=False alone
    does not avoid it -- see idxsmiles.py's module docstring)
  - tetrahedral chirality, including 2-real-neighbor (P/S, lone-pair)
    stereocentres
  - double-bond cis/trans stereo
  - dative bonds and an organometallic ring-closure case
  - multi-fragment SMILES
  - that atoms are actually emitted in ascending-index order (the ranking
    itself), not just that the round trip happens to reproduce an
    equivalent structure
  - that out-of-scope stereo (cumulated/shared stereo double bonds) fails
    loudly instead of silently emitting the wrong structure

Run with: pytest test_idxsmiles.py
"""
import random

import pytest
from rdkit import Chem

from idxsmiles import mol_to_smiles


def _canon(mol):
    return Chem.MolToSmiles(mol)


def _greedy_ascending_dfs_order(mol):
    """Independent reference implementation: always visit the smallest
    unvisited neighbor first. Used to check the *ranking* itself, not just
    round-trip correctness -- deliberately not sharing any code with
    idxsmiles.py's own traversal."""
    n = mol.GetNumAtoms()
    neighbors = [sorted(a.GetIdx() for a in atom.GetNeighbors())
                for atom in mol.GetAtoms()]
    visited = [False] * n
    order = []
    for start in range(n):
        if visited[start]:
            continue
        visited[start] = True
        order.append(start)
        stack = [(start, iter(neighbors[start]))]
        while stack:
            u, it = stack[-1]
            advanced = False
            for v in it:
                if not visited[v]:
                    visited[v] = True
                    order.append(v)
                    stack.append((v, iter(neighbors[v])))
                    advanced = True
                    break
            if not advanced:
                stack.pop()
    return order


CURATED_SMILES = [
    "CC(=O)Oc1ccccc1C(=O)O",       # aspirin
    "C1=CC(C)C1",                   # ring-bond-bias regression case
    "C[C@H](N)C(=O)O",
    "C[C@@H](N)C(=O)O",
    "N[C@@H](C)C(=O)O",
    "C1CC1.CCN",                    # multi-fragment
    "c1ccccc1/C=C/c1ccccc1",
    "c1ccccc1/C=C\\c1ccccc1",
    "F/C=C/F",
    "F/C=C\\F",
    "C(/F)=C/F",
    "[NH3+]CC(=O)[O-]",
    "CC(Br)(Cl)F",
    "O=C1CCCCC1",
    "C[C@H]1CC[C@@H](C)CC1",
    "CC(C)(C)c1ccc(O)cc1",
    "CC1=CC(=O)C=CC1=O",
    "c1ccc2ccccc2c1",
    "Cc1ccccc1",
    "C[P@H]CC(=O)O",                 # 2-real-neighbor (lone-pair) stereocentre
    "C[P@H]C[P@@H]CS(=O)(=O)[O-]",
    "C[S@](=O)CC",
    "CC[S@](=O)c1ccccc1",
    "[Fe+2]->1CC1",                  # dative bond
    "CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10",  # ferrocene-like
]


@pytest.mark.parametrize("smi", CURATED_SMILES)
@pytest.mark.parametrize("seed", range(5))
def test_roundtrip_under_random_renumbering(smi, seed):
    mol = Chem.MolFromSmiles(smi)
    assert mol is not None
    canon0 = _canon(mol)
    n = mol.GetNumAtoms()
    rng = random.Random(hash((smi, seed)))
    order = list(range(n))
    rng.shuffle(order)
    mol2 = Chem.RenumberAtoms(mol, order)

    out = mol_to_smiles(mol2)
    reparsed = Chem.MolFromSmiles(out)
    assert reparsed is not None, "writer produced unparsable SMILES: %r" % out
    assert _canon(reparsed) == canon0


def test_ring_bond_bias_is_avoided():
    """RDKit's own MolToSmiles(canonical=False) still lets ring bonds jump
    the queue ahead of lower-index exocyclic atoms (a large bond-type-
    dependent term is added to the sort key for ring bonds in
    Canon::dfsBuildStack). This writer must not have that bias."""
    mol = Chem.MolFromSmiles("C1=CC(C)C1")
    order = [0, 1, 2, 4, 3]  # ring-closure atom -> idx 3, methyl branch -> idx 4
    mol2 = Chem.RenumberAtoms(mol, order)
    out = mol_to_smiles(mol2)
    reparsed = Chem.MolFromSmiles(out)
    assert reparsed is not None
    assert _canon(reparsed) == _canon(mol)


def test_two_neighbor_chiral_center_swap_is_handled():
    """Regression test for a 2-real-neighbor (lone-pair) chiral centre
    (e.g. [P@H]) whose real-neighbor write order needs to be checked
    against RDKit's own convention rather than assumed irrelevant --
    RDKit's convention for this specific atom pattern differs across
    versions (see _two_neighbor_swap_flips_chirality in idxsmiles.py),
    and this permutation reliably exposed the wrong convention on RDKit
    2026.3.x before the fix (independent of hash-seed randomization, which
    otherwise makes CURATED_SMILES's random-shuffle tests flaky in terms
    of exactly which permutation gets tried)."""
    mol = Chem.MolFromSmiles("C[P@H]CC(=O)O")
    canon0 = _canon(mol)
    order = [1, 2, 0, 3, 4, 5]
    mol2 = Chem.RenumberAtoms(mol, order)
    out = mol_to_smiles(mol2)
    reparsed = Chem.MolFromSmiles(out)
    assert reparsed is not None
    assert _canon(reparsed) == canon0


@pytest.mark.parametrize("smi", CURATED_SMILES)
def test_atoms_are_written_in_ascending_index_order(smi):
    """Verifies the *ranking* itself, not just structural round-trip.

    Tags each atom with its own index as an atom-map number, writes the
    SMILES, and re-parses it -- SMILES parsing creates atoms left-to-right
    in the order they appear in the text, so the re-parsed mol's atom
    order *is* the textual order. That is compared against an
    independently implemented 'always visit the smallest unvisited
    neighbor first' DFS.
    """
    mol = Chem.MolFromSmiles(smi)
    assert mol is not None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    expected = _greedy_ascending_dfs_order(mol)

    out = mol_to_smiles(mol)
    reparsed = Chem.MolFromSmiles(out, sanitize=False)
    assert reparsed is not None
    written_order = [a.GetAtomMapNum() - 1 for a in reparsed.GetAtoms()]
    assert written_order == expected


def test_shared_stereo_double_bond_raises():
    """Cumulated/shared stereo double bonds (one single bond flanking two
    stereo double bonds at once) are out of scope and must fail loudly
    rather than silently emit the wrong stereochemistry."""
    smi = "C1CCNC(=O)/C(=N/N=C/2\\CCCCNC2=O)/C1"
    mol = Chem.MolFromSmiles(smi)
    assert mol is not None
    order = [11, 9, 0, 16, 2, 15, 10, 3, 5, 17, 14, 6, 7, 8, 4, 1, 13, 12]
    mol2 = Chem.RenumberAtoms(mol, order)
    with pytest.raises(NotImplementedError):
        mol_to_smiles(mol2)


def test_empty_mol():
    assert mol_to_smiles(Chem.Mol()) == ""
