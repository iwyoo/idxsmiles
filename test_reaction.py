"""Tests for idxsmiles.reaction.align_reaction.

Covers:
  - structural round-trip preservation of reactants, agents, and product
    (canonical-SMILES equivalence, order of '.'-separated fragments and
    atom-map numbers stripped away)
  - the actual *alignment* property: atoms shared between a reactant
    fragment and the product are written in the same relative order in
    both -- checked by tagging atoms with isotopes (survives atom-map
    stripping) instead of hardcoding expected output strings
  - a reactant fragment that only bonds to the rest of the molecule *in*
    the product (e.g. a second reagent) is still written as one
    contiguous block, even where the product interleaves it via a branch
    -- a documented limitation, not a bug
  - leaving-group (unmapped) atoms are pushed after mapped atoms within
    their own fragment
  - agents pass through with atom maps stripped but order untouched
  - '>' inside a '->' dative bond is not mistaken for a reaction-level
    separator
  - malformed input (wrong number of reaction separators, unparsable
    fragment) raises ValueError

Run with: pytest test_reaction.py
"""
from rdkit import Chem

import pytest

from idxsmiles.reaction import align_reaction, _split_reaction_smiles


def _canon(smi: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))


def _isotope_tagged(mapped_rxn_smiles: str) -> str:
    """Copy of `mapped_rxn_smiles` with each mapped atom's atom-map number
    also stamped on as an isotope. `align_reaction` only strips atom-map
    numbers, so the isotope survives into the unmapped output and lets a
    test recover, from the *output* alone, which input atom each output
    atom was -- without relying on (or duplicating) align_reaction's own
    ranking logic."""
    reactants, agents, products = _split_reaction_smiles(mapped_rxn_smiles)
    tagged = []
    for part in (reactants, agents, products):
        if not part:
            tagged.append(part)
            continue
        mol = Chem.MolFromSmiles(part)
        assert mol is not None
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum():
                atom.SetIsotope(atom.GetAtomMapNum())
        tagged.append(Chem.MolToSmiles(mol))
    return ">".join(tagged)


def _written_map_nums(smiles_side: str) -> list[int]:
    """Isotope numbers of a (possibly multi-fragment) SMILES side, in the
    order its atoms are written -- RDKit parses left-to-right, so this is
    literally the textual order."""
    mol = Chem.MolFromSmiles(smiles_side, sanitize=False)
    assert mol is not None
    return [atom.GetIsotope() for atom in mol.GetAtoms() if atom.GetIsotope()]


# (label, mapped reaction SMILES, expected unmapped reactants, expected
# unmapped product) -- expected sides are checked as whole-string canonical
# SMILES, which RDKit normalizes independently of '.'-fragment order.
CURATED_REACTIONS = [
    ("esterification",
     "[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][CH2:6][OH:7]"
     ">>[CH3:1][C:2](=[O:3])[O:7][CH2:6][CH3:5]",
     "CC(=O)O.CCO", "CC(=O)OCC"),
    ("williamson_ether",
     "[cH:1]1[cH:2][cH:3][c:4]([Br:5])[cH:6][cH:7]1.[CH3:8][OH:9]"
     ">>[cH:1]1[cH:2][cH:3][c:4]([O:9][CH3:8])[cH:6][cH:7]1",
     "Brc1ccccc1.CO", "COc1ccccc1"),
    ("amide_coupling",
     "[CH3:1][C:2](=[O:3])[OH:4].[NH2:5][CH3:6]"
     ">>[CH3:1][C:2](=[O:3])[NH:5][CH3:6]",
     "CC(=O)O.NC", "CC(=O)NC"),
    ("sn2_amine",
     "[CH3:1][CH2:2][Br:3].[NH2:4][CH3:5]>>[CH3:1][CH2:2][NH:4][CH3:5]",
     "CCBr.NC", "CCNC"),
    ("chiral_ester",
     "[CH3:1][C@@H:2]([NH2:3])[C:4](=[O:5])[OH:6].[CH3:7][OH:8]"
     ">>[CH3:1][C@@H:2]([NH2:3])[C:4](=[O:5])[O:8][CH3:7]",
     "C[C@@H](N)C(=O)O.CO", "C[C@@H](N)C(=O)OC"),
    ("dative_bond_spectator",
     "[Fe+2]->1CC1.[CH3:1][OH:2]>>[Fe+2]->1CC1.[CH3:1][O:2]C",
     "[Fe+2]->1CC1.CO", "[Fe+2]->1CC1.COC"),
]


@pytest.mark.parametrize("name,rxn,expected_reactants,expected_product",
                         CURATED_REACTIONS)
def test_structural_roundtrip(name, rxn, expected_reactants, expected_product):
    result = align_reaction(rxn)
    reactants, _agents, product = _split_reaction_smiles(result)
    assert _canon(reactants) == _canon(expected_reactants)
    assert _canon(product) == _canon(expected_product)


@pytest.mark.parametrize("name,rxn,_er,_ep", CURATED_REACTIONS)
def test_no_atom_maps_in_output(name, rxn, _er, _ep):
    result = align_reaction(rxn)
    assert ":" not in result  # ':' only appears in atom-map notation


def test_agents_pass_through_unreordered_with_maps_stripped():
    rxn = ("[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][CH2:6][OH:7]"
          ">[Na+].[OH-:9]>[CH3:1][C:2](=[O:3])[O:7][CH2:6][CH3:5]")
    result = align_reaction(rxn)
    _reactants, agents, _product = _split_reaction_smiles(result)
    assert _canon(agents) == _canon("[Na+].[OH-]")
    assert ":" not in agents


@pytest.mark.parametrize("name,rxn,_er,_ep", CURATED_REACTIONS)
def test_mapped_atoms_follow_product_order_within_each_fragment(
        name, rxn, _er, _ep):
    """The real guarantee `align_reaction` makes: within any one reactant
    connected component, atoms shared with the product come out in the
    same relative order the product visits them in. (Cross-fragment order
    is *not* guaranteed -- see test_disconnected_fragment_limitation.)"""
    tagged_rxn = _isotope_tagged(rxn)
    result = align_reaction(tagged_rxn)
    reactants, _agents, product = _split_reaction_smiles(result)

    product_rank = {map_num: rank for rank, map_num in
                    enumerate(_written_map_nums(product))}

    for fragment in reactants.split("."):
        written = [m for m in _written_map_nums(fragment)
                  if m in product_rank]
        assert written == sorted(written, key=lambda m: product_rank[m]), (
            f"{name}: fragment {fragment!r} does not follow product order")


def test_disconnected_fragment_limitation():
    """Documents the known limitation: the methanol fragment bonds into
    the *middle* of the aryl ring's own traversal in the product (as a
    branch between ring atoms 4 and 6), but since it is a separate
    molecule in the reactants, it cannot be interleaved there -- it must
    still come out as one contiguous fragment, so the *global* atom order
    does not fully match the product's, even though each fragment's
    internal order does (see the general order test above)."""
    rxn = ("[cH:1]1[cH:2][cH:3][c:4]([Br:5])[cH:6][cH:7]1.[CH3:8][OH:9]"
          ">>[cH:1]1[cH:2][cH:3][c:4]([O:9][CH3:8])[cH:6][cH:7]1")
    tagged_rxn = _isotope_tagged(rxn)
    result = align_reaction(tagged_rxn)
    reactants, _agents, product = _split_reaction_smiles(result)

    reactant_order = _written_map_nums(reactants)
    product_order = _written_map_nums(product)
    shared = set(reactant_order) & set(product_order)
    global_reactant_order = [m for m in reactant_order if m in shared]
    global_product_order = [m for m in product_order if m in shared]
    assert global_reactant_order != global_product_order


def test_leaving_group_is_ordered_after_mapped_atoms_in_its_fragment():
    rxn = "[CH3:1][CH2:2][Br:3].[NH2:4][CH3:5]>>[CH3:1][CH2:2][NH:4][CH3:5]"
    tagged_rxn = _isotope_tagged(rxn)
    result = align_reaction(tagged_rxn)
    reactants, _agents, _product = _split_reaction_smiles(result)

    bromide_fragment = next(f for f in reactants.split(".") if "Br" in f)
    written = _written_map_nums(bromide_fragment)
    # map 3 (Br) has no counterpart in the product -- it must be last.
    assert written[-1] == 3
    assert written[:-1] == sorted(written[:-1])


@pytest.mark.parametrize("bad_rxn", [
    "no separators at all",
    "reactants>only_one_separator",
    "a>b>c>d",
])
def test_malformed_separator_count_raises(bad_rxn):
    with pytest.raises(ValueError):
        align_reaction(bad_rxn)


@pytest.mark.parametrize("bad_rxn", [
    "not_a_smiles>>C",
    "C>>not_a_smiles",
    "C>not_a_smiles>C",
])
def test_unparsable_fragment_raises(bad_rxn):
    with pytest.raises(ValueError):
        align_reaction(bad_rxn)


def test_dative_bond_not_mistaken_for_separator():
    """'->' inside a fragment must survive splitting -- regression test for
    the naive `.split('>')` bug where a dative bond's '>' was mistaken for
    a reaction-level separator (idxsmiles supports '->'/'<-' dative bonds,
    see mol_to_smiles's module docstring)."""
    rxn = "[Fe+2]->1CC1.[CH3:1][OH:2]>>[Fe+2]->1CC1.[CH3:1][O:2]C"
    result = align_reaction(rxn)
    reactants, _agents, product = _split_reaction_smiles(result)
    assert _canon(reactants) == _canon("[Fe+2]->1CC1.CO")
    assert _canon(product) == _canon("[Fe+2]->1CC1.COC")
