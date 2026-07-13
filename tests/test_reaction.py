"""Tests for idxsmiles.reaction.align_reaction.

Covers:
  - structural round-trip preservation of reactants, agents, and product
    (canonical-SMILES equivalence, order of '.'-separated fragments and
    atom-map numbers stripped away)
  - the actual *alignment* property, in both directions (`align_to`):
    atoms shared between the "moving" side and the "fixed" side are
    written in the same relative order in both -- checked by tagging
    atoms with isotopes (survives atom-map stripping) instead of
    hardcoding expected output strings
  - the fixed side (whichever `align_to` picks) is serialized exactly as
    idxsmiles' own writer would serialize the unmodified input mol -- it
    is never reordered
  - known, structurally-unavoidable limitations of the alignment: a
    moving-side fragment disconnected from the rest of the moving side
    can't be interleaved mid-traversal the way it may be on the fixed
    side, and more generally a DFS can't visit a child before its parent
    even if the fixed side's order would ask for that
  - leaving-group (unmapped) atoms are pushed after mapped atoms within
    their own fragment
  - agents pass through with atom maps stripped but order untouched
  - `random=True` is reproducible given the same `seed`, can produce
    different (but still structurally correct) output across seeds, and
    `align_to` validation
  - '>' inside a '->' dative bond is not mistaken for a reaction-level
    separator
  - malformed input (wrong number of reaction separators, unparsable
    fragment) raises ValueError
  - warnings (rather than failures) when atom-map numbers can't drive an
    alignment: no overlap between reactants and product, or a duplicate
    atom-map number on one side

Run with: pytest tests/test_reaction.py
"""
import warnings

from rdkit import Chem

import pytest

from idxsmiles import mol_to_smiles
from idxsmiles.reaction import align_reaction, _split_reaction_smiles


def _canon(smi: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))


def _isotope_tagged(mapped_rxn_smiles: str) -> str:
    """Copy of `mapped_rxn_smiles` with each mapped atom's atom-map number
    also stamped on as an isotope. `align_reaction` only strips atom-map
    numbers, so the isotope survives into the unmapped output and lets a
    test recover, from the *output* alone, which input atom each output
    atom was -- without relying on (or duplicating) align_reaction's own
    ranking logic. Re-serialized with idxsmiles' own `mol_to_smiles`
    (rather than RDKit's canonical writer) so the input's atom order is
    preserved rather than silently re-canonicalized before it ever reaches
    align_reaction."""
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
        tagged.append(mol_to_smiles(mol))
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


@pytest.mark.parametrize("align_to", ["reactant", "product"])
@pytest.mark.parametrize("name,rxn,expected_reactants,expected_product",
                         CURATED_REACTIONS)
def test_structural_roundtrip(name, rxn, expected_reactants, expected_product,
                              align_to):
    result = align_reaction(rxn, align_to=align_to)
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


# (name, align_to) pairs known to hit the DFS parent-before-child
# limitation documented in test_dfs_parent_before_child_limitation: in
# each, the fixed side's rank order puts a bridging atom's own child (the
# far branch reached only through it) at a *lower* rank than the bridging
# atom itself, which no DFS can honor. Not a bug -- see that test.
_KNOWN_PARENT_BEFORE_CHILD_LIMITATION = {
    ("esterification", "reactant"),
    ("williamson_ether", "reactant"),
    ("chiral_ester", "reactant"),
}


@pytest.mark.parametrize("align_to", ["reactant", "product"])
@pytest.mark.parametrize("name,rxn,_er,_ep", CURATED_REACTIONS)
def test_moving_side_follows_fixed_side_order_within_each_fragment(
        name, rxn, _er, _ep, align_to):
    """The real guarantee `align_reaction` makes: within any one
    moving-side connected component, atoms shared with the fixed side come
    out in the same relative order the fixed side visits them in.
    (Cross-fragment order is *not* guaranteed -- see
    test_disconnected_fragment_limitation and
    test_dfs_parent_before_child_limitation.)"""
    if (name, align_to) in _KNOWN_PARENT_BEFORE_CHILD_LIMITATION:
        pytest.skip("known DFS parent-before-child limitation, see "
                    "test_dfs_parent_before_child_limitation")
    tagged_rxn = _isotope_tagged(rxn)
    result = align_reaction(tagged_rxn, align_to=align_to)
    reactants, _agents, product = _split_reaction_smiles(result)
    fixed_side, moving_side = (
        (product, reactants) if align_to == "product" else (reactants, product))

    fixed_rank = {m: r for r, m in enumerate(_written_map_nums(fixed_side))}

    for fragment in moving_side.split("."):
        written = [m for m in _written_map_nums(fragment) if m in fixed_rank]
        assert written == sorted(written, key=lambda m: fixed_rank[m]), (
            f"{name}/align_to={align_to}: fragment {fragment!r} does not "
            "follow the fixed side's order")


@pytest.mark.parametrize("name,rxn,_er,_ep", CURATED_REACTIONS)
def test_fixed_side_is_unreordered_when_aligning_to_product(name, rxn, _er, _ep):
    _reactants, _agents, products = _split_reaction_smiles(rxn)
    product_mol = Chem.MolFromSmiles(products)
    for atom in product_mol.GetAtoms():
        atom.SetAtomMapNum(0)
    expected = mol_to_smiles(product_mol)

    result = align_reaction(rxn, align_to="product")
    _r, _a, actual_product = _split_reaction_smiles(result)
    assert actual_product == expected


@pytest.mark.parametrize("name,rxn,_er,_ep", CURATED_REACTIONS)
def test_fixed_side_is_unreordered_when_aligning_to_reactant(name, rxn, _er, _ep):
    reactants, _agents, _products = _split_reaction_smiles(rxn)
    reactant_mol = Chem.MolFromSmiles(reactants)
    for atom in reactant_mol.GetAtoms():
        atom.SetAtomMapNum(0)
    expected = mol_to_smiles(reactant_mol)

    result = align_reaction(rxn)  # default align_to="reactant"
    actual_reactants, _a, _p = _split_reaction_smiles(result)
    assert actual_reactants == expected


def test_disconnected_fragment_limitation():
    """Documents a known limitation of align_to='product': the methanol
    fragment bonds into the *middle* of the aryl ring's own traversal in
    the product (as a branch between ring atoms 4 and 6), but since it is
    a separate molecule in the reactants, it cannot be interleaved there
    -- it must still come out as one contiguous fragment, so the *global*
    atom order does not fully match the product's, even though each
    fragment's internal order does (see the general order test above)."""
    rxn = ("[cH:1]1[cH:2][cH:3][c:4]([Br:5])[cH:6][cH:7]1.[CH3:8][OH:9]"
          ">>[cH:1]1[cH:2][cH:3][c:4]([O:9][CH3:8])[cH:6][cH:7]1")
    tagged_rxn = _isotope_tagged(rxn)
    result = align_reaction(tagged_rxn, align_to="product")
    reactants, _agents, product = _split_reaction_smiles(result)

    reactant_order = _written_map_nums(reactants)
    product_order = _written_map_nums(product)
    shared = set(reactant_order) & set(product_order)
    global_reactant_order = [m for m in reactant_order if m in shared]
    global_product_order = [m for m in product_order if m in shared]
    assert global_reactant_order != global_product_order


def test_dfs_parent_before_child_limitation():
    """Documents a more general limitation, not specific to disconnected
    fragments: with align_to='reactant', the product's ester oxygen
    (map 7) bridges two reactant-derived branches -- the ethanol-derived
    C6-C5 chain and the acid-derived C2(=O3)-C1 chain. The fixed
    (reactant) side's rank order would require visiting O3 (map 3) before
    its own parent C2 (map 2), which no DFS can do -- align_reaction still
    visits C2 before O3, but picks that atom's lowest-target-rank child
    first (O3 before C1), the closest achievable approximation."""
    rxn = ("[OH:7][CH2:6][CH3:5].[OH:4][C:3](=[O:2])[CH3:1]"
          ">>[CH3:1][C:2](=[O:3])[O:7][CH2:6][CH3:5]")
    tagged_rxn = _isotope_tagged(rxn)
    result = align_reaction(tagged_rxn)  # default align_to="reactant"
    _reactants, _agents, product = _split_reaction_smiles(result)
    assert _written_map_nums(product) == [7, 6, 5, 2, 3, 1]


def test_leaving_group_is_ordered_after_mapped_atoms_in_its_fragment():
    rxn = "[CH3:1][CH2:2][Br:3].[NH2:4][CH3:5]>>[CH3:1][CH2:2][NH:4][CH3:5]"
    tagged_rxn = _isotope_tagged(rxn)
    result = align_reaction(tagged_rxn, align_to="product")
    reactants, _agents, _product = _split_reaction_smiles(result)

    bromide_fragment = next(f for f in reactants.split(".") if "Br" in f)
    written = _written_map_nums(bromide_fragment)
    # map 3 (Br) has no counterpart in the product -- it must be last.
    assert written[-1] == 3
    assert written[:-1] == sorted(written[:-1])


def test_invalid_align_to_raises():
    with pytest.raises(ValueError, match="align_to"):
        align_reaction("CC(=O)O.CCO>>CC(=O)OCC", align_to="bogus")


def test_random_seed_is_reproducible():
    rxn = ("[cH:1]1[cH:2][cH:3][c:4]([Br:5])[cH:6][cH:7]1.[CH3:8][OH:9]"
          ">>[cH:1]1[cH:2][cH:3][c:4]([O:9][CH3:8])[cH:6][cH:7]1")
    assert (align_reaction(rxn, random=True, seed=42)
           == align_reaction(rxn, random=True, seed=42))


def test_random_true_can_produce_multiple_distinct_orderings():
    rxn = ("[cH:1]1[cH:2][cH:3][c:4]([Br:5])[cH:6][cH:7]1.[CH3:8][OH:9]"
          ">>[cH:1]1[cH:2][cH:3][c:4]([O:9][CH3:8])[cH:6][cH:7]1")
    results = {align_reaction(rxn, random=True, seed=s) for s in range(20)}
    assert len(results) > 1


@pytest.mark.parametrize("align_to", ["reactant", "product"])
@pytest.mark.parametrize("name,rxn,expected_reactants,expected_product",
                         CURATED_REACTIONS)
def test_random_still_structurally_correct(
        name, rxn, expected_reactants, expected_product, align_to):
    for seed in (0, 1, 2):
        result = align_reaction(rxn, align_to=align_to, random=True, seed=seed)
        reactants, _agents, product = _split_reaction_smiles(result)
        assert _canon(reactants) == _canon(expected_reactants)
        assert _canon(product) == _canon(expected_product)


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


def test_warns_when_no_atom_maps_at_all():
    with pytest.warns(UserWarning, match="no atom-map numbers shared"):
        align_reaction("CC(=O)O.CCO>>CC(=O)OCC")


def test_warns_when_only_one_side_is_mapped():
    with pytest.warns(UserWarning, match="no atom-map numbers shared"):
        align_reaction("[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][CH2:6][OH:7]"
                       ">>CC(=O)OCC")
    with pytest.warns(UserWarning, match="no atom-map numbers shared"):
        align_reaction("CC(=O)O.CCO>>[CH3:1][C:2](=[O:3])[O:7][CH2:6][CH3:5]")


def test_warns_on_duplicate_atom_map_number():
    rxn = ("[CH3:1][C:2](=[O:3])[OH:1].[CH3:5][CH2:6][OH:7]"
          ">>[CH3:1][C:2](=[O:3])[O:7][CH2:6][CH3:5]")
    with pytest.warns(UserWarning, match="duplicate atom-map number"):
        align_reaction(rxn)


@pytest.mark.parametrize("align_to", ["reactant", "product"])
@pytest.mark.parametrize("name,rxn,_er,_ep", CURATED_REACTIONS)
def test_well_formed_input_does_not_warn(name, rxn, _er, _ep, align_to):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        align_reaction(rxn, align_to=align_to)
