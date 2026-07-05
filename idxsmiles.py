"""Custom SMILES writer that visits atoms in strict ascending atom-index order.

Background
----------
RDKit's MolToSmiles(canonical=True) picks the next atom to visit using
canonical invariant ranks (Canon::rankMolAtoms), which have nothing to do
with the atom's index in the Mol object, so the SMILES atom order can look
arbitrary.

MolToSmiles(canonical=False) looks like it should fix this: internally the
`ranks` array is filled with plain atom indices (`std::iota`) instead of the
canonical invariant. But `Canon::dfsBuildStack` (Code/GraphMol/Canon.cpp)
does NOT sort candidate atoms by that rank alone -- if the candidate bond is
part of a ring, a large bond-type-dependent term gets added to the sort key:

    auto rank = ranks[otherIdx];
    if (bond is in a ring) rank += (MAX_BONDTYPE - bondType) * MAX_NATOMS**2;

so a lower-index ring neighbor can be visited *after* a higher-index
non-ring neighbor (verified: renumbering 'C1=CC(C)C1' so the ring-closure
atom gets index 3 and the exocyclic methyl gets index 4 makes
canonical=False still emit atom 4 before atom 3). Since that bias lives in
the compiled C++ writer and can't be disabled from Python, this module
implements its own DFS-based writer: atoms are always visited in strict
ascending atom-index order, no exceptions for ring bonds.

RDKit is still used for:
  - per-atom / per-bond text (Atom.GetSmarts() / Bond.GetSmarts() happen to
    reuse RDKit's own SMILES atom/bond writer -- see
    SmilesWrite::GetAtomSmiles / GetBondSmiles), and
  - reproducing the chirality- and double-bond-stereo bookkeeping from
    Canon::canonicalizeFragment / Canon::canonicalizeDoubleBond, so that
    @/@@ and /,\\ come out correct for *our* traversal order rather than
    RDKit's.

Only standard tetrahedral chirality (CHI_TETRAHEDRAL_CW/CCW) and simple
double-bond cis/trans are handled. Non-tetrahedral stereo (square planar,
trigonal bipyramidal, octahedral) and cumulated/shared stereo double bonds
are out of scope and raise NotImplementedError.
"""

import heapq

from rdkit import Chem

_FLIP_DIR = {"/": "\\", "\\": "/"}


def _count_swaps(probe, ref):
    """Parity (# of transpositions) needed to turn `probe` into `ref`'s
    order. Mirrors RDKit's Atom::getPerturbationOrder."""
    order = [ref.index(x) for x in probe]
    swaps = 0
    for i in range(len(order)):
        while order[i] != i:
            j = order[i]
            order[i], order[j] = order[j], order[i]
            swaps += 1
    return swaps


def _atom_has_fourth_valence(atom):
    return atom.GetNumExplicitHs() == 1 or atom.GetImplicitValence() == 1


def _atom_is_unsaturated(atom):
    return any(b.GetBondTypeAsDouble() > 1 for b in atom.GetBonds())


def _chiral_needs_tag_inversion(atom, is_atom_first, num_closures):
    """Port of Canon::chiralAtomNeedsTagInversion: accounts for the implicit
    H of a 3-explicit-neighbor chiral atom being written as [C@H] instead of
    RDKit's internal 'H is the last neighbor' convention."""
    if atom.GetDegree() != 3:
        return False
    if is_atom_first and atom.GetNumExplicitHs() == 1:
        return True
    return (not _atom_has_fourth_valence(atom) and num_closures == 1 and
            not _atom_is_unsaturated(atom))


_two_neighbor_swap_flips_chirality_cached = None


def _two_neighbor_swap_flips_chirality():
    """Whether the installed RDKit treats the two real neighbors of a
    2-real-neighbor (lone-pair) stereocentre with an explicit/implicit H
    (e.g. a chiral [P@H] or [S@H]) as order-sensitive for @/@@ meaning.

    This is not a fixed RDKit convention: verified empirically that RDKit
    <= 2025.9.x parses '[P@H](C)CC(=O)O' and '[P@H](CC(=O)O)C' (same
    stereocentre, real neighbors listed in opposite order) as the *same*
    chirality, while RDKit >= 2026.3.x parses them as *opposite* -- a
    behavior change in RDKit's own SMILES handling of this atom pattern,
    not a bug in this module's traversal. Rather than hardcode either
    answer (and silently produce wrong output on whichever RDKit version
    guessed wrong), detect the installed RDKit's actual behavior once
    against this reference molecule and cache the result.
    """
    global _two_neighbor_swap_flips_chirality_cached
    if _two_neighbor_swap_flips_chirality_cached is None:
        neighbors_first = Chem.MolToSmiles(Chem.MolFromSmiles("[P@H](C)CC(=O)O"))
        neighbors_swapped = Chem.MolToSmiles(Chem.MolFromSmiles("[P@H](CC(=O)O)C"))
        _two_neighbor_swap_flips_chirality_cached = (
            neighbors_first != neighbors_swapped)
    return _two_neighbor_swap_flips_chirality_cached


def _permutes_real_bonds_for_chirality(atom):
    """Whether reordering this chiral atom's 2 real bonds (relative to
    their native Atom::getBonds() order) should flip its @/@@ tag.

    True for the common case (3 or 4 real neighbors, chirality defined
    directly by their listed order, matching Canon::canonicalizeFragment's
    getPerturbationOrder). For a 2-real-neighbor stereocentre with an
    explicit/implicit H, deferred to the installed RDKit's own observed
    convention -- see _two_neighbor_swap_flips_chirality.
    """
    if atom.GetDegree() != 2 or not _atom_has_fourth_valence(atom):
        return True
    return _two_neighbor_swap_flips_chirality()


def mol_to_smiles(mol: Chem.Mol) -> str:
    """Serialize `mol` to SMILES, visiting atoms in ascending atom-index
    order at every branch point (instead of RDKit's canonical or
    ring-biased non-canonical order). Re-parsing the result reproduces the
    same molecular structure, including stereochemistry."""
    work = Chem.RWMol(mol)
    n = work.GetNumAtoms()
    if n == 0:
        return ""

    if not work.HasProp("_StereochemDone"):
        Chem.AssignStereochemistry(work, cleanIt=False, force=False)

    # Bond.GetSmarts() renders '/'/'\\' straight from Bond::getBondDir(),
    # regardless of our traversal direction; clear any pre-existing marks so
    # only the ones *we* compute below (stereo_bond_char) show up.
    for bond in work.GetBonds():
        bond.SetBondDir(Chem.BondDir.NONE)

    neighbors = [sorted(a.GetIdx() for a in atom.GetNeighbors())
                for atom in work.GetAtoms()]

    visited = [False] * n
    parent_atom = [-1] * n
    parent_bond = [-1] * n
    visit_order = [-1] * n
    children = [[] for _ in range(n)]
    order_counter = 0
    fragment_roots = []

    # ---- pass 1: atom-idx-ordered DFS spanning forest ----------------
    for start in range(n):
        if visited[start]:
            continue
        fragment_roots.append(start)
        visited[start] = True
        visit_order[start] = order_counter
        order_counter += 1
        stack = [(start, iter(neighbors[start]))]
        while stack:
            u, it = stack[-1]
            advanced = False
            for v in it:
                if not visited[v]:
                    bond = work.GetBondBetweenAtoms(u, v)
                    parent_atom[v] = u
                    parent_bond[v] = bond.GetIdx()
                    children[u].append(v)
                    visited[v] = True
                    visit_order[v] = order_counter
                    order_counter += 1
                    stack.append((v, iter(neighbors[v])))
                    advanced = True
                    break
            if not advanced:
                stack.pop()

    tree_bonds = {b for b in parent_bond if b >= 0}
    ring_bonds = [bond.GetIdx() for bond in work.GetBonds()
                 if bond.GetIdx() not in tree_bonds]

    # ---- pass 2: ring-closure digit assignment -----------------------
    def endpoints_in_visit_order(bond_idx):
        bond = work.GetBondWithIdx(bond_idx)
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        return (a, b) if visit_order[a] < visit_order[b] else (b, a)

    ring_open_close = {}
    events = []
    for bidx in ring_bonds:
        opener, closer = endpoints_in_visit_order(bidx)
        ring_open_close[bidx] = (opener, closer)
        events.append((visit_order[opener], 0, bidx))
        events.append((visit_order[closer], 1, bidx))
    events.sort(key=lambda e: (e[0], e[1]))

    free_digits = []
    next_digit = 1
    digit_of_bond = {}
    for _, kind, bidx in events:
        if kind == 0:
            d = heapq.heappop(free_digits) if free_digits else next_digit
            if not free_digits or d == next_digit:
                next_digit = max(next_digit, d + 1)
            digit_of_bond[bidx] = d
        else:
            heapq.heappush(free_digits, digit_of_bond[bidx])

    opens_at = {a: [] for a in range(n)}
    closes_at = {a: [] for a in range(n)}
    for bidx, (opener, closer) in ring_open_close.items():
        opens_at[opener].append(bidx)
        closes_at[closer].append(bidx)
    for lst in list(opens_at.values()) + list(closes_at.values()):
        lst.sort()

    def ring_digit_text(digit):
        if digit < 10:
            return str(digit)
        if digit < 100:
            return "%" + str(digit)
        return "%(" + str(digit) + ")"

    # ---- pass 3: stereo bookkeeping -----------------------------------
    # 3a. tetrahedral chirality: figure out, per atom, the bond-index order
    # we will actually emit (ring closures first, then children ascending),
    # compare it to the atom's native bond order, and invert @ / @@ on our
    # working copy if the permutation parity is odd (matching
    # Canon::canonicalizeFragment).
    fragment_root_set = set(fragment_roots)
    for atom in work.GetAtoms():
        if atom.GetChiralTag() not in (Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                                       Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                raise NotImplementedError(
                    "only tetrahedral chirality is supported, atom %d has %s"
                    % (atom.GetIdx(), atom.GetChiralTag()))
            continue
        idx = atom.GetIdx()
        true_order = []
        if parent_bond[idx] >= 0:
            true_order.append(parent_bond[idx])
        for bidx in opens_at[idx] + closes_at[idx]:
            true_order.append(bidx)
        for child in children[idx]:
            true_order.append(work.GetBondBetweenAtoms(idx, child).GetIdx())
        # the parent bond (if any) is written *before* the atom symbol in
        # SMILES ("...P(A)..."), so for permutation purposes it must be
        # moved to wherever it naturally sits in the neighbor list -- it is
        # simplest to just treat it as coming first, matching RDKit's
        # convention of "the bond we arrived on is neighbor 0".
        ref_order = [b.GetIdx() for b in atom.GetBonds()]
        n_swaps = (_count_swaps(true_order, ref_order)
                  if _permutes_real_bonds_for_chirality(atom) else 0)
        is_first = idx in fragment_root_set
        num_closures = len(opens_at[idx]) + len(closes_at[idx])
        if _chiral_needs_tag_inversion(atom, is_first, num_closures):
            n_swaps += 1
        if n_swaps % 2:
            atom.InvertChirality()

    # 3b. double-bond cis/trans: for each stereo double bond, pick one
    # flanking (non-double) bond per side and compute the '/'/'\\'
    # character to emit for it, given *our* traversal order.
    stereo_bond_char = {}
    for bond in work.GetBonds():
        if (bond.GetBondType() != Chem.BondType.DOUBLE or
                bond.GetStereo() in (Chem.BondStereo.STEREONONE,
                                     Chem.BondStereo.STEREOANY)):
            continue
        stereo_atoms = list(bond.GetStereoAtoms())
        if len(stereo_atoms) != 2:
            continue
        # GetStereoAtoms()[0] is a neighbor of the begin atom, [1] a neighbor
        # of the end atom. a1 must be whichever double-bond atom *we* visit
        # first, which need not be the begin atom.
        begin_idx, end_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ref_begin, ref_end = stereo_atoms
        if visit_order[begin_idx] < visit_order[end_idx]:
            a1, a2, ref_a, ref_b = begin_idx, end_idx, ref_begin, ref_end
        else:
            a1, a2, ref_a, ref_b = end_idx, begin_idx, ref_end, ref_begin
        same_side = bond.GetStereo() in (Chem.BondStereo.STEREOZ,
                                         Chem.BondStereo.STEREOCIS)

        def pick_flanking_bond(sp2_idx, other_sp2_idx):
            candidates = []
            if parent_atom[sp2_idx] >= 0 and parent_atom[sp2_idx] != other_sp2_idx:
                candidates.append(parent_bond[sp2_idx])
            for bidx in opens_at[sp2_idx] + closes_at[sp2_idx]:
                b = work.GetBondWithIdx(bidx)
                other = b.GetOtherAtomIdx(sp2_idx)
                if other != other_sp2_idx:
                    candidates.append(bidx)
            for child in children[sp2_idx]:
                if child != other_sp2_idx:
                    candidates.append(work.GetBondBetweenAtoms(sp2_idx, child).GetIdx())
            plain = [c for c in candidates if c not in ring_open_close]
            return plain[0] if plain else (candidates[0] if candidates else None)

        bond_a = pick_flanking_bond(a1, a2)
        bond_b = pick_flanking_bond(a2, a1)
        if bond_a is None or bond_b is None:
            continue

        # The two sides use opposite comparison directions (this asymmetry
        # matches Canon::canonicalizeDoubleBond's isFirstFromAtom1/2Flipped):
        # on the a1 (first-visited) side we ask "is a1 before its anchor?";
        # on the a2 (second-visited) side we ask "is the anchor before a2?".
        def anchor_and_flip(sp2_idx, flank_bond_idx, sp2_first):
            b = work.GetBondWithIdx(flank_bond_idx)
            anchor = b.GetOtherAtomIdx(sp2_idx)
            is_ring = flank_bond_idx in ring_open_close
            if sp2_first:
                ordered_before = visit_order[sp2_idx] < visit_order[anchor]
            else:
                ordered_before = visit_order[anchor] < visit_order[sp2_idx]
            flipped = ordered_before != is_ring
            return anchor, flipped

        anchor_a, flip_a = anchor_and_flip(a1, bond_a, True)
        anchor_b, flip_b = anchor_and_flip(a2, bond_b, False)

        # parity(anchor == reference atom) on each side: if the flanking
        # atom we actually used differs from bond.GetStereoAtoms()'s
        # reference, it sits on the opposite branch, which flips the
        # relationship once.
        parity_a = 0 if anchor_a == ref_a else 1
        parity_b = 0 if anchor_b == ref_b else 1
        effective_same_side = same_side == (parity_a == parity_b)

        char_a = "/" if not flip_a else "\\"
        char_b_base = "/" if not effective_same_side else "\\"
        char_b = char_b_base if not flip_b else _FLIP_DIR[char_b_base]

        if bond_a in stereo_bond_char and stereo_bond_char[bond_a] != char_a:
            raise NotImplementedError(
                "conflicting bond direction for shared/cumulated stereo bond %d"
                % bond_a)
        if bond_b in stereo_bond_char and stereo_bond_char[bond_b] != char_b:
            raise NotImplementedError(
                "conflicting bond direction for shared/cumulated stereo bond %d"
                % bond_b)
        stereo_bond_char[bond_a] = char_a
        stereo_bond_char[bond_b] = char_b

    # ---- pass 4: text emission -----------------------------------------
    out = []

    def bond_text(bond_idx, from_idx):
        # Bond.GetSmarts() always renders a DATIVE bond as "->" (it assumes
        # the bond's own begin atom is to the left); flip it if we are
        # actually writing the bond starting from the end atom.
        bond = work.GetBondWithIdx(bond_idx)
        text = bond.GetSmarts()
        if (bond.GetBondType() == Chem.BondType.DATIVE and
                from_idx != bond.GetBeginAtomIdx()):
            text = "<-" if text == "->" else "->"
        return text

    def emit_bond(bond_idx, from_idx):
        if bond_idx in stereo_bond_char:
            out.append(stereo_bond_char[bond_idx])
        else:
            out.append(bond_text(bond_idx, from_idx))

    def emit_atom(idx):
        out.append(work.GetAtomWithIdx(idx).GetSmarts())
        for bidx in opens_at[idx] + closes_at[idx]:
            if bidx in stereo_bond_char:
                out.append(stereo_bond_char[bidx])
            elif bidx in opens_at[idx]:
                out.append(bond_text(bidx, idx))
            out.append(ring_digit_text(digit_of_bond[bidx]))
        kids = children[idx]
        for i, child in enumerate(kids):
            bidx = work.GetBondBetweenAtoms(idx, child).GetIdx()
            last = (i == len(kids) - 1)
            if not last:
                out.append("(")
            emit_bond(bidx, idx)
            emit_atom(child)
            if not last:
                out.append(")")

    frag_smiles = []
    for root in fragment_roots:
        out = []
        emit_atom(root)
        frag_smiles.append("".join(out))
    return ".".join(frag_smiles)
