import re
from typing import Dict, List

from .mol_utils import canonicalize_smiles, inchify, get_longest_smiles


# Reaction type constants
STANDARD = "STANDARD"
STANDARD_AGENT = "STANDARD_agent"
STANDARD_WITH_TILDE = "STANDARD_WITH_TILDE"
EXTENDED = "EXTENDED"
UNKNOWN = "UNKNOWN"


# Regex pattern to extract the fragment info from the extended info of reaction SMILES
_EXTENDED_FRAGMENT_REGEX = re.compile(r"f:[\d\.,]+")
# Regex pattern to extract the fragment groups from the fragment info
_FRAGMENT_GROUP_REGEX = re.compile(r"(\d+(?:\.\d+)*)")


def detect_rxn_type(rxn):
    """
    Auto-detect the reaction SMILES format.

    Returns one of:
        EXTENDED              - |f:0.2,3.4.5| fragment notation
        STANDARD_agent        - A.B>C>D  (with reagents/agents)
        STANDARD_WITH_TILDE   - A~B>C~D (tilde-separated fragments)
        STANDARD              - A.B>>D  (simple)
    """
    if not isinstance(rxn, str) or not rxn:
        return UNKNOWN

    # 1) EXTENDED: anything in |...| brackets is a fragment-group marker
    if re.search(r"\|[^|]*\|", rxn):
        return EXTENDED

    # 2) Tilde-separated fragments (only meaningful when tildes appear
    #    alongside the > delimiters)
    if "~" in rxn:
        return STANDARD_WITH_TILDE

    # 3) Use split count, not character count, because '>' also appears
    #    inside SMILES atom names (e.g. "O>") - we want the number of
    #    delimiter-separated parts, not the number of '>' chars.
    if rxn_type_from_split(rxn) == "agent":
        return STANDARD_AGENT
    if rxn_type_from_split(rxn) == "standard":
        return STANDARD

    return UNKNOWN


def rxn_type_from_split(rxn):
    """Inspect split('>') parts to distinguish STANDARD vs STANDARD_agent.

    Returns 'agent' for 3 non-empty parts, 'standard' for 2 parts (or 3 parts
    where the middle/reagents part is empty, which corresponds to '>>' in the
    source string), 'invalid' otherwise.
    """
    parts = rxn.split(">")
    if len(parts) == 3:
        # An empty middle part means the source contained '>>' (no reagents)
        if parts[1] == "":
            return "standard"
        return "agent"
    if len(parts) == 2:
        return "standard"
    return "invalid"


def _split_parts(rxn_type, rxn):
    """
    Split a reaction SMILES into (reactants_str, reagents_str, products_str).
    Strip EXTENDED fragment markers and translate tildes to dots so that the
    rest of the pipeline can treat all input uniformly.
    """
    if rxn_type == EXTENDED:
        # Remove all |...| markers
        cleaned = re.sub(r"\|[^|]*\|", "", rxn)
        parts = cleaned.split(">")
    elif rxn_type == STANDARD_WITH_TILDE:
        # Tildes are used to group fragment alternatives inside each part;
        # flatten them into the same dot-separated list for downstream code.
        parts = rxn.replace("~", ".").split(">")
    else:
        parts = rxn.split(">")

    if len(parts) == 3:
        reactants, reagents, products = parts
    elif len(parts) == 2:
        # STANDARD has no agents column -> treat as reactants, empty reagents, products
        reactants, products = parts
        reagents = ""
    else:
        raise ValueError(f"Invalid reaction SMILES: {rxn}")

    return reactants, reagents, products


def process_reaction(rxn, with_agent=True, rxn_type=None):
    """
    Process and canonicalize a reaction SMILES.

    Supports four reaction SMILES formats and auto-detects which one is
    supplied. The output is always normalised to the STANDARD format:
        reactants>>product
    Reagents/catalysts that also appear in products are removed from both
    sides (treated as agents).

    Args:
        rxn:       The reaction SMILES string.
        with_agent: If True (default), reagents/agents are merged into the
                    precursors and any molecule appearing in both precursors
                    and products is removed from both.
        rxn_type:  Optional. One of STANDARD / STANDARD_agent /
                   STANDARD_WITH_TILDE / EXTENDED. If None, the type is
                   detected automatically from the string.

    Returns:
        Canonical reaction SMILES in "precursors>>product" form, or "" on
        failure.
    """
    if rxn_type is None:
        rxn_type = detect_rxn_type(rxn)

    if rxn_type == UNKNOWN:
        return ""

    try:
        reactants, reagents, products = _split_parts(rxn_type, rxn)
    except ValueError:
        return ""

    try:
        if with_agent:
            precursors = [inchify(canonicalize_smiles(r)) for r in reactants.split(".")]
            if len(reagents) > 0:
                precursors += [
                    inchify(canonicalize_smiles(r)) for r in reagents.split(".")
                ]
            products = [inchify(canonicalize_smiles(p)) for p in products.split(".")]
            precursors_set = set(precursors)
            products_set = set(products)
            common_set = precursors_set & products_set
            products_set = products_set - common_set
            precursors = list(precursors_set)
            products = list(products_set)
        else:
            precursors = [inchify(canonicalize_smiles(r)) for r in reactants.split(".")]
            products = [inchify(canonicalize_smiles(p)) for p in products.split(".")]
    except Exception:
        return ""

    joined_precursors = ".".join(sorted(precursors, key=lambda s: len(s), reverse=True))
    joined_products = get_longest_smiles(
        ".".join(sorted(products, key=lambda s: len(s), reverse=True))
    )
    return f"{joined_precursors}>>{joined_products}"


# ---------------------------------------------------------------------------
# Extended / Tilde reaction SMILES conversion
# ---------------------------------------------------------------------------
#
# The functions below convert between the four reaction SMILES formats
# (STANDARD, STANDARD_agent, STANDARD_WITH_TILDE, EXTENDED) and the
# structured (reactants, agents, products) representation.  Design follows
# the conventions used by `rxn4chemistry/rxn-chemutils`,
# `MolecularAI/reaction_utils` and `rxn4chemistry/rxn-reaction-preprocessing`
# (see README "Acknowledgements").
# ---------------------------------------------------------------------------


def split_smiles_and_fragment_info(rxn: str):
    """
    Split an (extended) reaction SMILES into the pure SMILES part and the
    fragment info suffix.

    Example:
        "CC.O.[Na+].[Cl-]>>CCO |f:2.3|"
        -> ("CC.O.[Na+].[Cl-]>>CCO", "|f:2.3|")
    """
    m = re.search(r"^(\S+)\s*(.*)$", rxn.strip())
    if m is None:
        return rxn, ""
    return m.group(1), m.group(2)


def determine_fragment_groups(extended_info: str):
    """
    From a fragment info string like ``|f:0.2,5.6|`` return the list of
    groups of indices that must be merged, e.g. ``[[0, 2], [5, 6]]``.
    """
    match = _EXTENDED_FRAGMENT_REGEX.search(extended_info)
    if match is None:
        return []
    group_matches = _FRAGMENT_GROUP_REGEX.findall(match.group(0))
    return [[int(i) for i in g.split(".")] for g in group_matches]


def _split_reaction_smiles(rxn: str):
    """Split a reaction SMILES at the three reaction arrows ``>``, ``>``, ``>``.

    Lone ``>`` characters inside a SMILES (e.g. dative bonds like ``Cl-]``)
    are kept inside the current side and are NOT treated as reaction arrows.
    ``>>`` counts as two arrows and yields an empty agents side in the middle.
    """
    sides = []
    cur = []
    i = 0
    arrow_count = 0
    while i < len(rxn):
        ch = rxn[i]
        if ch == ">":
            if i + 1 < len(rxn) and rxn[i + 1] == ">":
                # '>>' ends the reactants side and starts an empty agents side
                sides.append("".join(cur))
                cur = []
                sides.append("")
                i += 2
                arrow_count = 2
                continue
            if arrow_count == 2:
                # End of agents, start of products
                sides.append("".join(cur))
                cur = []
                i += 1
                arrow_count = 3
                continue
            # Lone '>' in the SMILES (e.g. dative bond) - keep it
            cur.append(ch)
            i += 1
            continue
        cur.append(ch)
        i += 1
    sides.append("".join(cur))
    return sides


def _expand_dots(side: str):
    """Expand a side (reactants/reagents/products) into a list of fragments."""
    return [frag for frag in side.split(".") if frag]


def tilde_to_standard(rxn: str):
    """
    Convert a ``STANDARD_WITH_TILDE`` reaction SMILES to ``STANDARD`` form
    by replacing every ``~`` with ``.``.

    Example:
        ``CC.O.[Na+]~[Cl-]>>CCO`` -> ``CC.O.[Na+].[Cl-]>>CCO``
    """
    return rxn.replace("~", ".")


def standard_to_tilde(rxn: str):
    """
    Convert a ``STANDARD`` reaction SMILES to ``STANDARD_WITH_TILDE`` form
    by replacing ``.`` with ``~`` inside each side of the reaction
    (reactants, reagents, products are separated by '>').

    Example:
        ``CC.O.[Na+].[Cl-]>>CCO`` -> ``CC~O~[Na+]~[Cl-]>>CCO``
    """
    sides = _split_reaction_smiles(rxn)
    converted = ["~".join(_expand_dots(side)) for side in sides]
    return ">".join(converted)


def parse_tilde_reaction_smiles(rxn: str, remove_atom_maps: bool = True):
    """
    Parse a ``STANDARD_WITH_TILDE`` reaction SMILES into three lists of
    compounds.  Within each side, ``~`` joins fragments of the same
    compound, while ``.`` separates different compounds.  Multi-fragment
    compounds are placed at the end of each side, matching
    ``parse_extended_reaction_smiles``.

    Example:
        ``CC.O.[Na+]~[Cl-]>>CCO``
        -> (['CC', 'O', '[Na+].[Cl-]'], [], ['CCO'])
    """
    pure_smiles, _ = split_smiles_and_fragment_info(rxn)
    if remove_atom_maps:
        pure_smiles = re.sub(r":\d+\]", "]", pure_smiles)
    sides = _split_reaction_smiles(pure_smiles)
    parsed = []
    for side in sides:
        if not side:
            parsed.append([])
            continue
        compounds = side.split(".")
        merged = []
        singletons = []
        for compound in compounds:
            if "~" in compound:
                merged.append(compound.replace("~", "."))
            else:
                if compound:
                    singletons.append(compound)
        parsed.append(singletons + merged)
    while len(parsed) < 3:
        parsed.append([])
    return parsed[0], parsed[1], parsed[2]


def tilde_to_extended(rxn: str):
    """
    Convert a ``STANDARD_WITH_TILDE`` reaction SMILES to ``EXTENDED`` form.

    Example:
        ``CC.O.[Na+]~[Cl-]>>CCO``
        -> ``CC.O.[Na+].[Cl-]>>CCO |f:2.3|``
    """
    reactants, agents, products = parse_tilde_reaction_smiles(rxn)
    return to_extended_reaction_smiles(reactants, agents, products)


def extended_to_tilde(rxn: str):
    """
    Convert an ``EXTENDED`` reaction SMILES to ``STANDARD_WITH_TILDE`` form.

    Example:
        ``CC.O.[Na+].[Cl-]>>CCO |f:2.3|``
        -> ``CC.O.[Na+]~[Cl-]>>CCO``
    """
    reactants, agents, products = parse_extended_reaction_smiles(rxn)

    def _side_to_tilde(side):
        parts = []
        for compound in side:
            parts.append('~'.join(compound.split('.')))
        return '.'.join(parts)

    sides = [_side_to_tilde(reactants), _side_to_tilde(agents), _side_to_tilde(products)]
    return ">".join(sides)


def parse_extended_reaction_smiles(rxn: str, remove_atom_maps: bool = True):
    """
    Parse an EXTENDED reaction SMILES into a list of three sides
    (reactants, agents, products). Each side is a list of compounds, where
    compounds that share a fragment index are merged back into a single
    dot-separated string.

    Example:
        ``CC.O.[Na+].[Cl-]>>CCO |f:2.3|``
        -> (['CC', 'O', '[Na+].[Cl-]'], [], ['CCO'])
    """
    pure_smiles, fragment_info = split_smiles_and_fragment_info(rxn)
    if remove_atom_maps:
        pure_smiles = re.sub(r":\d+\]", "]", pure_smiles)

    sides = _split_reaction_smiles(pure_smiles)
    fragments_per_side = [_expand_dots(s) for s in sides]
    fragment_groups = determine_fragment_groups(fragment_info)

    return _merge_fragments(fragments_per_side, fragment_groups)


def to_extended_reaction_smiles(reactants, agents, products):
    """
    Build an EXTENDED reaction SMILES from three lists of compounds.
    Compounds that contain '.' (multi-fragment compounds, e.g. ``[Na+].[Cl-]``)
    are recorded in the fragment info at the end.

    Example:
        (['CC', 'O', '[Na+].[Cl-]'], [], ['CCO'])
        -> ``CC.O.[Na+].[Cl-]>>CCO |f:2.3|``
    """
    sides = [list(reactants), list(agents), list(products)]
    groups = []        # indices of fragments that belong to the same compound
    smiles_groups = [] # one dot-separated SMILES per side

    offset = 0
    for side in sides:
        side_fragments = []
        for compound in side:
            fragments = compound.split(".") if compound else [""]
            side_fragments.extend(fragments)
            if len(fragments) > 1:
                groups.append(list(range(offset, offset + len(fragments))))
            offset += len(fragments)
        smiles_groups.append(".".join(side_fragments))

    pure_smiles = ">".join(smiles_groups)
    if not groups:
        return pure_smiles
    fragment_info = "|f:" + ",".join(".".join(str(i) for i in g) for g in groups) + "|"
    return f"{pure_smiles} {fragment_info}"


def extended_to_standard(rxn: str):
    """
    Convert an EXTENDED reaction SMILES back to STANDARD form (just dropping
    the trailing fragment info suffix).

    Example:
        ``CC.O.[Na+].[Cl-]>>CCO |f:2.3|``
        -> ``CC.O.[Na+].[Cl-]>>CCO``
    """
    pure_smiles, _ = split_smiles_and_fragment_info(rxn)
    return pure_smiles


def standard_to_extended(rxn: str):
    """
    Build an EXTENDED reaction SMILES from a STANDARD reaction SMILES.
    Compounds that already contain '.' (i.e. the caller knows which fragments
    belong to the same compound) are recorded in the fragment info.

    Example:
        ``CC.O.[Na+].[Cl-]>>CCO``
        -> ``CC.O.[Na+].[Cl-]>>CCO |f:2.3|``
    """
    sides = _split_reaction_smiles(rxn)
    # Pad to 3 sides: reactants, agents, products
    while len(sides) < 3:
        sides.append("")

    groups = []        # indices of fragments belonging to the same compound
    smiles_groups = [] # one dot-separated SMILES per side
    offset = 0

    for side_str in sides:
        side_fragments = []
        for compound in side_str.split("."):
            if not compound:
                continue
            fragments = compound.split(".")
            side_fragments.extend(fragments)
            if len(fragments) > 1:
                groups.append(list(range(offset, offset + len(fragments))))
            offset += len(fragments)
        smiles_groups.append(".".join(side_fragments))

    pure_smiles = ">".join(smiles_groups)
    if not groups:
        return pure_smiles
    fragment_info = "|f:" + ",".join(".".join(str(i) for i in g) for g in groups) + "|"
    return f"{pure_smiles} {fragment_info}"


def _merge_fragments(fragments_per_side, fragment_groups):
    """
    Merge fragments according to fragment groups.

    For each side (a list of fragments), fragments whose indices appear
    together in a group are joined back into a single dot-separated string
    representing one compound.

    Returns a list of three lists: merged reactants, merged agents, merged
    products.  In each list, multi-fragment compounds are placed **after**
    single-fragment compounds, matching the convention used by
    ``rxn-chemutils``.
    """
    offset = 0
    merged_sides = []
    for fragments in fragments_per_side:
        n = len(fragments)
        grouped = {}   # relative index -> joined compound
        used = set()
        # Apply groups whose indices fall in this side
        for group in fragment_groups:
            relative = [i - offset for i in group]
            if all(0 <= r < n for r in relative):
                grouped[relative[0]] = ".".join(fragments[r] for r in relative)
                used.update(relative)
        # Build final list, putting grouped (multi-fragment) compounds at the end
        singletons = [fragments[i] for i in range(n) if i not in used]
        multi = [grouped[r] for r in sorted(grouped)]
        merged_sides.append(singletons + multi)
        offset += n
    return merged_sides


# ---------------------------------------------------------------------------
# Validity checks
# ---------------------------------------------------------------------------


def is_valid_reaction(rxn, rxn_type=None, allow_empty_products=False, is_strict=False):
    """
    Return ``True`` if ``rxn`` is a structurally well-formed reaction
    SMILES, ``False`` otherwise.

    A reaction is considered valid when:

    1. The string is non-empty (unless ``allow_empty_products`` is set
       for the products side).
    2. Its format is one of the four supported types (auto-detected if
       ``rxn_type`` is not given) — ``UNKNOWN`` always returns ``False``.
    3. After stripping EXTENDED fragment markers and expanding tildes,
       the SMILES splits into either 2 sides (``A.B>>D``) or 3 sides
       (``A.B>C>D``).
    4. Every individual fragment is a valid molecule SMILES according to
       :func:`canonicalize.mol_utils.is_valid_molecule`.

    When ``is_strict=True`` the following additional structural problems
    are also rejected:

    * Empty fragments (e.g. ``CCO..CCN``).
    * Duplicate compounds within a single side
      (e.g. ``CCO.CCO>>CCOCCO`` — both are ``CCO``).
    * Duplicate compounds across the three sides
      (e.g. ``CCO>>CCO``, ``CCO>>CCO.CCO`` or the same reagent on
      both sides of the arrow).
    * Same compound appearing as both a reactant and a product in the
      same reaction (no-op / self-redox), unless it is the solvent that
      legitimately appears on both sides (e.g. ``CCO>CCO>CCO``).

    Args:
        rxn: A reaction SMILES string.
        rxn_type: Optional. One of STANDARD / STANDARD_agent /
            STANDARD_WITH_TILDE / EXTENDED.  If ``None``, the type is
            detected automatically.
        allow_empty_products: If ``True``, the products side may be
            empty (which is sometimes legitimate for ``>reagents>``
            partial reaction templates).
        is_strict: If ``True``, perform the extra structural checks
            listed above.  Defaults to ``False`` (lenient mode).

    Returns:
        ``True`` if the reaction is valid, ``False`` otherwise.

    Examples:
        >>> is_valid_reaction("CCO.CCO>>CCOCCO")
        True
        >>> is_valid_reaction("CCO.CCO>>CCOCCO", is_strict=True)
        False
        >>> is_valid_reaction("CCO.CCO>O>CCOCCO")
        True
        >>> is_valid_reaction("CC.O.[Na+]~[Cl-]>>CCO")
        True
        >>> is_valid_reaction("CC.O.[Na+].[Cl-]>>CCO |f:2.3|")
        True
        >>> is_valid_reaction("not_a_reaction")
        False
    """
    # Local import to avoid circular import at module load time
    from .mol_utils import is_valid_molecule

    if not isinstance(rxn, str) or not rxn.strip():
        return False

    if rxn_type is None:
        rxn_type = detect_rxn_type(rxn)
    if rxn_type == UNKNOWN:
        return False

    # Strip fragment info and normalise tildes to dots so that we can
    # validate each fragment uniformly.
    try:
        reactants, reagents, products = _split_parts(rxn_type, rxn)
    except ValueError:
        return False

    sides = [reactants, reagents, products]
    if len(sides) < 2:
        return False

    # The reactants side must be non-empty for any reaction.
    if not sides[0]:
        return False
    if not sides[-1]:
        if not allow_empty_products:
            return False

    fragments_by_side: List[List[str]] = []
    for side in sides:
        if not side:
            fragments_by_side.append([])
            continue
        fragments = side.split(".")
        # Empty fragments (consecutive dots or leading/trailing dots)
        if any(not f for f in fragments):
            return False
        if not is_valid_molecule(list(fragments)):
            return False
        fragments_by_side.append(fragments)

    if not is_strict:
        return True

    # ------------------------------------------------------------------
    # Strict checks
    # ------------------------------------------------------------------
    canonical = lambda frags: [canonicalize_smiles(f) for f in frags]  # noqa: E731

    # 1. Duplicate compounds within a single side.
    for fragments in fragments_by_side:
        if not fragments:
            continue
        canon = canonical(fragments)
        if len(set(canon)) != len(canon):
            return False

    # 2. Duplicate compounds across the three sides — any fragment that
    #    appears more than once across (reactants, reagents, products)
    #    is rejected, *unless* it is a legitimate solvent shared by all
    #    three sides (e.g. ``CCO>CCO>CCO``).  The standard exception is
    #    when the same molecule appears in the three sides with the
    #    *exact same role* (solvent/mediator only).
    canon_sides = [canonical(frags) for frags in fragments_by_side]
    all_canon = canon_sides[0] + canon_sides[1] + canon_sides[-1]
    counts: Dict[str, int] = {}
    for c in all_canon:
        if not c:
            continue
        counts[c] = counts.get(c, 0) + 1
    # Allow the "all three sides are the same molecule" case (solvent
    # only) — that contributes 3 occurrences and is the only time a
    # count of 3 should be permitted.  Any other duplicate (count >= 2
    # that is not a multiple-of-3 solvent) is invalid.
    for c, n in counts.items():
        if n == 1:
            continue
        if n == 3 and all(
            c in side
            for side in canon_sides
        ):
            continue
        return False

    return True


