# RxKit

Chemical reaction and molecule SMILES canonicalization toolkit powered by RDKit.

## Acknowledgements / References

The reaction SMILES handling in this package (in particular the conversion
between `STANDARD`, `STANDARD_AGENT`, `STANDARD_WITH_TILDE` and `EXTENDED`
formats, the `|f:...|` fragment-info parser, and the multi-component compound
merging logic) is inspired by and adapted from the following open-source
projects:

- [`rxn4chemistry/rxn-chemutils`](https://github.com/rxn4chemistry/rxn-chemutils) — reaction SMILES utilities, including extended-reaction SMILES parsing and conversion.
- [`MolecularAI/reaction_utils`](https://github.com/MolecularAI/reaction_utils/tree/main) — reaction preprocessing utilities from Molecular AI (AstraZeneca).
- [`rxn4chemistry/rxn-reaction-preprocessing`](https://github.com/rxn4chemistry/rxn-reaction-preprocessing/tree/main) — reaction preprocessing pipeline used by the IBM RXN project.

Thanks to the maintainers and contributors of these projects.

## Features

- **SMILES canonicalization** — normalize molecule SMILES strings to canonical forms
- **InChI-based canonicalization** (`inchify`) — canonicalize via InChI round-trip for tautomer-insensitive matching
- **Reaction processing** (`process_reaction`) — canonicalize reaction SMILES, deduplicate agents, and extract product
- **SMILES tokenizer** — tokenize SMILES strings for sequence models
- **Augmented SMILES** — generate randomized SMILES for data augmentation

## Installation

### From PyPI (recommended)

```bash
pip install rxkit
```

### From source

```bash
git clone https://github.com/CreamyLong/RxKit.git
cd RxKit
pip install -e .
```

**Note**: RDKit is a required dependency. If you encounter RDKit installation issues, install it via conda first:

```bash
conda install -c conda-forge rdkit
pip install rxkit
```

## Quick Start

```python
from rxkit.canonicalize import canonicalize_smiles, inchify, process_reaction, smi_tokenizer

# Canonicalize a single molecule
canonicalize_smiles("CCO")
# 'CCO'

# Tautomer-insensitive canonicalization via InChI
inchify("CC(=O)O")
# 'CC(=O)O'

# Process a reaction
process_reaction("C=C.CC(=O)O>>CC(=O)OC")
# 'C=C.CC(=O)O>>COC(C)=O'

# Tokenize for sequence models
smi_tokenizer("CC(=O)O")
# 'C C ( = O ) O'
```

## API Reference

### `mol_utils`

| Function | Description |
|---|---|
| `canonicalize_smiles(smiles)` | Canonicalize a SMILES string, stripping atom map numbers. Returns empty string on failure. |
| `inchify(smi, extended_tautomer_check=False)` | Canonicalize via InChI round-trip. If the result contains disconnected fragments (`.`) — indicating broken ionic bonds — the original SMILES is returned. Set `extended_tautomer_check=True` for stricter tautomer handling. |
| `smi_tokenizer(smi)` | Tokenize a SMILES string into space-separated tokens. |
| `augm_smile(smi)` | Generate a randomized SMILES string for data augmentation. |
| `getNumHeavyAtoms(smi)` | Return the number of heavy (non-hydrogen) atoms. |
| `get_longest_smiles(smis)` | Given dot-separated SMILES, return the longest fragment (excluding PPh₃). |
| `remove_isotope_information(rxn)` | Strip isotope mass numbers from a reaction SMILES (e.g. `[13CH3]` → `[CH3]`). |
| `is_valid_molecule(smi, allow_empty=False)` | Return `True` if `smi` parses to a non-empty RDKit molecule. Accepts a list of SMILES (all must validate). `allow_empty=True` treats empty strings as valid. |

### `rxn_utils`

| Function | Description |
|---|---|
| `detect_rxn_type(rxn)` | Auto-detect reaction SMILES format. Returns `EXTENDED`, `STANDARD_agent`, `STANDARD_WITH_TILDE`, `STANDARD`, or `UNKNOWN`. |
| `process_reaction(rxn, with_agent=True, rxn_type=None)` | Canonicalize a reaction SMILES. When `with_agent=True`, reagents that also appear in products are removed as agents. `rxn_type` is auto-detected if not provided. The output is always normalised to STANDARD format (`A.B>>D`). |
| `is_valid_reaction(rxn, rxn_type=None, allow_empty_products=False, is_strict=False)` | Return `True` if `rxn` is a structurally well-formed reaction SMILES — non-empty, recognized format, and every fragment is a valid molecule SMILES. Set `allow_empty_products=True` to allow empty product sides (e.g. partial templates). Set `is_strict=True` to additionally reject empty fragments (`..`), duplicate compounds within a side (`CCO.CCO>>...`), and the same molecule appearing on multiple sides (`CCO>>CCO`, `CCO>Cc1ccccc1>...`). Solvent-only templates (`CCO>CCO>CCO`) remain valid. |
| `tilde_to_standard(rxn)` | Replace every `~` in a tilde reaction SMILES with `.` (e.g. `CC.O.[Na+]~[Cl-]>>CCO` → `CC.O.[Na+].[Cl-]>>CCO`). |
| `standard_to_tilde(rxn)` | Naive string-level rewrite: replace `.` with `~` inside each side of a STANDARD reaction SMILES. |
| `parse_tilde_reaction_smiles(rxn)` | Parse a `STANDARD_WITH_TILDE` reaction SMILES into `(reactants, agents, products)` lists of compounds. |
| `tilde_to_extended(rxn)` | Convert a `STANDARD_WITH_TILDE` reaction SMILES to `EXTENDED` form. |
| `parse_extended_reaction_smiles(rxn)` | Parse an `EXTENDED` reaction SMILES into `(reactants, agents, products)` lists of compounds, merging fragments that share a group index back into a single dot-separated compound. |
| `to_extended_reaction_smiles(reactants, agents, products)` | Build an `EXTENDED` reaction SMILES from three lists of compounds. Compounds that already contain `.` are recorded in the trailing `\|f:...\|` fragment info. |
| `extended_to_standard(rxn)` | Drop the trailing `\|f:...\|` fragment info from an EXTENDED reaction SMILES. |
| `standard_to_extended(rxn)` | Build an EXTENDED reaction SMILES from a STANDARD reaction SMILES (each `.` becomes a separate fragment; no fragment groups are recorded because they cannot be inferred from the string). |
| `extended_to_tilde(rxn)` | Convert an `EXTENDED` reaction SMILES to `STANDARD_WITH_TILDE` form by joining fragments in each group with `~` (e.g. `CC.O.[Na+].[Cl-]>>CCO \|f:2.3\|` → `CC.O.[Na+]~[Cl-]>>CCO`). |
| `split_smiles_and_fragment_info(rxn)` | Split an (extended) reaction SMILES into the pure SMILES part and the fragment info suffix. |
| `determine_fragment_groups(extended_info)` | Parse a fragment info string like `\|f:0.2,5.6\|` into `[[0, 2], [5, 6]]`. |

**Supported reaction SMILES types:**

| Type | Example |
|---|---|
| `STANDARD` | `C=C.CC(=O)O>>CC(=O)OC` |
| `STANDARD_agent` | `C=C.CC(=O)O>O>CC(=O)OC` |
| `STANDARD_WITH_TILDE` | `CCO~C>>CCO.CC` |
| `EXTENDED` | `\|f:0.1,2.3\|C=C.CC(=O)O>>CC(=O)OC` |

The output of `process_reaction` is always normalised to `STANDARD` form (`A.B>>D`), with `~` and `|f:...|` markers stripped.

## Extended / Tilde Reaction SMILES

When a compound is a multi-component system (e.g. an ionic solid like NaCl written as `[Na+].[Cl-]`), the EXTENDED format uses a trailing fragment info marker to record which `.`-separated fragments must be grouped into the same compound. In `STANDARD_WITH_TILDE` form the same grouping is expressed inline with `~`.

```python
from rxkit.canonicalize import extended_to_tilde, tilde_to_extended, parse_extended_reaction_smiles

extended = "CC.O.[Na+].[Cl-]>>CCO |f:2.3|"

# EXTENDED -> TILDE
extended_to_tilde(extended)
# 'CC.O.[Na+]~[Cl-]>>CCO'

# TILDE -> EXTENDED
tilde_to_extended("CC.O.[Na+]~[Cl-]>>CCO")
# 'CC.O.[Na+].[Cl-]>>CCO |f:2.3|'

# Parse to a structured representation
parse_extended_reaction_smiles(extended)
# (['CC', 'O', '[Na+].[Cl-]'], [], ['CCO'])
```

Fragment indices in the `|f:...|` block are **global** (i.e. they index into the dot-separated fragments of the whole reaction SMILES, going through reactants → reagents → products in order).  Multi-fragment compounds are placed **after** single-fragment compounds in each parsed list, matching the convention used by [`rxn-chemutils`](https://github.com/rxn4chemistry/rxn-chemutils).

## Reaction Processing Details

`process_reaction` takes a reaction SMILES in `reactants>reagents>products` format:

- **`with_agent=True`** (default): Reagents are merged with reactants as precursors. Any molecule appearing in both precursors and products is treated as an agent and removed from both sides. This is useful for cleaning up reactions where catalysts/solvents are listed as reagents.
- **`with_agent=False`**: Reactants and products are canonicalized independently. Reagents are ignored.

## Run Tests

```bash
pip install pytest
pytest tests/
```

## License

MIT
