"""Tests for mol_utils module."""

import pytest
from rxkit.canonicalize import canonicalize_smiles, inchify, smi_tokenizer, augm_smile, getNumHeavyAtoms, get_longest_smiles


class TestCanonicalizeSmiles:
    def test_simple_molecule(self):
        assert canonicalize_smiles("CCO") == "CCO"

    def test_isomeric_smiles(self):
        result = canonicalize_smiles("C/C=C/C")
        assert result in ("C/C=C/C", "CC=CC")

    def test_strips_atom_map(self):
        result = canonicalize_smiles("[CH3:1][OH:2]")
        assert result == "CO"

    def test_invalid_smiles(self):
        assert canonicalize_smiles("invalid") == ""

    def test_benzene(self):
        assert canonicalize_smiles("c1ccccc1") == "c1ccccc1"

    def test_same_molecule_different_representation(self):
        a = canonicalize_smiles("C(=O)O")
        b = canonicalize_smiles("OC=O")
        c = canonicalize_smiles("O=C(O)")
        assert a == b == c


class TestInchify:
    def test_simple(self):
        result = inchify("CCO")
        assert result == "CCO"

    def test_tautomer(self):
        """Tautomers are canonicalized via InChI."""
        keto = inchify("CC(=O)C")
        enol = inchify("CC(O)=C")
        # keto and enol forms are distinct; inchify returns canonical form for each
        assert keto == "CC(C)=O"
        assert enol == "C=C(C)O"

    def test_ion_bond_preserved(self):
        """Ionic bond SMILES should be preserved as-is."""
        result = inchify("[Na+].[Cl-]")
        assert "Na+" in result and "Cl-" in result
        assert result != ""

    def test_extended_tautomer(self):
        result = inchify("CC(=O)O", extended_tautomer_check=True)
        assert result == "CC(=O)O"


class TestSmiTokenizer:
    def test_simple(self):
        tokens = smi_tokenizer("CCO")
        assert tokens == "C C O"

    def test_with_brackets(self):
        tokens = smi_tokenizer("[Na+].[Cl-]")
        assert tokens == "[Na+] . [Cl-]"

    def test_with_branches(self):
        tokens = smi_tokenizer("CC(=O)O")
        assert tokens == "C C ( = O ) O"

    def test_roundtrip(self):
        smi = "CC(=O)O"
        tokens = smi_tokenizer(smi)
        reconstructed = tokens.replace(" ", "")
        assert reconstructed == smi

    def test_with_reaction(self):
        tokens = smi_tokenizer("C=C.CC(O)=O>>CC(=O)OC")
        assert ">>" in tokens


class TestAugmSmile:
    def test_returns_valid_smiles(self):
        smi = "CC(=O)O"
        augmented = augm_smile(smi)
        mol = __import__('rdkit').Chem.MolFromSmiles(augmented)
        assert mol is not None

    def test_benzene_augmentation(self):
        smi = "c1ccccc1"
        augmented = augm_smile(smi)
        mol = __import__('rdkit').Chem.MolFromSmiles(augmented)
        assert mol is not None


class TestGetNumHeavyAtoms:
    def test_ethanol(self):
        assert getNumHeavyAtoms("CCO") == 3

    def test_benzene(self):
        assert getNumHeavyAtoms("c1ccccc1") == 6

    def test_water(self):
        assert getNumHeavyAtoms("O") == 1


class TestGetLongestSmiles:
    def test_single_molecule(self):
        assert get_longest_smiles("CCO") == "CCO"

    def test_multi_molecule(self):
        result = get_longest_smiles("C.CCO.CC")
        assert result == "CCO"

    def test_pph3_removed(self):
        result = get_longest_smiles("c1ccc(cc1)P(c1ccccc1)c1ccccc1.CCO")
        assert result == "CCO"

    def test_empty(self):
        assert get_longest_smiles("") == ""

    def test_all_same_length(self):
        """When all fragments have same length, returns one of them."""
        result = get_longest_smiles("C.C")
        assert result == "C"


class TestTildeStandardConversion:
    """Tests for tilde <-> standard reaction SMILES conversion."""

    def test_tilde_to_standard(self):
        from canonicalize import tilde_to_standard
        assert tilde_to_standard("CC.O.[Na+]~[Cl-]>>CCO") == "CC.O.[Na+].[Cl-]>>CCO"

    def test_standard_to_tilde_naive(self):
        """``standard_to_tilde`` is a string-level rewrite: ``.`` -> ``~`` everywhere on each side."""
        from canonicalize import standard_to_tilde
        assert standard_to_tilde("CC.O.[Na+].[Cl-]>>CCO") == "CC~O~[Na+]~[Cl-]>>CCO"

    def test_no_tilde_unchanged(self):
        from canonicalize import tilde_to_standard
        assert tilde_to_standard("CC.O>>CCO") == "CC.O>>CCO"


class TestExtendedTildeConversion:
    """The user requirement: EXTENDED -> STANDARD_WITH_TILDE conversion.

    Reference: rxn-chemutils' extended_reaction_smiles / reaction_smiles
    modules use ``|f:idx1.idx2,...|`` to record which fragments belong to
    the same compound.  In STANDARD_WITH_TILDE form, those fragments are
    joined with ``~`` instead.
    """

    def test_split_smiles_and_fragment_info(self):
        from canonicalize import split_smiles_and_fragment_info
        pure, info = split_smiles_and_fragment_info(
            "CC.O.[Na+].[Cl-]>>CCO |f:2.3|"
        )
        assert pure == "CC.O.[Na+].[Cl-]>>CCO"
        assert info == "|f:2.3|"

    def test_determine_fragment_groups(self):
        from canonicalize import determine_fragment_groups
        assert determine_fragment_groups("|f:2.3|") == [[2, 3]]
        assert determine_fragment_groups("|f:0.2,5.6|") == [[0, 2], [5, 6]]
        assert determine_fragment_groups("") == []

    def test_user_requirement_extended_to_tilde(self):
        """The exact example from the user's request."""
        from canonicalize import extended_to_tilde
        src = "CC.O.[Na+].[Cl-]>>CCO |f:2.3|"
        expected = "CC.O.[Na+]~[Cl-]>>CCO"
        assert extended_to_tilde(src) == expected

    def test_roundtrip_tilde_to_extended(self):
        from canonicalize import tilde_to_extended
        result = tilde_to_extended("CC.O.[Na+]~[Cl-]>>CCO")
        assert result == "CC.O.[Na+].[Cl-]>>CCO |f:2.3|"

    def test_parse_extended(self):
        from canonicalize import parse_extended_reaction_smiles
        r, a, p = parse_extended_reaction_smiles(
            "CC.O.[Na+].[Cl-]>>CCO |f:2.3|"
        )
        # singletons first, multi-fragment compounds last
        assert r == ["CC", "O", "[Na+].[Cl-]"]
        assert a == []
        assert p == ["CCO"]

    def test_parse_extended_groups_on_both_sides(self):
        from canonicalize import parse_extended_reaction_smiles
        r, a, p = parse_extended_reaction_smiles(
            "CC.O.[Na+].[Cl-]>>CCO.[Na+].[Cl-] |f:2.3,5.6|"
        )
        assert r == ["CC", "O", "[Na+].[Cl-]"]
        assert p == ["CCO", "[Na+].[Cl-]"]

    def test_to_extended_reaction_smiles(self):
        from canonicalize import to_extended_reaction_smiles
        result = to_extended_reaction_smiles(
            ["CC", "O", "[Na+].[Cl-]"], [], ["CCO"]
        )
        assert result == "CC.O.[Na+].[Cl-]>>CCO |f:2.3|"

    def test_to_extended_no_groups_returns_standard(self):
        """If no compound has multiple fragments, the EXTENDED form is just STANDARD."""
        from canonicalize import to_extended_reaction_smiles
        result = to_extended_reaction_smiles(["CC", "O"], [], ["CCO"])
        assert result == "CC.O>>CCO"

    def test_parse_tilde_reaction_smiles(self):
        from canonicalize import parse_tilde_reaction_smiles
        r, a, p = parse_tilde_reaction_smiles("CC.O.[Na+]~[Cl-]>>CCO")
        assert r == ["CC", "O", "[Na+].[Cl-]"]
        assert a == []
        assert p == ["CCO"]

    def test_extended_to_standard_drops_fragment_info(self):
        from canonicalize import extended_to_standard
        assert extended_to_standard(
            "CC.O.[Na+].[Cl-]>>CCO |f:2.3|"
        ) == "CC.O.[Na+].[Cl-]>>CCO"


class TestIsValidMolecule:
    """Tests for the is_valid_molecule helper."""

    def test_valid_simple(self):
        from canonicalize import is_valid_molecule
        assert is_valid_molecule("CCO") is True

    def test_valid_aromatic(self):
        from canonicalize import is_valid_molecule
        assert is_valid_molecule("c1ccccc1") is True

    def test_valid_ionic(self):
        from canonicalize import is_valid_molecule
        assert is_valid_molecule("[Na+].[Cl-]") is True

    def test_invalid_garbage(self):
        from canonicalize import is_valid_molecule
        assert is_valid_molecule("not_a_smiles") is False

    def test_invalid_unbalanced_bracket(self):
        from canonicalize import is_valid_molecule
        assert is_valid_molecule("C[C") is False

    def test_empty_invalid_by_default(self):
        from canonicalize import is_valid_molecule
        assert is_valid_molecule("") is False

    def test_empty_allowed_when_flag_set(self):
        from canonicalize import is_valid_molecule
        assert is_valid_molecule("", allow_empty=True) is True

    def test_none_is_invalid(self):
        from canonicalize import is_valid_molecule
        assert is_valid_molecule(None) is False  # type: ignore[arg-type]

    def test_non_string_is_invalid(self):
        from canonicalize import is_valid_molecule
        assert is_valid_molecule(123) is False  # type: ignore[arg-type]

    def test_list_all_valid(self):
        from canonicalize import is_valid_molecule
        assert is_valid_molecule(["CCO", "c1ccccc1", "[Na+]"]) is True

    def test_list_one_invalid(self):
        from canonicalize import is_valid_molecule
        assert is_valid_molecule(["CCO", "not_a_smiles"]) is False

    def test_list_with_empty_invalid_by_default(self):
        from canonicalize import is_valid_molecule
        assert is_valid_molecule(["CCO", ""]) is False

    def test_list_with_empty_allowed(self):
        from canonicalize import is_valid_molecule
        assert is_valid_molecule(["CCO", ""], allow_empty=True) is True

