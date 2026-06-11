"""Tests for rxn_utils module."""

import pytest
from rxkit.canonicalize import (
    process_reaction,
    detect_rxn_type,
    STANDARD,
    STANDARD_AGENT,
    STANDARD_WITH_TILDE,
    EXTENDED,
    UNKNOWN,
)


class TestDetectRxnType:
    def test_standard(self):
        assert detect_rxn_type("C=C.CC(=O)O>>CC(=O)OC") == STANDARD

    def test_standard_agent(self):
        assert detect_rxn_type("C=C.CC(=O)O>O>CC(=O)OC") == STANDARD_AGENT

    def test_standard_with_tilde(self):
        assert detect_rxn_type("C~CC(=O)O>C~CC(=O)O") == STANDARD_WITH_TILDE

    def test_extended(self):
        assert detect_rxn_type("|f:0.2,3.4.5|C=C.CC(=O)O>>CC(=O)OC") == EXTENDED

    def test_unknown(self):
        assert detect_rxn_type("not_a_reaction") == UNKNOWN

    def test_empty(self):
        assert detect_rxn_type("") == UNKNOWN

    def test_none(self):
        assert detect_rxn_type(None) == UNKNOWN

    def test_extended_priority_over_tilde(self):
        """If both markers appear, EXTENDED wins (|...| is more specific)."""
        assert detect_rxn_type("|f:0|C~CC>>D") == EXTENDED


class TestProcessReactionStandard:
    def test_simple_reaction(self):
        """STANDARD format: A.B>>D"""
        result = process_reaction("C=C.CC(=O)O>>CC(=O)OC")
        assert ">>" in result
        assert result != ""

    def test_simple_reaction_without_agent(self):
        result = process_reaction("C=C.CC(=O)O>>CC(=O)OC", with_agent=False)
        assert ">>" in result
        assert result != ""

    def test_no_reagents(self):
        result = process_reaction("C=C.CC(=O)O>>CC(=O)OC")
        assert result != ""

    def test_output_format(self):
        result = process_reaction("C=C.CC(=O)O>>CC(=O)OC", with_agent=False)
        parts = result.split(">>")
        assert len(parts) == 2

    def test_amide_formation(self):
        result = process_reaction("CC(=O)O.CN>>CNC(C)=O.O", with_agent=False)
        assert result != ""


class TestProcessReactionStandardAgent:
    def test_with_explicit_agent(self):
        result = process_reaction("C=C.CC(=O)O>O>CC(=O)OC")
        assert ">>" in result
        assert result != ""

    def test_agent_dedup(self):
        """Catalyst appears in both reagents and products → removed from both."""
        result = process_reaction("c1ccccc1.CCO>CC(=O)O>c1ccccc1")
        assert result != ""
        # Catalyst benzene removed from products side
        products = result.split(">>")[1]
        assert "c1ccccc1" != products

    def test_explicit_type(self):
        result = process_reaction(
            "C=C.CC(=O)O>O>CC(=O)OC", rxn_type=STANDARD_AGENT
        )
        assert result != ""


class TestProcessReactionTilde:
    def test_tilde_fragments(self):
        """Tilde-separated alternatives on each side."""
        result = process_reaction("CCO~C>>CCO.CC")
        assert ">>" in result
        assert result != ""

    def test_tilde_expanded_to_dots(self):
        """Output uses dots, not tildes, after processing."""
        result = process_reaction("CCO~C>>CCO")
        assert "~" not in result
        assert ">>" in result

    def test_explicit_type(self):
        result = process_reaction("CCO~C>>CCO", rxn_type=STANDARD_WITH_TILDE)
        assert result != ""


class TestProcessReactionExtended:
    def test_fragment_markers_stripped(self):
        """EXTENDED |f:...| markers are stripped before processing."""
        result = process_reaction("|f:0.2,3.4.5|C=C.CC(=O)O>>CC(=O)OC")
        assert "|" not in result
        assert ">>" in result

    def test_fragment_with_agent(self):
        result = process_reaction("|f:0.1|C=C.CC(=O)O>O>CC(=O)OC")
        assert "|" not in result
        assert ">>" in result

    def test_explicit_type(self):
        result = process_reaction(
            "|f:0.1|C=C.CC(=O)O>>CC(=O)OC", rxn_type=EXTENDED
        )
        assert result != ""


class TestProcessReactionErrorHandling:
    def test_invalid_input(self):
        result = process_reaction("invalid>stuff>more")
        assert result == ""

    def test_empty_input(self):
        result = process_reaction("")
        assert result == ""

    def test_unknown_type_returns_empty(self):
        result = process_reaction("garbage_no_delimiters")
        assert result == ""


class TestProcessReactionOutputNormalization:
    """Output should always be STANDARD format: A.B>>D"""

    @pytest.mark.parametrize("rxn", [
        "C=C.CC(=O)O>>CC(=O)OC",                                   # STANDARD
        "C=C.CC(=O)O>O>CC(=O)OC",                                  # STANDARD_agent
        "CCO~C>>CCO.CC",                                           # TILDE
        "|f:0.1,2.3|C=C.CC(=O)O>>CC(=O)OC",                       # EXTENDED
    ])
    def test_output_has_standard_separator(self, rxn):
        result = process_reaction(rxn)
        # No tildes, no fragment markers
        assert "~" not in result
        assert "|" not in result
        # Has the >> separator and exactly 2 parts when split
        assert result.count(">>") == 1


class TestIsValidReaction:
    """Tests for the is_valid_reaction helper."""

    def test_valid_standard(self):
        from rxkit.canonicalize import is_valid_reaction
        assert is_valid_reaction("CCO.CCO>>CCOCCO") is True

    def test_valid_standard_agent(self):
        from rxkit.canonicalize import is_valid_reaction
        assert is_valid_reaction("C=C.CC(=O)O>O>CC(=O)OC") is True

    def test_valid_tilde(self):
        from rxkit.canonicalize import is_valid_reaction
        assert is_valid_reaction("CC.O.[Na+]~[Cl-]>>CCO") is True

    def test_valid_extended(self):
        from rxkit.canonicalize import is_valid_reaction
        assert is_valid_reaction("CC.O.[Na+].[Cl-]>>CCO |f:2.3|") is True

    def test_invalid_garbage(self):
        from rxkit.canonicalize import is_valid_reaction
        assert is_valid_reaction("not_a_reaction") is False

    def test_invalid_garbage_fragments(self):
        from rxkit.canonicalize import is_valid_reaction
        assert is_valid_reaction("CCO.garbage>>CCOCCO") is False

    def test_empty_is_invalid(self):
        from rxkit.canonicalize import is_valid_reaction
        assert is_valid_reaction("") is False
        assert is_valid_reaction("   ") is False
        assert is_valid_reaction(None) is False  # type: ignore[arg-type]

    def test_empty_products_invalid_by_default(self):
        from rxkit.canonicalize import is_valid_reaction
        assert is_valid_reaction("CCO>>") is False

    def test_empty_products_allowed_when_flag(self):
        from rxkit.canonicalize import is_valid_reaction
        # Empty products are allowed when the flag is explicitly set
        assert is_valid_reaction("CCO>>", allow_empty_products=True) is True
        # Same for STANDARD_agent-style templates with empty products
        assert is_valid_reaction("CCO>O>", allow_empty_products=True) is True

    def test_empty_reactants_invalid(self):
        from rxkit.canonicalize import is_valid_reaction
        assert is_valid_reaction(">>CCO") is False

    def test_explicit_rxn_type(self):
        from rxkit.canonicalize import is_valid_reaction, STANDARD
        # Pass an EXTENDED SMILES but tell the validator it is STANDARD.
        # The result should be False because the |f:...| marker is not a
        # valid molecule SMILES.
        rxn = "CCO.CCO>>CCOCCO |f:0.1|"
        assert is_valid_reaction(rxn, rxn_type=STANDARD) is False
        # Same string with the correct auto-detected type is valid:
        assert is_valid_reaction(rxn) is True

    def test_too_many_sides(self):
        from rxkit.canonicalize import is_valid_reaction
        # 4 '>' characters yields 5 sides
        assert is_valid_reaction("CCO>CCO>CCO>CCO>CCO") is False


class TestIsValidReactionStrict:
    """Tests for the ``is_strict`` mode of is_valid_reaction."""

    def test_duplicate_in_reactants_raises(self):
        from rxkit.canonicalize import is_valid_reaction
        # Lenient mode: passes.  Strict mode: rejects duplicate reactant.
        assert is_valid_reaction("CCO.CCO>>CCOCCO") is True
        assert is_valid_reaction("CCO.CCO>>CCOCCO", is_strict=True) is False

    def test_duplicate_in_products_raises(self):
        from rxkit.canonicalize import is_valid_reaction
        assert is_valid_reaction("CCOCCO>>CCO.CCO") is True
        assert is_valid_reaction("CCOCCO>>CCO.CCO", is_strict=True) is False

    def test_duplicate_in_reagents_raises(self):
        from rxkit.canonicalize import is_valid_reaction
        assert is_valid_reaction("CCN>c1ccccc1.c1ccccc1>CCNC(=O)c1ccccc1") is True
        assert (
            is_valid_reaction(
                "CCN>c1ccccc1.c1ccccc1>CCNC(=O)c1ccccc1",
                is_strict=True,
            )
            is False
        )

    def test_same_reactant_and_product_raises(self):
        from rxkit.canonicalize import is_valid_reaction
        # No-op reaction: same molecule on both sides.
        assert is_valid_reaction("CCO>>CCO") is True
        assert is_valid_reaction("CCO>>CCO", is_strict=True) is False

    def test_reactant_in_products_raises(self):
        from rxkit.canonicalize import is_valid_reaction
        # ``CCO`` appears as both reactant and product.  Not a pure
        # solvent-only reaction, so strict mode rejects it.
        assert is_valid_reaction("CCO.CCN>>CCO.CCN") is True
        assert is_valid_reaction("CCO.CCN>>CCO.CCN", is_strict=True) is False

    def test_reagent_appears_in_reactants_raises(self):
        from rxkit.canonicalize import is_valid_reaction
        # Benzene appears as both reagent and reactant.
        assert is_valid_reaction("c1ccccc1.CC(=O)O>c1ccccc1>CC(=O)OC") is True
        assert (
            is_valid_reaction(
                "c1ccccc1.CC(=O)O>c1ccccc1>CC(=O)OC",
                is_strict=True,
            )
            is False
        )

    def test_solvent_all_three_sides_is_allowed(self):
        from rxkit.canonicalize import is_valid_reaction
        # ``CCO`` is a solvent appearing on all three sides.  This is
        # the canonical "solvent-only" template and must be allowed
        # even in strict mode.
        assert is_valid_reaction("CCO>CCO>CCO", is_strict=True) is True
        assert (
            is_valid_reaction("c1ccccc1>c1ccccc1>c1ccccc1", is_strict=True) is True
        )

    def test_empty_fragment_raises(self):
        from rxkit.canonicalize import is_valid_reaction
        # Consecutive dots
        assert is_valid_reaction("CCO..CCN>>CCNCCO") is False
        assert is_valid_reaction("CCO..CCN>>CCNCCO", is_strict=True) is False

    def test_leading_dot_fragment_raises(self):
        from rxkit.canonicalize import is_valid_reaction
        # Leading dot — also an empty fragment.
        assert is_valid_reaction(".CCO>>CCO") is False
        assert is_valid_reaction(".CCO>>CCO", is_strict=True) is False

    def test_tilde_reaction_with_duplicate(self):
        from rxkit.canonicalize import is_valid_reaction
        # Tilde reaction where two fragment groups are actually the
        # same compound.
        assert is_valid_reaction("CCO~CCO>>CCOCCO") is True
        assert is_valid_reaction("CCO~CCO>>CCOCCO", is_strict=True) is False

    def test_extended_reaction_with_duplicate(self):
        from rxkit.canonicalize import is_valid_reaction
        # Same molecule in two different groups — also a duplicate.
        assert is_valid_reaction("CCO.CCO>>CCOCCO |f:0.1|") is True
        assert (
            is_valid_reaction("CCO.CCO>>CCOCCO |f:0.1|", is_strict=True) is False
        )

    def test_strict_passes_for_clean_reaction(self):
        from rxkit.canonicalize import is_valid_reaction
        # A clean, well-formed reaction should still pass in strict mode.
        assert is_valid_reaction("CCO.CCN>>CCNCCO", is_strict=True) is True
        assert (
            is_valid_reaction("c1ccccc1.CC(=O)O>O>CC(=O)OC", is_strict=True) is True
        )

    def test_strict_canonicalisation_isolates_different_writings(self):
        from rxkit.canonicalize import is_valid_reaction
        # ``CCO`` and ``OCC`` are the same molecule written differently.
        assert is_valid_reaction("CCO.OCC>>CCOCCO") is True
        assert is_valid_reaction("CCO.OCC>>CCOCCO", is_strict=True) is False

