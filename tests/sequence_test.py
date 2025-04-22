import unittest

from sequal.modification import Modification
from sequal.sequence import ModdedSequenceGenerator, Sequence

nsequon = Modification(
    "HexNAc", regex_pattern="N[^P][S|T]", mod_type="variable", labile=True
)
osequon = Modification(
    "Mannose", regex_pattern="[S|T]", mod_type="variable", labile=True
)
sulfation = Modification(
    "Sulfation", regex_pattern="S", mod_type="variable", labile=True
)
carbox = Modification(
    "Carboxylation", regex_pattern="E", mod_type="variable", labile=True
)
carbox2 = Modification(
    "Carboxylation2", regex_pattern="E", mod_type="variable", labile=True, mass=43.98983
)
propiona = Modification("Propionamide", regex_pattern="C", mod_type="static")


class TestAASequence(unittest.TestCase):
    def test_normal_sequence(self):
        seq = Sequence("TESTEST")

    def test_mod_rightseq(self):
        seq = Sequence("TEN[HexNAc]ST")

    def test_two_mod_rightseq(self):
        seq = Sequence("TEN[HexNAc][HexNAc]ST")

    def test_mod_leftseq(self):
        seq = Sequence("TE[HexNAc]NST", mod_position="left")
        for i in seq.seq:
            print(i, i.mods)

    def test_two_mod_leftseq(self):
        seq = Sequence("TE[HexNAc][HexNAc]NST", mod_position="left")
        for i in seq.seq:
            print(i, i.mods)

    def test_custom_string(self):
        seq = Sequence("TENST")
        a = {1: "tes", 2: ["1", "200"]}
        print(
            seq.to_string_customize(
                a,
                individual_annotation_enclose=False,
                individual_annotation_separator=".",
            )
        )


class TestModdedSequence(unittest.TestCase):
    def test_variable_mod_generator(self):
        seq = "TESNSTT"
        mods = [nsequon, osequon, carbox]
        g = ModdedSequenceGenerator(seq, mods, [])
        print(g.variable_map.mod_position_dict)
        for i in g.generate():
            print(i)

    def test_static_mod_generator(self):
        seq = "TECSNTT"
        mods = [propiona]
        g = ModdedSequenceGenerator(seq, static_mods=mods)
        for i in g.generate():
            print(i)

    def test_static_and_variable_mod_generator(self):
        seq = "TECSNTT"
        static_mods = [propiona]
        variable_mods = [nsequon, osequon, carbox, carbox2]
        g = ModdedSequenceGenerator(seq, variable_mods, static_mods)
        for i in g.generate():
            print(i)


class TestProForma(unittest.TestCase):
    def test_basic_peptide_with_modification(self):
        """Test a simple peptide with a modification."""
        proforma = "PEP[Phospho]TIDE"
        seq = Sequence.from_proforma(proforma)

        # Check sequence and modification
        assert seq.to_stripped_string() == "PEPTIDE"
        assert seq[2].mods[0].value == "Phospho"

        # Check roundtrip conversion
        assert seq.to_proforma() == proforma

    def test_mass_shift_notation(self):
        """Test mass shift notation."""
        proforma = "PEP[+79.966]TIDE"
        seq = Sequence.from_proforma(proforma)

        # Check sequence and modification
        assert seq.to_stripped_string() == "PEPTIDE"
        assert seq[2].mods[0].value == "+79.966"
        assert abs(seq[2].mods[0].mass - 79.966) < 0.0001

        # Check roundtrip (note: may not be identical due to internal representation)
        assert seq.to_proforma() == "PEP[+79.966]TIDE"

    def test_multiple_modifications(self):
        """Test multiple modifications on a single residue."""
        proforma = "PEPS[Phospho][Acetyl]TIDE"
        seq = Sequence.from_proforma(proforma)

        # Check sequence and modifications
        assert seq.to_stripped_string() == "PEPSTIDE"
        assert len(seq[3].mods) == 2
        assert seq[3].mods[0].value == "Phospho"
        assert seq[3].mods[1].value == "Acetyl"

        # Check roundtrip
        assert seq.to_proforma() == proforma

    def test_terminal_modifications(self):
        """Test N-terminal and C-terminal modifications."""
        proforma = "[Acetyl]-PEPTIDE-[Amidated]"
        seq = Sequence.from_proforma(proforma)

        # Check sequence
        assert seq.to_stripped_string() == "PEPTIDE"

        # Check N-terminal modification
        assert -1 in seq.mods
        assert seq.mods[-1][0].value == "Acetyl"

        # Check C-terminal modification
        assert -2 in seq.mods
        assert seq.mods[-2][0].value == "Amidated"

        # Check roundtrip
        assert seq.to_proforma() == proforma

    def test_ambiguous_modifications(self):
        """Test ambiguous modification locations."""
        proforma = "PEPS{Phospho}TIDE"
        seq = Sequence.from_proforma(proforma)

        # Check sequence
        assert seq.to_stripped_string() == "PEPSTIDE"

        # Check ambiguous modification at position 3
        assert seq[3].mods[0].value == "Phospho"
        assert seq[3].mods[0].mod_type == "ambiguous"

        # Check roundtrip
        assert seq.to_proforma() == proforma

    def test_complex_sequence(self):
        """Test a complex sequence with multiple features."""
        proforma = "[Acetyl]-PEP[Phospho]T{Oxidation}IDE-[Amidated]"
        seq = Sequence.from_proforma(proforma)

        # Check sequence
        assert seq.to_stripped_string() == "PEPTIDE"

        # Check N-terminal modification
        assert -1 in seq.mods
        assert seq.mods[-1][0].value == "Acetyl"

        # Check residue modifications
        assert seq[2].mods[0].value == "Phospho"
        assert seq[3].mods[0].value == "Oxidation"
        assert seq[3].mods[0].mod_type == "ambiguous"

        # Check C-terminal modification
        assert -2 in seq.mods
        assert seq.mods[-2][0].value == "Amidated"

        # Check roundtrip
        assert seq.to_proforma() == proforma

    def test_negative_mass_shift(self):
        """Test negative mass shift notation."""
        proforma = "PEP[-17.027]TIDE"
        seq = Sequence.from_proforma(proforma)

        assert seq.to_stripped_string() == "PEPTIDE"
        assert seq[2].mods[0].value == "-17.027"
        assert abs(seq[2].mods[0].mass + 17.027) < 0.0001

    def test_conversion_from_sequence_to_proforma(self):
        """Test converting a pre-built sequence to ProForma."""
        # Create a sequence manually
        seq = Sequence("PEPTIDE")
        seq[2].add_modification(Modification("Phospho"))

        # Convert to ProForma
        proforma = seq.to_proforma()
        assert proforma == "PEP[Phospho]TIDE"

        # Round trip
        seq2 = Sequence.from_proforma(proforma)
        assert seq2.to_stripped_string() == "PEPTIDE"
        assert seq2[2].mods[0].value == "Phospho"

    def test_inter_chain_crosslinks(self):
        """Test inter-chain crosslinks with // separator."""
        proforma = "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[#XL1]SESPEK"
        seq = Sequence.from_proforma(proforma)

        # Check multi-chain properties
        assert seq.is_multi_chain
        assert len(seq.chains) == 2

        # Check first chain
        assert seq.chains[0].to_stripped_string() == "SEKUENCE"
        assert seq.chains[0].seq[2].mods[0].value == "02001"
        assert seq.chains[0].seq[2].mods[0].crosslink_id == "XL1"

        # Check second chain
        print(seq.chains[1])
        assert seq.chains[1].to_stripped_string() == "EMEVTKSESPEK"
        assert seq.chains[1].seq[5].mods[0].is_crosslink_ref
        assert seq.chains[1].seq[5].mods[0].crosslink_id == "XL1"

        # Verify roundtrip conversion
        assert seq.to_proforma() == proforma

    def test_disulfide_bonds(self):
        """Test disulfide bond representation."""
        proforma = "EVTSEKC[XLMOD:00034#XL1]LEMSC[#XL1]EFD"
        seq = Sequence.from_proforma(proforma)

        # Check sequence
        assert seq.to_stripped_string() != "EVTSEKCLEMMSCEFF"

        # Check disulfide bond
        assert seq.seq[6].mods[0].value == "00034"
        assert seq.seq[6].mods[0].source == "XLMOD"
        assert seq.seq[6].mods[0].crosslink_id == "XL1"
        assert seq.seq[11].mods[0].is_crosslink_ref

        # Verify roundtrip
        assert seq.to_proforma() == proforma

    def test_branched_peptides(self):
        """Test branched peptide representation."""
        # Test basic branched peptide
        proforma = "ETFGD[MOD:00093#BRANCH]LEMSEFD"
        seq = Sequence.from_proforma(proforma)

        # Check sequence
        assert seq.to_stripped_string() == "ETFGDLEMSEFD"

        # Check branch modification - verify source and value separately
        assert seq.seq[4].mods[0].source == "MOD"
        assert seq.seq[4].mods[0].value == "00093"
        assert seq.seq[4].mods[0]._is_branch is True

        # Verify roundtrip
        assert seq.to_proforma() == proforma

        # Test branched with reference
        proforma2 = "ETFGD[MOD:00093#BRANCH]LEMS[#BRANCH]EFD"
        seq2 = Sequence.from_proforma(proforma2)

        # Check branch reference
        assert seq2.seq[8].mods[0]._is_branch_ref is True

        # Verify roundtrip (ensure no duplicate branch references)
        assert seq2.to_proforma() == proforma2

    def test_gno_notation(self):
        proforma = "NEEYN[GNO:G59626AS]K"
        seq = Sequence.from_proforma(proforma)
        assert seq.to_stripped_string() == "NEEYNK"
        assert seq.seq[4].mods[0].mod_value.primary_value == "G59626AS"
        assert seq.seq[4].mods[0].mod_value.source == "GNO"
        assert seq.seq[4].mods[0].mod_value.pipe_values[0].is_valid_glycan

    def test_delta_mass_notation(self):
        """Test parsing and serialization of standard delta mass notations."""
        proforma_strings = [
            "EM[U:+15.9949]EVEES[U:+79.9663]PEK",
            "EM[U:+15.995]EVEES[U:+79.966]PEK",
            "EM[M:+15.9949]EVEES[R:+79.9663]PEK",
            "EM[Obs:+15.995]EVEES[Obs:+79.978]PEK",
        ]

        for proforma in proforma_strings:
            seq = Sequence.from_proforma(proforma)

            # Verify the sequence itself
            assert seq.to_stripped_string() == "EMEVEESPEK"

            # Check the first modification (position 1, M residue)
            mod1 = seq.seq[1].mods[0]
            assert mod1.source in ["U", "M", "X", "Obs"]
            assert mod1.value.startswith("+")
            assert float(mod1.value) > 15.9

            # Check the second modification (position 6, S residue)
            mod2 = seq.seq[6].mods[0]
            assert mod2.source in ["U", "R", "G", "Obs"]
            assert mod2.value.startswith("+")
            assert float(mod2.value) > 79.9

            # Verify round-trip conversion
            assert seq.to_proforma() == proforma

        # Test with "Obs" prefix and specific mass checking
        observed_proforma = "EM[Obs:+15.995]EVEES[Obs:+79.978]PEK"
        seq = Sequence.from_proforma(observed_proforma)

        # Check specific observed value with exact precision
        mod1 = seq.seq[1].mods[0]
        assert mod1.source == "Obs"
        assert mod1.value == "+15.995"

        mod2 = seq.seq[6].mods[0]
        assert mod2.source == "Obs"
        assert mod2.value == "+79.978"

        # Test mixed standard and observed modifications
        mixed_proforma = "EM[U:+15.9949]EVEES[Obs:+79.978]PEK"
        seq = Sequence.from_proforma(mixed_proforma)

        mod1 = seq.seq[1].mods[0]
        assert mod1.source == "U"
        assert mod1.value.startswith("+15.994")

        mod2 = seq.seq[6].mods[0]
        assert mod2.source == "Obs"
        assert mod2.value == "+79.978"

        # Verify round-trip conversion maintains precision
        assert seq.to_proforma() == mixed_proforma

    def test_gap_notation(self):
        """Test gaps of known mass."""
        proforma = "RTAAX[+367.0537]WT"
        seq = Sequence.from_proforma(proforma)

        # Check sequence contains 'X'
        assert seq.to_stripped_string() == "RTAAXWT"

        # Check gap modification
        assert seq.seq[4].value == "X"
        assert seq.seq[4].mods[0].mod_type == "gap"
        assert seq.seq[4].mods[0].value == "+367.0537"

        # Verify mass of the gap
        assert abs(seq.seq[4].mods[0].mass - 367.0537) < 0.0001

        # Verify roundtrip
        assert seq.to_proforma() == proforma

        # Test with a negative mass gap
        proforma2 = "PEPTX[-10.0]IDE"
        seq2 = Sequence.from_proforma(proforma2)
        assert seq2.to_stripped_string() == "PEPTXIDE"
        assert seq2.seq[4].mods[0].value == "-10.0"
        assert seq2.to_proforma() == proforma2

    def test_formula_notation(self):
        """Test modifications with chemical formulas."""
        # Valid formulas
        valid_proforma = [
            "SEQUEN[Formula:C12H20O2]CE",
            "SEQUEN[Formula:C12 H20 O2]CE",  # with spaces
            "SEQUEN[Formula:[13C2]CH6N]CE",  # with isotope
            "SEQUEN[Formula:HN-1O2]CE",  # with negative cardinality
            "SEQUEN[Formula:[13C2][12C-2]H2N]CE",  # complex isotope replacement
        ]

        for proforma in valid_proforma:
            seq = Sequence.from_proforma(proforma)

            # Verify sequence
            assert seq.to_stripped_string() == "SEQUENCE"

            # Check modification
            mod = seq.seq[5].mods[0]
            assert mod.source == "Formula"

            # Check pipe value validation status
            assert mod.mod_value is not None
            for pv in mod.mod_value:
                if pv.source == "Formula":
                    assert pv.is_valid_formula

            # Verify roundtrip conversion
            assert seq.to_proforma() == proforma

        # Invalid formulas
        invalid_proforma = [
            "SEQUEN[Formula:123]CE",  # Not starting with element
            "SEQUEN[Formula:C12H20O2+]CE",  # Invalid character
            "SEQUEN[Formula:C-0]CE",  # Zero cardinality
        ]

        for proforma in invalid_proforma:
            seq = Sequence.from_proforma(proforma)

            # Check modification is preserved but marked invalid
            mod = seq.seq[5].mods[0]
            assert mod.source == "Formula"

            # Check pipe value validation status
            assert mod.mod_value is not None
            for pv in mod.mod_value:
                if pv.source == "Formula":
                    assert (
                        not pv.is_valid_formula
                    ), f"Formula {pv.value} should be invalid"

            # Verify roundtrip conversion still works
            assert seq.to_proforma() == proforma

    def test_glycan_notation(self):
        """Test modifications with glycan notation."""
        proforma_strings = [
            "SEQUEN[Glycan:Hex2HexNAc]CE",
            "SEQUEN[Glycan:HexNAc1Hex2]CE",
            "SEQUEN[Glycan:NeuAc1Hex3]CE",
            "SEQUEN[Glycan:Fuc1dHex2Pen3]CE",
        ]

        for proforma in proforma_strings:
            seq = Sequence.from_proforma(proforma)

            # Verify the sequence itself
            assert seq.to_stripped_string() == "SEQUENCE"

            # Check modification
            mod = seq.seq[5].mods[0]
            assert mod.source == "Glycan"

            # Verify roundtrip conversion
            assert seq.to_proforma() == proforma

    def test_valid_labile_modifications(self):
        # Test single labile modification
        proforma_str = "{Glycan:Hex}EMEVNESPEK"
        seq = Sequence.from_proforma(proforma_str)
        assert seq.to_stripped_string() == "EMEVNESPEK"
        assert len(seq.mods[-3]) == 1
        assert seq.mods[-3][0].value == "Hex"
        assert seq.mods[-3][0].source == "Glycan"
        assert seq.to_proforma() == proforma_str
        # Test multiple labile modifications
        proforma_str = "{Glycan:Hex}{Glycan:NeuAc}EMEVNESPEK"
        seq = Sequence.from_proforma(proforma_str)
        assert seq.to_stripped_string() == "EMEVNESPEK"
        assert len(seq.mods[-3]) == 2
        assert seq.mods[-3][0].value == "Hex"
        assert seq.mods[-3][1].value == "NeuAc"
        assert seq.to_proforma() == proforma_str

    def test_unknown_position_modifications(self):
        # Single unknown position
        seq = Sequence.from_proforma("[Phospho]?EMEVNESPEK")
        assert seq.to_stripped_string() == "EMEVNESPEK"
        assert len(seq.mods[-4]) == 1
        assert seq.mods[-4][0].value == "Phospho"
        # Test roundtrip conversion
        assert seq.to_proforma() == "[Phospho]?EMEVNESPEK"

        # Multiple unknown positions with individual listing
        seq = Sequence.from_proforma("[Phospho][Phospho]?EMEVNESPEK")
        assert seq.to_stripped_string() == "EMEVNESPEK"
        assert len(seq.mods[-4]) == 2
        assert all(mod.value == "Phospho" for mod in seq.mods[-4])
        # Test roundtrip conversion
        assert seq.to_proforma() == "[Phospho]^2?EMEVNESPEK"

        # Multiple unknown positions with caret notation
        seq = Sequence.from_proforma("[Phospho]^2?EMEVNESPEK")
        assert seq.to_stripped_string() == "EMEVNESPEK"
        assert len(seq.mods[-4]) == 2
        assert all(mod.value == "Phospho" for mod in seq.mods[-4])
        # Test roundtrip conversion
        assert seq.to_proforma() == "[Phospho]^2?EMEVNESPEK"

        # With N-terminal modification
        seq = Sequence.from_proforma("[Phospho]^2?[Acetyl]-EMEVNESPEK")
        assert len(seq.mods[-4]) == 2
        assert seq.mods[-1][0].value == "Acetyl"
        # Test roundtrip conversion
        assert seq.to_proforma() == "[Phospho]^2?[Acetyl]-EMEVNESPEK"

    def test_ambiguity_groups(self):
        """Test ambiguity groups for modifications with multiple possible sites."""
        # Test with a simple ambiguity group
        proforma = "EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1]PEK"
        seq = Sequence.from_proforma(proforma)

        # Check sequence
        assert seq.to_stripped_string() == "EMEVTSESPEK"

        # Check modifications
        assert seq.seq[1].mods[0].value == "Oxidation"  # Regular mod

        # Check ambiguity group
        assert seq.seq[4].mods[0].is_ambiguity_ref  # T position
        assert seq.seq[4].mods[0].ambiguity_group == "g1"

        assert seq.seq[5].mods[0].is_ambiguity_ref  # First S position
        assert seq.seq[5].mods[0].ambiguity_group == "g1"

        assert seq.seq[7].mods[0].value == "Phospho"  # Preferred S position
        assert seq.seq[7].mods[0].ambiguity_group == "g1"
        assert not seq.seq[7].mods[0].is_ambiguity_ref

        # Check roundtrip
        assert seq.to_proforma() == proforma

        # Test with multiple ambiguity groups
        proforma2 = "EM[Oxidation#g1]E[#g1]VTS[Phospho#g2]ES[#g2]PEK"
        seq2 = Sequence.from_proforma(proforma2)

        # Check roundtrip with multiple groups
        assert seq2.to_proforma() == proforma2

    def test_range_modifications(self):
        """Test range notation for modifications with multiple possible sites."""
        # Simple range
        proforma = "PRT(ESFRMS)[+19.0523]ISK"
        seq = Sequence.from_proforma(proforma)

        # Check sequence
        assert seq.to_stripped_string() == "PRTESFRMSISK"

        # Check modification is applied to range
        for i in range(3, 9):  # Positions 3-8 (ESFRMS)
            assert any(
                mod.in_range and mod.value == "+19.0523" for mod in seq.seq[i].mods
            )

        # Check roundtrip conversion
        assert seq.to_proforma() == proforma

        # Range containing modification
        proforma2 = "PRT(EC[Carbamidomethyl]FRMS)[+19.0523]ISK"
        seq2 = Sequence.from_proforma(proforma2)

        # Check both modifications
        assert any(mod.value == "Carbamidomethyl" for mod in seq2.seq[4].mods)
        for i in range(3, 9):
            assert any(
                mod.in_range and mod.value == "+19.0523" for mod in seq2.seq[i].mods
            )

        # Check roundtrip
        assert seq2.to_proforma() == proforma2

    def test_localization_scores(self):
        """Test ambiguity groups with localization scores."""
        # Basic example from the spec
        proforma = "EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK"
        seq = Sequence.from_proforma(proforma)

        # Check localization scores
        assert seq.seq[4].mods[0].localization_score == 0.01  # T position
        assert seq.seq[5].mods[0].localization_score == 0.09  # S position
        assert seq.seq[7].mods[0].localization_score == 0.90  # Preferred S position

        # Check roundtrip conversion preserves scores
        assert seq.to_proforma() == proforma

        # Test unknown position with scores example from the spec
        proforma2 = (
            "[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]S[#s1(0.09)]ES[#s1(0.90)]PEK"
        )
        seq2 = Sequence.from_proforma(proforma2)

        # Check unknown position modification
        assert len(seq2.mods[-4]) == 1
        assert seq2.mods[-4][0].value == "Phospho"
        assert seq2.mods[-4][0].ambiguity_group == "s1"

        # Check all localization scores
        assert seq2.seq[4].mods[0].localization_score == 0.01
        assert seq2.seq[5].mods[0].localization_score == 0.09
        assert seq2.seq[7].mods[0].localization_score == 0.90

        # Check roundtrip
        assert seq2.to_proforma() == proforma2

    def test_range_with_localization_scores(self):
        """Test localization scores with range notation."""
        # Example 1: Range with ambiguity group and score
        proforma = "PRT(ESFRMS)[+19.0523#g1(0.01)]ISK[#g1(0.99)]"
        seq = Sequence.from_proforma(proforma)

        # Verify range modification
        range_start = 3
        range_end = 8
        for i in range(range_start, range_end + 1):
            mods = seq.seq[i].mods
            assert len(mods) == 1
            assert mods[0].value == "+19.0523"
            assert mods[0].ambiguity_group == "g1"
            assert mods[0].localization_score == 0.01
            assert mods[0].in_range == True

        # Check last position with reference to the same ambiguity group
        assert seq.seq[11].mods[0].is_ambiguity_ref == True
        assert seq.seq[11].mods[0].ambiguity_group == "g1"
        assert seq.seq[11].mods[0].localization_score == 0.99

        # Test roundtrip
        assert seq.to_proforma() == "PRT(ESFRMS)[+19.0523#g1(0.01)]ISK[#g1(0.99)]"

        # Example 2: More complex case
        proforma2 = "PR[#g1(0.91)]T(EC[Carbamidomethyl]FRMS)[+19.0523#g1(0.09)]ISK"
        seq2 = Sequence.from_proforma(proforma2)

        # Check individual position outside range
        assert seq2.seq[1].mods[0].is_ambiguity_ref == True
        assert seq2.seq[1].mods[0].ambiguity_group == "g1"
        assert seq2.seq[1].mods[0].localization_score == 0.91

        # Check nested modification inside range
        assert seq2.seq[4].mods[0].value == "Carbamidomethyl"

        # Check range modification with score
        for i in range(3, 9):
            if i != 4:  # Skip position with Carbamidomethyl which is handled separately
                mods = [m for m in seq2.seq[i].mods if m.in_range]
                assert len(mods) > 0
                assert mods[0].value == "+19.0523"
                assert mods[0].ambiguity_group == "g1"
                assert mods[0].localization_score == 0.09

        # Test roundtrip
        assert (
            seq2.to_proforma()
            == "PR[#g1(0.91)]T(EC[Carbamidomethyl]FRMS)[+19.0523#g1(0.09)]ISK"
        )

    def test_multiple_modifications_same_residue(self):
        """Test multiple modifications on the same amino acid or range."""
        # Example from the spec with multiple mods on a range
        proforma = "MPGLVDSNPAPPESQEKKPLK(PCCACPETKKARDACIIEKGEEHCGHLIEAHKECMRALGFKI)[Oxidation][Oxidation][half cystine][half cystine]"
        seq = Sequence.from_proforma(proforma)

        # Verify sequence structure
        assert (
            seq.to_stripped_string()
            == "MPGLVDSNPAPPESQEKKPLKPCCACPETKKARDACIIEKGEEHCGHLIEAHKECMRALGFKI"
        )

        # Check that all four modifications are applied to the range
        range_start = 22
        range_end = 62

        # Test that all positions in the range have the four modifications
        for i in range(range_start, range_end + 1):
            assert len(seq.seq[i].mods) == 4

            # Check each modification type
            mod_values = [mod.value for mod in seq.seq[i].mods]
            assert mod_values.count("Oxidation") == 2
            assert mod_values.count("half cystine") == 2

        # Test a simpler case with multiple modifications on a single residue
        proforma2 = "PEPTIDEK[Acetyl][Methyl]"
        seq2 = Sequence.from_proforma(proforma2)

        # Check modifications on K
        mods = seq2.seq[7].mods
        assert len(mods) == 2
        assert "Acetyl" in [m.value for m in mods]
        assert "Methyl" in [m.value for m in mods]

        # Verify roundtrip
        assert seq2.to_proforma() == "PEPTIDEK[Acetyl][Methyl]"

    def test_global_modifications(self):
        """Test global modifications in ProForma format (section 4.6)."""
        # Test Case 1: Isotope labeling
        isotope_tests = [
            ("<13C>ATPEILTVNSIGQLK", "13C", None),
            ("<15N>ATPEILTVNSIGQLK", "15N", None),
            ("<D>ATPEILTVNSIGQLK", "D", None),
            ("<13C><15N>ATPEILTVNSIGQLK", ["13C", "15N"], None),
        ]

        for proforma, expected_value, _ in isotope_tests:
            seq = Sequence.from_proforma(proforma)

            # Check base sequence
            assert seq.to_stripped_string() == "ATPEILTVNSIGQLK"

            # Check global mods
            if isinstance(expected_value, list):
                assert len(seq.global_mods) == len(expected_value)
                for i, mod in enumerate(seq.global_mods):
                    assert mod.value == expected_value[i]
                    assert mod.global_mod_type == "isotope"
                    assert mod.target_residues is None
            else:
                assert len(seq.global_mods) == 1
                assert seq.global_mods[0].value == expected_value
                assert seq.global_mods[0].global_mod_type == "isotope"
                assert seq.global_mods[0].target_residues is None

            # Test roundtrip
            assert seq.to_proforma() == proforma

        # Test Case 2: Fixed protein modifications with target residues
        fixed_mod_tests = [
            (
                "<[S-carboxamidomethyl-L-cysteine]@C>ATPEILTCNSIGCLK",
                "S-carboxamidomethyl-L-cysteine",
                ["C"],
            ),
            ("<[MOD:01090]@C>ATPEILTCNSIGCLK", "01090", ["C"]),
            ("<[Oxidation]@C,M>MTPEILTCNSIGCLK", "Oxidation", ["C", "M"]),
        ]

        for proforma, expected_value, expected_targets in fixed_mod_tests:
            seq = Sequence.from_proforma(proforma)

            # Check base sequence
            assert seq.to_stripped_string() == proforma.split(">")[1]

            # Check global mod
            assert len(seq.global_mods) == 1
            assert seq.global_mods[0].value == expected_value
            assert seq.global_mods[0].global_mod_type == "fixed"
            assert seq.global_mods[0].target_residues == expected_targets

            # Test roundtrip
            assert seq.to_proforma() == proforma

        # Test Case 3: Combined with other modifications
        combined_tests = [
            ("<[MOD:01090]@C>[Phospho]?EM[Oxidation]EVTSECSPEK"),
            ("<[MOD:01090]@C>[Acetyl]-EM[Oxidation]EVTSECSPEK"),
        ]

        for proforma in combined_tests:
            seq = Sequence.from_proforma(proforma)

            # Check global mod exists
            assert len(seq.global_mods) == 1
            assert seq.global_mods[0].source == "MOD"
            assert seq.global_mods[0].value == "01090"
            assert seq.global_mods[0].target_residues == ["C"]

            # Check other modifications
            if "?" in proforma:
                # Should have an ambiguous modification
                assert len(seq.mods.get(-4, [])) > 0
                assert any(mod.value == "Phospho" for mod in seq.mods.get(-4, []))
            elif "-" in proforma:
                # Should have an N-terminal modification
                assert len(seq.mods.get(-1, [])) > 0
                assert any(mod.value == "Acetyl" for mod in seq.mods.get(-1, []))

            # Check residue modification
            m_pos = 1  # Position of M in the sequence
            assert any(mod.value == "Oxidation" for mod in seq.seq[m_pos].mods)

            # Test roundtrip
            assert seq.to_proforma() == proforma

    def test_sequence_ambiguity(self):
        """Test representation of amino acid sequence ambiguity."""
        # Test simple ambiguity
        proforma = "(?DQ)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K"
        seq = Sequence.from_proforma(proforma)

        # Check sequence
        assert seq.to_stripped_string() == "NGTWEMESNENFEGYMK"

        # Check ambiguity
        assert len(seq.sequence_ambiguities) == 1
        assert seq.sequence_ambiguities[0].value == "DQ"
        assert seq.sequence_ambiguities[0].position == 0

        # Check modifications are still parsed correctly
        assert seq.seq[5].mods[0].value == "Oxidation"  # M position
        assert seq.seq[15].mods[0].value == "Oxidation"  # M position

        # Test roundtrip
        assert seq.to_proforma() == proforma

        # Test another example
        proforma2 = "(?N)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K"
        seq2 = Sequence.from_proforma(proforma2)

        assert seq2.sequence_ambiguities[0].value == "N"
        assert seq2.to_proforma() == proforma2

    def test_info_tags(self):
        """Test INFO tag support."""
        # Simple info tag
        proforma = "ELVIS[Phospho|INFO:newly discovered]K"
        seq = Sequence.from_proforma(proforma)

        # Check the modification has an info tag
        assert len(seq.seq[4].mods) == 1
        mod = seq.seq[4].mods[0]
        assert mod.value == "Phospho"
        print(mod.mod_value)
        assert len(mod.info_tags) == 1
        assert mod.info_tags[0] == "newly discovered"

        # Check roundtrip
        assert seq.to_proforma() == proforma

        # Multiple info tags
        proforma2 = "ELVIS[Phospho|INFO:newly discovered|INFO:Created on 2021-06]K"
        seq2 = Sequence.from_proforma(proforma2)

        mod2 = seq2.seq[4].mods[0]
        assert len(mod2.info_tags) == 2
        assert mod2.info_tags[0] == "newly discovered"
        assert mod2.info_tags[1] == "Created on 2021-06"

        # Check roundtrip
        assert seq2.to_proforma() == proforma2

    def test_info_tags_with_terminal_and_global_mods(self):
        """Test INFO tag support with terminal and global modifications."""
        # INFO tag with N-terminal modification
        proforma = "[Acetyl|INFO:Added during processing]-PEPTIDE"
        seq = Sequence.from_proforma(proforma)

        # Check N-terminal modification has info tag
        assert len(seq.mods[-1]) == 1  # N-terminal is position -1
        mod = seq.mods[-1][0]
        assert mod.value == "Acetyl"
        assert len(mod.info_tags) == 1
        assert mod.info_tags[0] == "Added during processing"

        # Check roundtrip
        assert seq.to_proforma() == proforma

        # INFO tag with C-terminal modification
        proforma = "PEPTIDE-[Amidated|INFO:Common C-terminal mod]"
        seq = Sequence.from_proforma(proforma)

        # Check C-terminal modification has info tag
        assert len(seq.mods[-2]) == 1  # C-terminal is position -2
        mod = seq.mods[-2][0]
        assert mod.value == "Amidated"
        assert len(mod.info_tags) == 1
        assert mod.info_tags[0] == "Common C-terminal mod"

        # Check roundtrip
        assert seq.to_proforma() == proforma

        # INFO tag with global isotope modification
        proforma = "<13C|INFO:Stable isotope labeling>PEPTIDE"
        seq = Sequence.from_proforma(proforma)

        # Check global modification has info tag
        assert len(seq.global_mods) == 1
        mod = seq.global_mods[0]
        assert mod.value == "13C"
        assert mod.global_mod_type == "isotope"
        assert len(mod.info_tags) == 1
        assert mod.info_tags[0] == "Stable isotope labeling"

        # Check roundtrip
        assert seq.to_proforma() == proforma

        # INFO tag with global fixed modification
        proforma = "<[Carbamidomethyl|INFO:Standard alkylation]@C>PEPTCDE"
        seq = Sequence.from_proforma(proforma)

        # Check global fixed modification has info tag
        assert len(seq.global_mods) == 1
        mod = seq.global_mods[0]
        assert mod.value == "Carbamidomethyl"
        assert mod.global_mod_type == "fixed"
        assert mod.target_residues == ["C"]
        assert len(mod.info_tags) == 1
        assert mod.info_tags[0] == "Standard alkylation"

        # Check roundtrip
        assert seq.to_proforma() == proforma

        # Complex case with multiple INFO tags and modifications
        proforma = "<[Carbamidomethyl|INFO:Standard alkylation]@C>[Acetyl|INFO:Added during processing]-PEPTCDE-[Amidated|INFO:Common C-terminal mod]"
        seq = Sequence.from_proforma(proforma)

        # Check all modifications have their respective info tags
        # Global mod
        assert len(seq.global_mods) == 1
        mod = seq.global_mods[0]
        assert mod.value == "Carbamidomethyl"
        assert len(mod.info_tags) == 1
        assert mod.info_tags[0] == "Standard alkylation"

        # N-terminal mod
        assert len(seq.mods[-1]) == 1
        mod = seq.mods[-1][0]
        assert mod.value == "Acetyl"
        assert len(mod.info_tags) == 1
        assert mod.info_tags[0] == "Added during processing"

        # C-terminal mod
        assert len(seq.mods[-2]) == 1
        mod = seq.mods[-2][0]
        assert mod.value == "Amidated"
        assert len(mod.info_tags) == 1
        assert mod.info_tags[0] == "Common C-terminal mod"

        # Check roundtrip
        assert seq.to_proforma() == proforma

    def test_joint_representation(self):
        """Test joint representation of experimental data and interpretation."""
        # Basic case with interpretation and mass
        proforma = "ELVIS[U:Phospho|+79.966331]K"
        seq = Sequence.from_proforma(proforma)
        mod = seq.seq[4].mods[0]
        assert mod.value == "Phospho"
        assert mod.source == "U"
        assert mod.mod_value.pipe_values[1].mass == 79.966331
        assert seq.to_proforma() == proforma

        proforma = "ELVIS[+79.966331]K"
        seq = Sequence.from_proforma(proforma)
        assert seq.to_proforma() == proforma

        # Case with observed mass
        proforma = "ELVIS[U:Phospho|Obs:+79.978]K"
        seq = Sequence.from_proforma(proforma)
        mod = seq.seq[4].mods[0]
        assert mod.value == "Phospho"
        assert mod.source == "U"
        assert mod.mod_value.pipe_values[1].observed_mass == 79.978
        assert seq.to_proforma() == proforma

        # Complex case with multiple synonyms
        proforma = "ELVIS[Phospho|O-phospho-L-serine|Obs:+79.966|INFO:Validated]K"
        seq = Sequence.from_proforma(proforma)
        mod = seq.seq[4].mods[0]
        assert mod.value == "Phospho"
        assert len(mod.synonyms) == 2
        assert mod.synonyms[1] == "O-phospho-L-serine"
        assert mod.mod_value.pipe_values[2].observed_mass == 79.966
        assert len(mod.info_tags) == 1
        assert mod.info_tags[0] == "Validated"

        # Check roundtrip
        assert seq.to_proforma() == proforma

    def test_crosslink_joint_representation(self):
        """Test joint representation of crosslinks with mass and other information."""
        # Crosslink with mass shift and info tag
        proforma = "PEPTK[XL:DSS#XL1|+138.068|INFO:reaction=NHS]IDE"
        seq = Sequence.from_proforma(proforma)

        mod = seq.seq[4].mods[0]
        assert mod.value == "DSS"
        assert mod.source == "XL"
        assert mod.crosslink_id == "XL1"
        assert abs(mod.mod_value.pipe_values[1].mass - 138.068) < 0.0001
        assert "reaction=NHS" in mod.info_tags

        # Check round-trip conversion
        assert seq.to_proforma() == proforma

        # Test with crosslink reference
        proforma = "PEPTK[XL:DSS#XL1|+138.068]IDEQR[#XL1]"
        seq = Sequence.from_proforma(proforma)

        # Check crosslink reference
        mod_ref = seq.seq[9].mods[0]
        assert mod_ref.is_crosslink_ref
        assert mod_ref.crosslink_id == "XL1"

        # Check round-trip
        assert seq.to_proforma() == proforma

    def test_branch_joint_representation(self):
        """Test joint representation of branches with other features."""
        # Branch with observed mass and synonym
        proforma = "PEPTK[DSS#BRANCH|Obs:+156.079|DSBU]IDE"
        seq = Sequence.from_proforma(proforma)

        mod = seq.seq[4].mods[0]
        assert mod.value == "DSS"
        assert mod.source == None
        assert mod.mod_value.pipe_values[0].is_branch
        assert mod.mod_value.pipe_values[1].observed_mass == 156.079
        assert "DSBU" in mod.synonyms

        # Check round-trip
        assert seq.to_proforma() == proforma

        # Branch with reference and info tag
        proforma = "PEPTK[XL:DSS#BRANCH|INFO:BranchData]IDE[#BRANCH]"
        seq = Sequence.from_proforma(proforma)

        # Check branch reference
        mod_ref = seq.seq[7].mods[0]
        assert mod_ref._is_branch_ref

        # Check round-trip
        assert seq.to_proforma() == proforma

    def test_ambiguity_group_joint_representation(self):
        """Test joint representation of ambiguity groups."""
        # Ambiguity group with localization score and mass
        proforma = "PEPT[U:Phospho#1(0.75)|+79.966|INFO:confidence=high]KIDE"
        seq = Sequence.from_proforma(proforma)

        mod = seq.seq[3].mods[0]
        assert mod.value == "Phospho"
        assert mod.source == "U"
        assert mod.ambiguity_group == "1"
        assert abs(mod.mod_value.pipe_values[0].localization_score - 0.75) < 0.0001
        assert abs(mod.mod_value.pipe_values[1].mass - 79.966) < 0.0001
        assert "confidence=high" in mod.info_tags

        # Check round-trip
        assert seq.to_proforma() == proforma

        # Ambiguity reference
        proforma = "PEPT[U:Phospho#1(0.75)]K[#1]IDE"
        seq = Sequence.from_proforma(proforma)

        # Check ambiguity reference
        mod_ref = seq.seq[4].mods[0]
        assert mod_ref.is_ambiguity_ref
        assert mod_ref.ambiguity_group == "1"

        # Check round-trip
        assert seq.to_proforma() == proforma

    def test_complex_multi_feature_representation(self):
        """Test complex joint representation with multiple features."""
        proforma = "PEP[U:Deamidation|+0.984]T[U:Phospho#1(0.75)|+79.966]K[XL:DSS#XL2|INFO:crosslinker]IDE[#BRANCH]R[#XL2]S[#1]"
        seq = Sequence.from_proforma(proforma)

        # Check deamidation
        mod1 = seq.seq[2].mods[0]
        assert mod1.value == "Deamidation"
        assert abs(mod1.mod_value.pipe_values[1].mass - 0.984) < 0.0001

        # Check phosphorylation with ambiguity group
        mod2 = seq.seq[3].mods[0]
        assert mod2.value == "Phospho"
        assert abs(mod2.mod_value.pipe_values[1].mass - 79.966) < 0.0001
        assert mod2.ambiguity_group == "1"
        assert abs(mod2.mod_value.pipe_values[0].localization_score - 0.75) < 0.0001

        # Check crosslink
        mod3 = seq.seq[4].mods[0]
        assert mod3.value == "DSS"
        assert mod3.crosslink_id == "XL2"
        assert "crosslinker" in mod3.info_tags

        # Check branch reference
        mod4 = seq.seq[9].mods[0]

        assert mod4.mod_value.pipe_values[0].is_ambiguity_ref

        # Check crosslink reference
        mod5 = seq.seq[8].mods[0]
        assert mod5.mod_value.pipe_values[0].is_crosslink_ref
        assert mod5.mod_value.pipe_values[0].crosslink_id == "XL2"

        # Check ambiguity reference
        mod6 = seq.seq[9].mods[0]
        assert mod6.mod_value.pipe_values[0].is_ambiguity_ref
        assert mod6.mod_value.pipe_values[0].ambiguity_group == "1"

        # Check round-trip
        assert seq.to_proforma() == proforma

    def test_duplicate_info_prevention(self):
        """Test prevention of duplicate information in joint representation."""
        # Test with duplicate modification
        proforma = "PEPT[U:Phospho|Phospho|+79.966|+79.966]IDE"
        seq = Sequence.from_proforma(proforma)

        # When converted back, duplicates should be removed
        assert seq.to_proforma() != "PEPT[U:Phospho|+79.966|Phospho]IDE"

        # Test with duplicate observed mass
        proforma = "PEPT[Phospho|+79.966|Obs:+79.968|Obs:+79.968]IDE"
        seq = Sequence.from_proforma(proforma)

        # When converted back, duplicates should be removed
        assert seq.to_proforma() == "PEPT[Phospho|+79.966|Obs:+79.968]IDE"


if __name__ == "__main__":
    unittest.main()
