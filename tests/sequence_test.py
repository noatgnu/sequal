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
        assert seq[2].mods[0].value == "Mass:+79.966"
        assert abs(seq[2].mods[0].mass - 79.966) < 0.0001

        # Check roundtrip (note: may not be identical due to internal representation)
        assert seq.to_proforma() == "PEP[Mass:+79.966]TIDE"

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
        assert seq[2].mods[0].value == "Mass:-17.027"
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

    def test_delta_mass_notation(self):
        """Test parsing and serialization of standard delta mass notations."""
        proforma_strings = [
            "EM[U:+15.9949]EVEES[U:+79.9663]PEK",
            "EM[U:+15.995]EVEES[U:+79.966]PEK",
            "EM[M:+15.9949]EVEES[R:+79.9663]PEK",
            "EM[X:+15.9949]EVEES[G:+79.9663]PEK",
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
        proforma_strings = [
            "SEQUEN[Formula:C12H20O2]CE",
            "SEQUEN[Formula:C12 H20 O2]CE",  # with spaces
            "SEQUEN[Formula:[13C2]CH6N]CE",  # with isotope
            "SEQUEN[Formula:HN-1O2]CE",  # with negative cardinality
            "SEQUEN[Formula:[13C2][12C-2]H2N]CE",  # complex isotope replacement
        ]

        for proforma in proforma_strings:
            seq = Sequence.from_proforma(proforma)

            # Verify the sequence itself
            assert seq.to_stripped_string() == "SEQUENCE"

            # Check modification
            mod = seq.seq[5].mods[0]
            assert mod.source == "Formula"

            # Verify roundtrip conversion
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


if __name__ == "__main__":
    unittest.main()
