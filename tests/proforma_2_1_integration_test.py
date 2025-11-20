"""
Comprehensive integration tests for ProForma 2.1 features.

Tests complex scenarios combining multiple ProForma 2.1 features:
- Named Entities (Section 8.2)
- Terminal Global Modifications (Section 11.3.2)
- Custom Monosaccharides (Section 10.2)
- Charged Formulas (Section 11.1)
- Ion Notation (Section 11.6)
- Placement Controls (Section 11.2)
"""

import unittest

from sequal.sequence import Sequence


class TestProForma21Integration(unittest.TestCase):
    """Comprehensive integration tests for ProForma 2.1 features."""

    # ========== Named Entities + Multiple Features ==========

    def test_peptidoform_name_with_terminal_global_mods_and_charge(self):
        """Test peptidoform name with terminal global mods and charge state."""
        seq = Sequence.from_proforma(
            "(>TMT-labeled peptide)<[TMT6plex]@K,N-term>PEPTIDEK/2"
        )

        self.assertEqual(seq.peptidoform_name, "TMT-labeled peptide")
        self.assertEqual(len(seq.global_mods), 1)
        self.assertIn("K", seq.global_mods[0].target_residues)
        self.assertIn(
            {"type": "terminal", "terminal": "N-term"},
            seq.global_mods[0].target_residues,
        )
        self.assertEqual(seq.charge, 2)

    def test_all_three_naming_levels_with_modifications(self):
        """Test all three naming levels together with modifications."""
        seq = Sequence.from_proforma(
            "(>>>Chimeric Spectrum 1234)(>>Precursor z=2)(>Phospho-peptide)<[Phospho]@S,T,Y>PEPT[Oxidation]IDES/2"
        )

        self.assertEqual(seq.compound_ion_name, "Chimeric Spectrum 1234")
        self.assertEqual(seq.peptidoform_ion_name, "Precursor z=2")
        self.assertEqual(seq.peptidoform_name, "Phospho-peptide")
        self.assertEqual(len(seq.global_mods), 1)
        self.assertEqual(seq.seq[3].mods[0].mod_value.primary_value, "Oxidation")
        self.assertEqual(seq.charge, 2)

    def test_named_entity_with_labile_glycans(self):
        """Test named entity with labile glycans."""
        seq = Sequence.from_proforma(
            "(>Glycopeptide){Glycan:{C8H13N1O5}1Hex2}PEPN[Glycan:HexNAc2Hex3]TIDE"
        )

        self.assertEqual(seq.peptidoform_name, "Glycopeptide")
        labile_mods = seq.mods.get(-3)
        self.assertIsNotNone(labile_mods)
        self.assertTrue(labile_mods[0].mod_value.pipe_values[0].is_valid_glycan)
        self.assertEqual(seq.seq[3].mods[0].mod_value.source, "Glycan")

    # ========== Terminal Global Modifications + Complex Scenarios ==========

    def test_multiple_terminal_global_mods_with_placement_controls(self):
        """Test multiple terminal global mods with placement controls."""
        seq = Sequence.from_proforma(
            "<[TMT6plex|Limit:1]@K,N-term><[Oxidation|Position:M,C]@M,C-term:G>MTPEILTCNSIGCLKG"
        )

        self.assertEqual(len(seq.global_mods), 2)
        self.assertIn("K", seq.global_mods[0].target_residues)
        self.assertIn(
            {"type": "terminal", "terminal": "N-term"},
            seq.global_mods[0].target_residues,
        )
        self.assertIn("M", seq.global_mods[1].target_residues)
        self.assertIn(
            {"type": "terminal_specific", "terminal": "C-term", "amino_acid": "G"},
            seq.global_mods[1].target_residues,
        )

    def test_terminal_global_mods_with_ambiguous_modifications(self):
        """Test terminal global mods with ambiguous modifications."""
        seq = Sequence.from_proforma(
            "<[Gln->pyro-Glu]@N-term:Q>QPEPTIDE[Phospho#1]S[#1]T"
        )

        self.assertEqual(len(seq.global_mods), 1)
        self.assertIn(
            {"type": "terminal_specific", "terminal": "N-term", "amino_acid": "Q"},
            seq.global_mods[0].target_residues,
        )
        self.assertEqual(seq.seq[7].mods[0].ambiguity_group, "1")
        self.assertTrue(seq.seq[8].mods[0].is_ambiguity_ref)

    def test_terminal_global_mods_with_charge_states(self):
        """Test terminal global mods with charge states."""
        seq = Sequence.from_proforma("<[Acetyl]@N-term><[TMT6plex]@K>PEPTIDEK/2")

        self.assertEqual(len(seq.global_mods), 2)
        self.assertIn(
            {"type": "terminal", "terminal": "N-term"},
            seq.global_mods[0].target_residues,
        )
        self.assertEqual(seq.global_mods[1].target_residues, ["K"])
        self.assertEqual(seq.charge, 2)

    # ========== Custom Monosaccharides + Complex Glycosylation ==========

    def test_mixed_custom_and_standard_monosaccharides(self):
        """Test mixed custom and standard monosaccharides."""
        seq = Sequence.from_proforma("PEPN[Glycan:HexNAc2{C8H13N1O5Na1:z+1}1Hex3]TIDE")

        self.assertEqual(seq.seq[3].mods[0].mod_value.source, "Glycan")
        self.assertTrue(seq.seq[3].mods[0].mod_value.pipe_values[0].is_valid_glycan)

    def test_labile_and_static_custom_monosaccharides_together(self):
        """Test labile and static custom monosaccharides together."""
        seq = Sequence.from_proforma(
            "{Glycan:{C11H17N1O9}1Hex2}PEPN[Glycan:{C8H13[15N1]O5}2Hex1]TIDE"
        )

        labile_mods = seq.mods.get(-3)
        self.assertIsNotNone(labile_mods)
        self.assertTrue(labile_mods[0].mod_value.pipe_values[0].is_valid_glycan)
        self.assertTrue(seq.seq[3].mods[0].mod_value.pipe_values[0].is_valid_glycan)
        self.assertIn("[15N1]", seq.seq[3].mods[0].mod_value.pipe_values[0].value)

    def test_custom_monosaccharides_in_sequences(self):
        """Test custom monosaccharides in sequences."""
        seq = Sequence.from_proforma("NSTPEPN[Glycan:{C8H13N1O5}1Hex2]TIDE")

        self.assertEqual(seq.seq[6].mods[0].mod_value.source, "Glycan")
        self.assertTrue(seq.seq[6].mods[0].mod_value.pipe_values[0].is_valid_glycan)

    # ========== Charged Formulas + Ion Notation ==========

    def test_charged_formula_with_ion_type_notation(self):
        """Test charged formula with ion type notation."""
        seq = Sequence.from_proforma("PEPTIDE[Formula:Zn1:z+2]-[b-type-ion]")

        self.assertEqual(seq.seq[6].mods[0].mod_value.pipe_values[0].charge_value, 2)
        c_term_mods = seq.mods.get(-2)
        self.assertIsNotNone(c_term_mods)
        self.assertTrue(c_term_mods[0].is_ion_type)

    def test_multiple_charged_formulas_with_different_charges(self):
        """Test multiple charged formulas with different charges."""
        seq = Sequence.from_proforma("PEPT[Formula:C2H3NO:z-1]IDE[Formula:Zn1:z+2]K")

        self.assertEqual(seq.seq[3].mods[0].mod_value.pipe_values[0].charge_value, -1)
        self.assertEqual(seq.seq[6].mods[0].mod_value.pipe_values[0].charge_value, 2)

    def test_charged_glycan_with_ion_notation(self):
        """Test charged glycan with ion notation."""
        seq = Sequence.from_proforma(
            "N[Glycan:{C8H13N1O5Na1:z+1}1Hex2]PEPTIDE-[y-type-ion]"
        )

        self.assertTrue(seq.seq[0].mods[0].mod_value.pipe_values[0].is_valid_glycan)
        c_term_mods = seq.mods.get(-2)
        self.assertTrue(c_term_mods[0].is_ion_type)

    # ========== Placement Controls + Advanced Features ==========

    def test_basic_placement_controls_on_fixed_position_mods(self):
        """Test basic placement controls on fixed position mods."""
        seq = Sequence.from_proforma("PEPT[Phospho|INFO:High confidence]IDE")

        self.assertEqual(seq.seq[3].mods[0].mod_value.primary_value, "Phospho")
        self.assertGreater(len(seq.seq[3].mods[0].info_tags), 0)

    def test_multiple_modifications_with_ambiguity(self):
        """Test multiple modifications with ambiguity."""
        seq = Sequence.from_proforma("PEPT[Phospho#G1]IDE[#G1]")

        self.assertEqual(seq.seq[3].mods[0].ambiguity_group, "G1")
        self.assertTrue(seq.seq[6].mods[0].is_ambiguity_ref)

    # ========== Chimeric Spectra with ProForma 2.1 ==========

    def test_chimeric_with_named_entities(self):
        """Test chimeric with named entities."""
        seq = Sequence.from_proforma(
            "(>>>Chimeric Scan 5678)(>Peptide1)PEPTIDE+SEQUENCEK"
        )

        self.assertEqual(seq.compound_ion_name, "Chimeric Scan 5678")
        self.assertEqual(seq.peptidoform_name, "Peptide1")
        self.assertTrue(seq.is_chimeric)

    # ========== Multi-Chain Sequences with ProForma 2.1 ==========

    def test_multi_chain_with_naming_level_on_first_chain(self):
        """Test multi-chain with naming level on first chain."""
        seq = Sequence.from_proforma(
            "(>>>Disulfide-linked)(>>Ion pair)(>Chain1)PEPTIDEC//CSEQUENCE"
        )

        self.assertEqual(seq.compound_ion_name, "Disulfide-linked")
        self.assertEqual(seq.peptidoform_ion_name, "Ion pair")
        self.assertTrue(seq.is_multi_chain)
        self.assertEqual(seq.chains[0].peptidoform_name, "Chain1")

    def test_multi_chain_with_different_terminal_global_mods(self):
        """Test multi-chain with different terminal global mods."""
        seq = Sequence.from_proforma(
            "<[Acetyl]@N-term>PEPTIDE//<[TMT6plex]@K,N-term>KSEQUENCEK"
        )

        self.assertTrue(seq.is_multi_chain)
        self.assertIn(
            {"type": "terminal", "terminal": "N-term"},
            seq.chains[0].global_mods[0].target_residues,
        )
        self.assertIn("K", seq.chains[1].global_mods[0].target_residues)
        self.assertIn(
            {"type": "terminal", "terminal": "N-term"},
            seq.chains[1].global_mods[0].target_residues,
        )

    # ========== Real-World Complex Scenarios ==========

    def test_tmt_labeled_phosphopeptide_with_multiple_ptms(self):
        """Test TMT-labeled phosphopeptide with multiple PTMs."""
        seq = Sequence.from_proforma(
            "(>TMT-labeled phosphopeptide)<[TMT6plex]@K,N-term><[Oxidation]@M>MTPEILTS[Phospho]CNSIGCLK/2"
        )

        self.assertEqual(seq.peptidoform_name, "TMT-labeled phosphopeptide")
        self.assertEqual(len(seq.global_mods), 2)
        self.assertIn(
            {"type": "terminal", "terminal": "N-term"},
            seq.global_mods[0].target_residues,
        )
        self.assertEqual(seq.global_mods[1].target_residues, ["M"])
        self.assertEqual(seq.seq[7].mods[0].mod_value.primary_value, "Phospho")
        self.assertEqual(seq.charge, 2)

    def test_glycopeptide_with_multiple_glycosylation_sites(self):
        """Test glycopeptide with multiple glycosylation sites."""
        seq = Sequence.from_proforma(
            "{Glycan:Hex1HexNAc1}N[Glycan:{C8H13N1O5}2Hex3HexNAc2]GSTN[Glycan:HexNAc2Hex3NeuAc1]VTPEPTIDE"
        )

        labile_mods = seq.mods.get(-3)
        self.assertIsNotNone(labile_mods)
        self.assertTrue(seq.seq[0].mods[0].mod_value.pipe_values[0].is_valid_glycan)
        self.assertTrue(seq.seq[4].mods[0].mod_value.pipe_values[0].is_valid_glycan)

    def test_complex_proteoform_with_multiple_modification_types(self):
        """Test complex proteoform with multiple modification types."""
        seq = Sequence.from_proforma(
            "(>Ubiquitinated protein)[Acetyl]-MRSGSHHHHHHGSPEPTM[Oxidation]IDEK[Ubiquitin#BRANCH]SEQUENCE-[Amidated]"
        )

        self.assertEqual(seq.peptidoform_name, "Ubiquitinated protein")
        n_term_mods = seq.mods.get(-1)
        self.assertEqual(n_term_mods[0].mod_value.primary_value, "Acetyl")
        self.assertEqual(seq.seq[17].mods[0].mod_value.primary_value, "Oxidation")
        self.assertTrue(seq.seq[21].mods[0].mod_value.is_branch)
        c_term_mods = seq.mods.get(-2)
        self.assertEqual(c_term_mods[0].mod_value.primary_value, "Amidated")

    def test_metal_binding_peptide_with_charged_formulas(self):
        """Test metal-binding peptide with charged formulas."""
        seq = Sequence.from_proforma(
            "(>Zinc finger peptide)CAQECGK[Formula:Zn1:z+2]SFTSALK[Formula:Zn1:z+2]SRHK/3"
        )

        self.assertEqual(seq.peptidoform_name, "Zinc finger peptide")
        self.assertEqual(seq.seq[6].mods[0].mod_value.pipe_values[0].charge_value, 2)
        self.assertEqual(seq.seq[13].mods[0].mod_value.pipe_values[0].charge_value, 2)
        self.assertEqual(seq.charge, 3)

    def test_peptide_with_unknown_position_mods(self):
        """Test peptide with unknown position mods."""
        seq = Sequence.from_proforma("[Phospho]^3[Oxidation]^2?MSTPEPTMSTY")

        unknown_mods = seq.mods.get(-4)
        self.assertEqual(len(unknown_mods), 5)
        self.assertEqual(unknown_mods[0].mod_value.primary_value, "Phospho")
        self.assertEqual(unknown_mods[3].mod_value.primary_value, "Oxidation")

    # ========== Edge Cases and Complex Combinations ==========

    def test_all_proforma_21_features_in_one_sequence(self):
        """Test all ProForma 2.1 features in one sequence."""
        seq = Sequence.from_proforma(
            "(>>>MS2 Spectrum)(>>Precursor z=2)(>Complex peptide)<[TMT6plex]@K,N-term><[Gln->pyro-Glu]@N-term:Q>"
            + "[Acetyl]-QN[Glycan:HexNAc2Hex3]PEPTM[Oxidation|INFO:High confidence]IDEK[Formula:Zn1:z+2]-[b-type-ion]/2"
        )

        self.assertEqual(seq.compound_ion_name, "MS2 Spectrum")
        self.assertEqual(seq.peptidoform_ion_name, "Precursor z=2")
        self.assertEqual(seq.peptidoform_name, "Complex peptide")
        self.assertEqual(len(seq.global_mods), 2)

        n_term_mods = seq.mods.get(-1)
        self.assertEqual(n_term_mods[0].mod_value.primary_value, "Acetyl")

        self.assertEqual(seq.seq[1].mods[0].mod_value.source, "Glycan")
        self.assertGreater(len(seq.seq[6].mods[0].info_tags), 0)
        self.assertEqual(seq.seq[10].mods[0].mod_value.pipe_values[0].charge_value, 2)

        c_term_mods = seq.mods.get(-2)
        self.assertTrue(c_term_mods[0].is_ion_type)

        self.assertEqual(seq.charge, 2)

    def test_round_trip_complex_integrated_sequences(self):
        """Test round-trip complex integrated sequences."""
        test_cases = [
            "(>Named peptide)<[TMT6plex]@K,N-term>PEPTIDEK/2",
            "N[Glycan:HexNAc2{C8H13N1O5}1]PEPTIDE",
            "<[Gln->pyro-Glu]@N-term:Q>QPEPTM[Oxidation]IDE-[b-type-ion]",
            "[Phospho]^2?STPEPTIDE",
        ]

        for original in test_cases:
            seq = Sequence.from_proforma(original)
            proforma = seq.to_proforma()
            seq2 = Sequence.from_proforma(proforma)
            proforma2 = seq2.to_proforma()

            self.assertEqual(proforma, original)
            self.assertEqual(proforma2, original)

    def test_terminal_modifications_with_global_terminal_mods(self):
        """Test terminal modifications with global terminal mods."""
        seq = Sequence.from_proforma("<[Acetyl]@N-term>[Carbamyl]-PEPTIDE-[Methyl]")

        self.assertEqual(len(seq.global_mods), 1)
        self.assertIn(
            {"type": "terminal", "terminal": "N-term"},
            seq.global_mods[0].target_residues,
        )

        n_term_mods = seq.mods.get(-1)
        self.assertEqual(n_term_mods[0].mod_value.primary_value, "Carbamyl")

        c_term_mods = seq.mods.get(-2)
        self.assertEqual(c_term_mods[0].mod_value.primary_value, "Methyl")

    def test_sequence_ambiguity_with_proforma_21_features(self):
        """Test sequence ambiguity with ProForma 2.1 features."""
        seq = Sequence.from_proforma(
            "(>Ambiguous peptide)<[TMT6plex]@K>PEPT(?IL)DE[Formula:Zn1:z+2]K/2"
        )

        self.assertEqual(seq.peptidoform_name, "Ambiguous peptide")
        self.assertEqual(len(seq.global_mods), 1)
        self.assertEqual(len(seq.sequence_ambiguities), 1)
        self.assertEqual(seq.sequence_ambiguities[0].value, "IL")
        self.assertEqual(seq.seq[5].mods[0].mod_value.pipe_values[0].charge_value, 2)
        self.assertEqual(seq.charge, 2)

    def test_isotope_labels_with_proforma_21_features(self):
        """Test isotope labels with ProForma 2.1 features."""
        seq = Sequence.from_proforma(
            "(>Heavy labeled)<[Label:13C(6)15N(2)]@K,R>PEPTIDER"
        )

        self.assertEqual(seq.peptidoform_name, "Heavy labeled")
        self.assertEqual(len(seq.global_mods), 1)
        self.assertIn("13C", seq.global_mods[0].mod_value.primary_value)
        self.assertIn("15N", seq.global_mods[0].mod_value.primary_value)

    # ========== Serialization Integrity ==========

    def test_maintain_order_of_multiple_global_modifications(self):
        """Test maintain order of multiple global modifications."""
        original = "<[TMT6plex]@K,N-term><[Oxidation]@M><[Phospho]@S,T,Y>MSTPEPTIDEK"
        seq = Sequence.from_proforma(original)
        proforma = seq.to_proforma()

        self.assertEqual(proforma, original)
        self.assertEqual(len(seq.global_mods), 3)

    def test_maintain_complex_nested_structures(self):
        """Test maintain complex nested structures."""
        original = "(>>>C1)(>>I1)(>P1)<[Mod]@N-term>{Glycan:Hex1}[Ac]-N[Gly:H1]S-[Am]/2"
        seq = Sequence.from_proforma(original)
        proforma = seq.to_proforma()

        self.assertEqual(seq.compound_ion_name, "C1")
        self.assertEqual(seq.peptidoform_ion_name, "I1")
        self.assertEqual(seq.peptidoform_name, "P1")
        self.assertEqual(proforma, original)

    # ========== Performance and Stress Tests ==========

    def test_very_long_sequence_with_multiple_modifications(self):
        """Test very long sequence with multiple modifications."""
        long_seq = "M" * 50
        modified_seq = "".join(
            f"{aa}[Oxidation]" if i % 5 == 0 else aa for i, aa in enumerate(long_seq)
        )
        seq = Sequence.from_proforma(f"<[Oxidation]@M>{modified_seq}")

        self.assertEqual(len(seq.global_mods), 1)
        self.assertEqual(len(seq.seq), 50)
        mods_count = sum(1 for aa in seq.seq if len(aa.mods) > 0)
        self.assertEqual(mods_count, 10)

    def test_many_unknown_position_modifications(self):
        """Test many unknown position modifications."""
        seq = Sequence.from_proforma("[Phospho|Position:S,T,Y]^10?SYTSTPEPTIDESTSTY")

        unknown_mods = seq.mods.get(-4)
        self.assertEqual(len(unknown_mods), 10)

    # ========== Glycan Count Validation ==========

    def test_validate_glycan_count_requirements(self):
        """Test glycan count validation rules."""
        valid_cases = [
            ("HexNAc", "at end, count 1 implied"),
            ("HexNAc2Hex", "HexNAc has count 2, Hex at end"),
            ("Hex(3)HexNAc2", "explicit counts"),
            ("{C8H13N1O5}1Hex2", "custom with count, Hex with count"),
        ]

        for glycan, reason in valid_cases:
            with self.subTest(glycan=glycan, reason=reason):
                seq = Sequence.from_proforma(f"N[Glycan:{glycan}]K")
                self.assertTrue(
                    seq.seq[0].mods[0].mod_value.pipe_values[0].is_valid_glycan,
                    f"Expected {glycan} to be valid ({reason})",
                )

    def test_reject_invalid_glycan_counts(self):
        """Test rejection of invalid glycan counts."""
        invalid_cases = [
            ("HexNAcHex", "HexNAc not at end, missing count"),
            ("HexNAc0", "zero count not allowed"),
            ("{C8H13N1O5}Hex2", "custom not at end, missing count"),
        ]

        for glycan, reason in invalid_cases:
            with self.subTest(glycan=glycan, reason=reason):
                seq = Sequence.from_proforma(f"N[Glycan:{glycan}]K")
                self.assertFalse(
                    seq.seq[0].mods[0].mod_value.pipe_values[0].is_valid_glycan,
                    f"Expected {glycan} to be invalid ({reason})",
                )


if __name__ == "__main__":
    unittest.main()
