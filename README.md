# SEQUAL / seq=
----
![Test](https://github.com/noatgnu/sequal/actions/workflows/test.yml/badge.svg)

Sequal is a Python package for in-silico generation of modified sequences from a sequence input and modifications. It is designed to assist in protein engineering, mass spectrometry analysis, drug design, and other bioinformatics research.

## Features

- Full support for ProForma 2.1 standard for proteoform notation
- Generate all possible sequences with static and variable modifications
- Flexible modification handling with support for:
  - Formula validation for chemical modifications
  - Charged formulas (positive and negative charges)
  - Glycan structure validation
  - Custom monosaccharides in glycan compositions
  - Info tags and metadata
  - Observed mass recording
- Named entities (peptidoform, peptidoform ion, and compound ion names)
- Terminal global modifications (N-term, C-term, and terminal-specific)
- Placement controls for unknown position modifications (Position, Limit, CoMKP, CoMUP)
- Ion notation with semantic detection (a-type, b-type, c-type, x-type, y-type, z-type)
- Indexing and slicing for convenient access to modification values
- Support for custom modification annotations
- Sequence ambiguity representation
- Utilities for mass spectrometry fragment generation
- Labile and non-labile ion simulation

## Installation

To install Sequal, use pip:

```sh
pip install sequal
```

## Usage

### ProForma 2.1 Support

Sequal supports the [ProForma 2.1 standard](https://github.com/HUPO-PSI/ProForma) for proteoform notation, which provides a standardized way to represent protein sequences with their modifications.

#### Parsing ProForma notation

```python
from sequal.sequence import Sequence

# Basic ProForma notation with modification
seq = Sequence.from_proforma("ELVIS[Phospho]K")
print(seq.seq[4].value)  # S
print(seq.seq[4].mods[0].synonyms[0])  # Phospho

# ProForma with terminal modifications
seq = Sequence.from_proforma("[Acetyl]-PEPTIDE-[Amidated]")
print(seq.mods[-1][0].synonyms[0])  # Acetyl (N-terminal)
print(seq.mods[-2][0].synonyms[0])  # Amidated (C-terminal)

# ProForma with global modifications
seq = Sequence.from_proforma("<[Carbamidomethyl]@C>PEPTCDE")
print(seq.global_mods[0].synonyms[0])  # Carbamidomethyl
print(seq.global_mods[0].target_residues)  # ['C']

# ProForma with sequence ambiguity
seq = Sequence.from_proforma("PEPT(?DE|ID)E")
print(seq.sequence_ambiguities[0].sequence)  # DE|ID
print(seq.sequence_ambiguities[0].position)  # 4
```

#### Working with information tags

```python
from sequal.sequence import Sequence

# ProForma with info tags
seq = Sequence.from_proforma("ELVIS[Phospho|INFO:newly discovered]K")
mod = seq.seq[4].mods[0]
print(mod.synonyms[0])  # Phospho
print(mod.info_tags[0])  # newly discovered

# Multiple info tags
seq = Sequence.from_proforma("PEPTIDE-[Amidated|INFO:Common C-terminal mod|INFO:Added manually]")
mod = seq.mods[-2][0]  # C-terminal modification
print(mod.synonyms[0])  # Amidated
print(mod.info_tags)  # ['Common C-terminal mod', 'Added manually']
```

#### Joint representation of experimental data and interpretation

```python
from sequal.sequence import Sequence

# ProForma with joint interpretation and mass
seq = Sequence.from_proforma("ELVIS[U:Phospho|+79.966331]K")
mod = seq.seq[4].mods[0]
print(mod.mod_value.pipe_values[0].value)  # Phospho
print(mod.mod_value.pipe_values[0].source)  # U
print(mod.mod_value.pipe_values[1].mass)  # 79.966331

# ProForma with observed mass
seq = Sequence.from_proforma("ELVIS[Phospho|Obs:+79.978]K")
mod = seq.seq[4].mods[0]
print(mod.synonyms[0])  # Phospho
print(mod.mod_value.pipe_values[1].observed_mass)  # 79.978

# Complex case with synonyms, observed mass and info tags
seq = Sequence.from_proforma("ELVIS[Phospho|O-phospho-L-serine|Obs:+79.966|INFO:Validated]K")
mod = seq.seq[4].mods[0]
print(mod.synonyms[0])  # Phospho
print(mod.synonyms[1])  # O-phospho-L-serine
print(mod.mod_value.pipe_values[3].observed_mass)  # 79.966
print(mod.info_tags[0])  # Validated
```

#### Accessing mod_value with indexing and slicing

```python
from sequal.sequence import Sequence

# Using indexing to access pipe values
seq = Sequence.from_proforma("ELVIS[Unimod:21|Phospho|INFO:Validated]K")
mod = seq.seq[4].mods[0]
print(mod.mod_value[0].value)  # Primary pipe value - Unimod:21
print(mod.mod_value[1].value)  # Second pipe value - Phospho
print(mod.mod_value[2].type)   # PipeValue.INFO_TAG

# Using slicing to access multiple pipe values
pipe_values = mod.mod_value[1:3]  # Get second and third pipe values
for pv in pipe_values:
    print(f"{pv.type}: {pv.value}")

# Iterating through all pipe values
for i, pv in enumerate(mod.mod_value):
    print(f"Pipe value {i}: {pv.value}")

# Getting the length of pipe values
print(f"Number of pipe values: {len(mod.mod_value)}")
```

#### Working with formula and glycan modifications

```python
from sequal.sequence import Sequence

# Working with formula modifications
seq = Sequence.from_proforma("PEPTIDE[Formula:C2H3NO]")
mod = seq.seq[6].mods[0]
print(f"Formula: {mod.mod_value[0].value}")
print(f"Is valid formula: {mod.mod_value[0].is_valid}")

# Invalid formula - still parsed but marked invalid
seq = Sequence.from_proforma("PEPTIDE[Formula:123]")
mod = seq.seq[6].mods[0]
print(f"Invalid formula: {mod.mod_value[0].value}")
print(f"Is valid formula: {mod.mod_value[0].is_valid}")  # False

# Working with glycan modifications
seq = Sequence.from_proforma("PEPTID[Glycan:HexNAc1Hex3]E")
mod = seq.seq[5].mods[0]
print(f"Glycan: {mod.mod_value[0].value}")
print(f"Is valid glycan: {mod.mod_value[0].is_valid}")
print(f"Glycan pipe value type: {mod.mod_value[0].type}")  # PipeValue.GLYCAN

# Invalid glycan - still parsed but marked invalid and stored as SYNONYM type
seq = Sequence.from_proforma("PEPTID[Glycan:Invalid123]E")
mod = seq.seq[5].mods[0]
print(f"Invalid glycan: {mod.mod_value[0].value}")
print(f"Is valid glycan: {mod.mod_value[0].is_valid}")  # False
print(f"Invalid glycan pipe value type: {mod.mod_value[0].type}")  # PipeValue.SYNONYM
```

#### Named entities (ProForma 2.1)

```python
from sequal.sequence import Sequence

# Peptidoform name
seq = Sequence.from_proforma("(>Tryptic peptide)SEQUEN[Phospho]CE")
print(seq.peptidoform_name)  # Tryptic peptide

# Peptidoform ion name
seq = Sequence.from_proforma("(>>Precursor ion z=2)SEQUENCE")
print(seq.peptidoform_ion_name)  # Precursor ion z=2

# Compound ion name (for chimeric spectra)
seq = Sequence.from_proforma("(>>>MS2 Scan 1234)PEPTIDE/2+SEQUENCE/3")
print(seq.compound_ion_name)  # MS2 Scan 1234

# All three naming levels together
seq = Sequence.from_proforma("(>>>MS2)(>>Precursor)(>Albumin)PEPTIDE")
print(seq.compound_ion_name)      # MS2
print(seq.peptidoform_ion_name)   # Precursor
print(seq.peptidoform_name)       # Albumin
```

#### Charged formulas (ProForma 2.1)

```python
from sequal.sequence import Sequence

# Positive charge
seq = Sequence.from_proforma("SEQUEN[Formula:Zn1:z+2]CE")
mod = seq.seq[5].mods[0]
print(mod.mod_value[0].value)  # Zn1
print(mod.mod_value[0].charge)  # z+2
print(mod.mod_value[0].charge_value)  # 2

# Negative charge
seq = Sequence.from_proforma("PEPTIDE[Formula:C2H3NO:z-1]")
mod = seq.seq[6].mods[0]
print(mod.mod_value[0].charge)  # z-1
print(mod.mod_value[0].charge_value)  # -1
```

#### Custom monosaccharides (ProForma 2.1)

```python
from sequal.sequence import Sequence

# Custom monosaccharide in glycan composition
seq = Sequence.from_proforma("SEQUEN[Glycan:{C8H13N1O5}1Hex2]CE")
mod = seq.seq[5].mods[0]
print(mod.mod_value[0].value)  # {C8H13N1O5}1Hex2
print(mod.mod_value[0].is_valid)  # True

# Charged custom monosaccharide
seq = Sequence.from_proforma("N[Glycan:{C8H13N1O5Na1:z+1}1Hex2HexNAc2]")
mod = seq.seq[0].mods[0]
print(mod.mod_value[0].value)  # {C8H13N1O5Na1:z+1}1Hex2HexNAc2

# Multiple custom monosaccharides
seq = Sequence.from_proforma("N[Glycan:{C8H13N1O5}2{C6H10O5}1Hex3]")
mod = seq.seq[0].mods[0]
print(mod.mod_value[0].is_valid)  # True
```

#### Terminal global modifications (ProForma 2.1)

```python
from sequal.sequence import Sequence

# N-terminal global modification
seq = Sequence.from_proforma("<[Acetyl]@N-term>PEPTIDE")
print(seq.global_mods[0].synonyms[0])  # Acetyl
print(seq.global_mods[0].target_residues)  # ['N-term']

# C-terminal global modification
seq = Sequence.from_proforma("<[Amidated]@C-term>PEPTIDE")
print(seq.global_mods[0].target_residues)  # ['C-term']

# Terminal-specific: only N-terminal glutamine
seq = Sequence.from_proforma("<[Gln->pyro-Glu]@N-term:Q>QATPEILMCNSIGCLMG")
print(seq.global_mods[0].target_residues)  # {'N-term': ['Q']}

# Multiple terminal global modifications
seq = Sequence.from_proforma("<[TMT6plex]@K,N-term><[Oxidation]@M,C-term:G>MTPEILTCNSIGCLK")
print(len(seq.global_mods))  # 2
```

#### Placement controls (ProForma 2.1)

```python
from sequal.sequence import Sequence

# Position constraint: limit where modifications can occur
seq = Sequence.from_proforma("[Oxidation|Position:M]^4?PEPTIDEMETCM")
print(seq.to_proforma())

# Limit per position: allow multiple modifications at same position
seq = Sequence.from_proforma("[Oxidation|Limit:2]^4?PEPTIDE")
print(seq.to_proforma())

# CoMKP: allow colocalization with known position modifications
seq = Sequence.from_proforma("[Oxidation|CoMKP]?PEPT[Phospho]IDE")
print(seq.to_proforma())

# CoMUP: allow colocalization with unknown position modifications
seq = Sequence.from_proforma("PEPTIDE[Dioxidation|CoMUP][Oxidation|CoMUP]")
mod1 = seq.seq[6].mods[0]
mod2 = seq.seq[6].mods[1]
print(mod1.colocalize_unknown)  # True
print(mod2.colocalize_unknown)  # True

# Combined placement controls
seq = Sequence.from_proforma("[Oxidation|Position:M,C|Limit:2|CoMKP]^4?PEPTIDE")
print(seq.to_proforma())
```

#### Ion notation (ProForma 2.1)

```python
from sequal.sequence import Sequence

# b-type ion
seq = Sequence.from_proforma("PEPTIDE-[b-type-ion]")
c_term_mods = seq.mods[-2]
print(c_term_mods[0].is_ion_type)  # True
print(c_term_mods[0].value)  # b-type-ion

# a-type ion
seq = Sequence.from_proforma("PEPTIDE-[a-type-ion]")
c_term_mods = seq.mods[-2]
print(c_term_mods[0].is_ion_type)  # True

# Unimod ion type reference
seq = Sequence.from_proforma("PEPTIDE-[UNIMOD:2132]")  # b-type-ion
c_term_mods = seq.mods[-2]
print(c_term_mods[0].is_ion_type)  # True

# Complex ion notation with formula
seq = Sequence.from_proforma("PEPTID[Formula:H-1C-1O-2|Info:d-ion]-[a-type-ion]")
print(seq.to_proforma())

# Ion notation with charge
seq = Sequence.from_proforma("SFFLYSKLTV-[b-type-ion]/2")
print(seq.charge)  # 2
print(seq.mods[-2][0].is_ion_type)  # True
```

#### Converting to ProForma format

```python
from sequal.sequence import Sequence

# Parse and convert back to ProForma
proforma = "ELVIS[Phospho|INFO:newly discovered]K"
seq = Sequence.from_proforma(proforma)
print(seq.to_proforma())  # ELVIS[Phospho|INFO:newly discovered]K

# Complex example with multiple modification types
proforma = "<[Carbamidomethyl]@C>[Acetyl]-PEPTCDE-[Amidated]"
seq = Sequence.from_proforma(proforma)
print(seq.to_proforma())  # <[Carbamidomethyl]@C>[Acetyl]-PEPTCDE-[Amidated]

# ProForma 2.1 features combined
proforma = "(>Tryptic peptide)<[TMT6plex]@K,N-term>SEQUEN[Formula:Zn1:z+2]CE-[b-type-ion]/2"
seq = Sequence.from_proforma(proforma)
print(seq.to_proforma())  # Perfect roundtrip preservation
```

### Sequence comprehension

#### Using Sequence Object with Unmodified Protein Sequence

```python
from sequal.sequence import Sequence
#Using Sequence object with unmodified protein sequence

seq = Sequence("TESTEST")
print(seq.seq) #should print "TESTEST"
print(seq[0:2]) #should print "TE"
```

#### Using Sequence Object with Modified Protein Sequence

```python
from sequal.sequence import Sequence
#Using Sequence object with modified protein sequence. []{}() could all be used as modification annotation.

seq = Sequence("TEN[HexNAc]ST")
for i in seq.seq:
    print(i, i.mods) #should print N [HexNAc] on the 3rd amino acid

seq = Sequence("TEN[HexNAc][HexNAc]ST")
for i in seq.seq:
    print(i, i.mods) #should print N [HexNAc, HexNAc] on the 3rd amino acid

# .mods property provides access to an arrays of all modifications at this amino acid

seq = Sequence("TE[HexNAc]NST", mod_position="left") #mod_position left indicate that the modification should be on the left of the amino acid instead of default which is right
for i in seq.seq:
    print(i, i.mods) #should print N [HexNAc] on the 3rd amino acid
```

#### Custom Annotation Formatting

```python
from sequal.sequence import Sequence
#Format sequence with custom annotation
seq = Sequence("TENST")
a = {1:"tes", 2:["1", "200"]}
print(seq.to_string_customize(a, individual_annotation_enclose=False, individual_annotation_separator="."))
# By supplying .to_string_customize with a dictionary of position on the sequence that you wish to annotate
# The above would print out TE[tes]N[1.200]ST
```

### Modification

#### Creating a Modification Object

```python
from sequal.modification import Modification

# Create a modification object and try to find all its possible positions using regex
mod = Modification("HexNAc", regex_pattern="N[^P][S|T]")
for ps, pe in mod.find_positions("TESNEST"):
    print(ps, pe)
    # this should print out the position 3 on the sequence as the start of the match and position 6 as the end of the match
```

### Generating Modified Sequences

Static Modification

```python
from sequal.sequence import ModdedSequenceGenerator
from sequal.modification import Modification

propiona = Modification("Propionamide", regex_pattern="C", mod_type="static")
seq = "TECSNTT"
mods = [propiona]
g = ModdedSequenceGenerator(seq, static_mods=mods)
for i in g.generate():
    print(i)  # should print {2: [Propionamide]}
```

Variable Modification

```python
from sequal.sequence import ModdedSequenceGenerator
from sequal.modification import Modification

nsequon = Modification("HexNAc", regex_pattern="N[^P][S|T]", mod_type="variable", labile=True)
osequon = Modification("Mannose", regex_pattern="[S|T]", mod_type="variable", labile=True)
carbox = Modification("Carboxylation", regex_pattern="E", mod_type="variable", labile=True)

seq = "TECSNTT"
mods = [nsequon, osequon, carbox]
g = ModdedSequenceGenerator(seq, mods, [])
print(g.variable_map.mod_position_dict)
# should print {'HexNAc0': [3], 'Mannose0': [0, 2, 4, 5, 6], 'Carboxylation0': [1]}

for i in g.generate():
    print(i)
    # should print all possible combinations of variable modifications
```

### Mass spectrometry utilities

Generating Non-Labile and Labile Ions

```python
from sequal.mass_spectrometry import fragment_non_labile, fragment_labile
from sequal.modification import Modification
from sequal.sequence import ModdedSequenceGenerator, Sequence

nsequon = Modification("HexNAc", regex_pattern="N[^P][S|T]", mod_type="variable", labile=True, labile_number=1, mass=203)
propiona = Modification("Propionamide", regex_pattern="C", mod_type="static", mass=71)

seq = "TECSNTT"
static_mods = [propiona]
variable_mods = [nsequon]

g = ModdedSequenceGenerator(seq, variable_mods, static_mods)
for i in g.generate():
    print(i)
    s = Sequence(seq, mods=i)
    for b, y in fragment_non_labile(s, "by"):
        print(b, "b{}".format(b.fragment_number))
        print(y, "y{}".format(y.fragment_number))

g = ModdedSequenceGenerator(seq, variable_mods, static_mods)
for i in g.generate():
    s = Sequence(seq, mods=i)
    ion = fragment_labile(s)
    if ion.has_labile:
        print(ion, "Y{}".format(ion.fragment_number))
        print(ion.mz_calculate(1))
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
