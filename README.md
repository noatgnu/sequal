# SEQUAL / seq=

Sequal is a Python package for in-silico generation of modified sequences from a sequence input and modifications. It is designed to assist in protein engineering, mass spectrometry analysis, drug design, and other bioinformatics research.

## Features

- Generate all possible sequences with static and variable modifications.
- Support for custom modification annotations.
- Utilities for mass spectrometry fragment generation.

## Installation

To install Sequal, use pip:

```sh
pip install sequal
```

## Usage

### ProForma 2.0 Support

Sequal supports the [ProForma 2.0 standard](https://www.psidev.info/proforma) for proteoform notation, which provides a standardized way to represent protein sequences with their modifications.

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

# .mods property provides an access to all amino acids at this amino acid

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
