# SEQUAL / seq=

Sequal is developed as a python package for in-silico generation of modified sequences from a sequence input and modifications.

## Dependencies

`None.`

## Usage


### Sequence comprehension


```python
from sequal.sequence import Sequence
#Using Sequence object with unmodified protein sequence

seq = Sequence("TESTEST")
print(seq.seq) #should print "TESTEST"
print(seq[0:2]) #should print "TE"
```

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

```python
from sequal.modification import Modification

# Create a modification object and try to find all its possible positions using regex
mod = Modification("HexNAc", regex_pattern="N[^P][S|T]")
for ps, pe in mod.find_positions("TESNEST"):
    print(ps, pe)
    # this should print out the position 3 on the sequence as the start of the match and position 6 as the end of the match
```

```python
from sequal.sequence import ModdedSequenceGenerator
from sequal.modification import Modification

# Examples for generation of modification combinations for a specific peptide

nsequon = Modification("HexNAc",regex_pattern="N[^P][S|T]", mod_type="variable", labile=True)
osequon = Modification("Mannose",regex_pattern="[S|T]", mod_type="variable", labile=True)
sulfation = Modification("Sulfation",regex_pattern="S", mod_type="variable", labile=True)
carbox = Modification("Carboxylation",regex_pattern="E", mod_type="variable", labile=True)
carbox2 = Modification("Carboxylation2", regex_pattern="E", mod_type="variable", labile=True, mass=43.98983)
propiona = Modification("Propionamide", regex_pattern="C", mod_type="static")

#Static modification 

seq = "TECSNTT"
mods = [propiona]
g = ModdedSequenceGenerator(seq, static_mods=mods)
for i in g.generate():
    print(i)
    # this would print out a dictionary with key being the position of modifications and values being an array of all modifications can be generated at that site
    # {2: [Propionamide]}

#Variable modification

mods = [nsequon, osequon, carbox]
g = ModdedSequenceGenerator(seq, mods, [])
print(g.variable_map.mod_position_dict)
#The object when supplied with an array of variable modifications would also create a dictionary of the modification and where they may be found
#{'HexNAc0': [3], 'Mannose0': [0, 2, 4, 5, 6], 'Carboxylation0': [1]}

for i in g.generate():
    print(i)
    #Similar to the static modification example, this will create all possible combinations of variable modifications however it would also include ones without the modification
    # {}
    # {1: [Carboxylation0]}
    # {6: [Mannose0]}
    # {6: [Mannose0], 1: [Carboxylation0]}
    # {5: [Mannose0]}
    # {5: [Mannose0], 1: [Carboxylation0]}
    # {5: [Mannose0], 6: [Mannose0]}
    # {5: [Mannose0], 6: [Mannose0], 1: [Carboxylation0]}
    # ...

```
### Mass spectrometry utilities

Here is an examples for usage of the `mass_spectrometry` module within `sequal` in combination with modified sequence generation

```python
from sequal.mass_spectrometry import fragment_non_labile, fragment_labile
from sequal.modification import Modification
from sequal.sequence import ModdedSequenceGenerator, Sequence

nsequon = Modification("HexNAc",regex_pattern="N[^P][S|T]", mod_type="variable", labile=True, labile_number=1, mass=203)
propiona = Modification("Propionamide", regex_pattern="C", mod_type="static", mass=71)


seq = "TECSNTT"
static_mods = [propiona]
variable_mods = [nsequon]

# Generating non labile b- and y- ions

g = ModdedSequenceGenerator(seq, variable_mods, static_mods)
for i in g.generate():
    print(i)
    # Print the combination of modifications
    s = Sequence(seq, mods=i)
    print(s)
    # Create new modified sequence object using the generated modifications and raw sequence string
    for b, y in fragment_non_labile(s, "by"):
        print(b, "b{}".format(b.fragment_number))
        print(y, "y{}".format(y.fragment_number))
        # Generate b- and y- non_labile ions objects based on the newly created sequence object

# Generating only labile ions based on the input modification

g = ModdedSequenceGenerator(seq, variable_mods, static_mods)
for i in g.generate():
    s = Sequence(seq, mods=i)
    ion = fragment_labile(s)
    if ion.has_labile:
        print(ion, "Y{}".format(ion.fragment_number))
        # TEC[Propionamide]SN[HexNAc]TT Y1
        print(ion.mz_calculate(1))
        # 1011.277047

```

