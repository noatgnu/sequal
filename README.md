# SEQUAL / seq=

Sequal is developed as a python package for in-silico generation of modified sequences from a sequence input and modifications.

## Dependencies

`None.`

## Usage
Sequence comprehension
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
    print(i, i.mods) #should print N [HexNAc, HexNAc] on the 3rd amino acid
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