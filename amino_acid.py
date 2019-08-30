from base_block import BaseBlock
from modification import Modification

AA_mass = {"A": 71.037114,
           "R":	156.101111,
           "N": 114.042927,
           "D": 115.026943,
           "C": 103.009185,
           "E": 129.042593,
           "Q": 128.058578,
           "G": 57.021464,
           "H": 137.058912,
           "I": 113.084064,
           "L": 113.084064,
           "K": 128.094963,
           "M": 131.040485,
           "F": 147.068414,
           "P": 97.052764,
           "S": 87.032028,
           "T": 101.047679,
           "U": 150.95363,
           "W": 186.079313,
           "Y": 163.06332,
           "V": 99.068414}


class AminoAcid(BaseBlock):
    def __init__(self, value, position=None, mass=None):
        super().__init__(value, position, branch=False, mass=mass)
        self.mods = []
        if not self.mass:
            if value in AA_mass:
                self.mass = AA_mass[value]

    def set_modification(self, i: Modification):
        self.mods.append(i)

    def __repr__(self):
        s = self.value
        for i in self.mods:
            s += str(i)
        return s

    def __str__(self):
        s = self.value
        for i in self.mods:
            s += str(i)
        return s


