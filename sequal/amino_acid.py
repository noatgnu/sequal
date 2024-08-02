from sequal.base_block import BaseBlock
from sequal.modification import Modification
from sequal.resources import AA_mass


class AminoAcid(BaseBlock):
    """
    Represents an amino acid block that can carry position, modifications, and amino acid value.

    Inherits from the BaseBlock class and adds functionality specific to amino acids, such as
    handling modifications and inferring mass from a predefined dictionary.

    :param value: str
        The name of the amino acid residue for this block.
    :param position: int, optional
        The position of the amino acid residue that this block belongs to (default is None).
    :param mass: float, optional
        The mass of the amino acid. If not specified, it will try to infer the mass from the internal
        hard-coded mass dictionary of amino acids (default is None).
    """

    def __init__(self, value, position=None, mass=None):
        """
        Initialize an AminoAcid object.

        :param value: str
            The name of the amino acid residue for this block.
        :param position: int, optional
            The position of the amino acid residue that this block belongs to (default is None).
        :param mass: float, optional
            The mass of the amino acid. If not specified, it will try to infer the mass from the internal
            hard-coded mass dictionary of amino acids (default is None).
        """
        super().__init__(value, position, branch=False, mass=mass)
        self.mods = []

        if not self.mass:
            if value in AA_mass:
                self.mass = AA_mass[value]

    def set_modification(self, i: Modification):
        """
        Add a modification to the list of modifications for this amino acid block.

        :param i: Modification
            The modification to be added.
        """
        self.mods.append(i)

    def __repr__(self):
        """
        Return a string representation of the amino acid block for debugging.

        :return: str
            The name of the amino acid with its modifications.
        """
        s = self.value
        for i in self.mods:
            s += "[{}]".format(i.value)
        return s

    def __str__(self):
        """
        Return a string representation of the amino acid block.

        :return: str
            The name of the amino acid with its modifications.
        """
        s = self.value
        for i in self.mods:
            s += "[{}]".format(i.value)
        return s
