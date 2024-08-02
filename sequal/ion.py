from sequal.sequence import Sequence
from sequal.mass import calculate_mass
from sequal.resources import proton

modifier = {
    "b": -18-19,
}

class Ion(Sequence):
    """
    Represents an ion fragment sequence object, inheriting properties from the Sequence class.

    This class is used to convert a Sequence object into an Ion fragment, with additional properties
    such as charge, ion type, and fragment number. It also handles modifications and labile groups
    within the sequence.

    :param seq: Sequence
        The Sequence object to be converted into an Ion fragment.
    :param charge: int, optional
        The charge of the ion (default is 1).
    :param ion_type: str, optional
        The name of the transition type (default is None).
    :param fragment_number: int, optional
        The number of the transition (default is None).
    """

    def __init__(self, seq, charge=1, ion_type=None, fragment_number=None):
        """
        Initialize an Ion object.

        :param seq: Sequence
            The Sequence object to be converted into an Ion fragment.
        :param charge: int, optional
            The charge of the ion (default is 1).
        :param ion_type: str, optional
            The name of the transition type (default is None).
        :param fragment_number: int, optional
            The number of the transition (default is None).
        """
        super().__init__(seq)
        self.charge = charge
        self.ion_type = ion_type
        self.fragment_number = fragment_number
        self.mods = {}
        self.has_labile = False
        # Iterating through each amino acid position and build a modification list for the ion
        for i, aa in enumerate(self.seq):
            for m in aa.mods:
                if i not in self.mods:
                    self.mods[i] = []
                self.mods[i].append(m)
                if m.labile:
                    self.has_labile = True

    def mz_calculate(self, charge=None, with_water=False, extra_mass=0):
        """
        Calculate the mass-to-charge ratio (m/z) of the ion.

        :param charge: int, optional
            The charge of the ion. If not specified, the object's charge is used.
        :param with_water: bool, optional
            Whether the mass will be calculated with or without water (default is False).
        :param extra_mass: float, optional
            Extra modification of mass that is not represented within the sequence (default is 0).

        :return: float
            The calculated m/z value of the ion.
        """
        if not charge:
            charge = self.charge
        m = calculate_mass(self.seq, with_water=with_water) + extra_mass

        # Charge is calculated with the hardcoded mass of protons
        mi = (m + charge * proton) / charge
        return mi


