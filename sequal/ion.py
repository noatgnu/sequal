from sequal.mass import calculate_mass
from sequal.resources import proton
from sequal.sequence import Sequence

modifier = {
    "b": -18 - 19,
}


class Ion(Sequence):
    """
    Represents an ion fragment sequence, inheriting from the Sequence class.

    This class converts a Sequence object into an ion fragment, adding properties
    such as charge, ion type, and fragment number. It also handles modifications
    and labile groups within the sequence.

    Attributes:
        charge (int): The charge of the ion.
        ion_type (str): The name of the transition type (e.g., 'b', 'y').
        fragment_number (int): The number of the transition.
        mods (dict): A dictionary of modifications, with keys as positions.
        has_labile (bool): True if the ion has labile modifications.
    """

    def __init__(self, seq, charge=1, ion_type=None, fragment_number=None):
        """
        Initializes an Ion object.

        Args:
            seq (Sequence): The Sequence object to be converted into an Ion fragment.
            charge (int, optional): The charge of the ion. Defaults to 1.
            ion_type (str, optional): The name of the transition type. Defaults to None.
            fragment_number (int, optional): The number of the transition. Defaults to None.
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
        Calculates the mass-to-charge ratio (m/z) of the ion.

        Args:
            charge (int, optional): The charge of the ion. If not specified, the
                object's charge is used. Defaults to None.
            with_water (bool, optional): Whether the mass will be calculated with
                or without water. Defaults to False.
            extra_mass (float, optional): Extra modification of mass that is not
                represented within the sequence. Defaults to 0.

        Returns:
            float: The calculated m/z value of the ion.
        """
        if not charge:
            charge = self.charge
        m = calculate_mass(self.seq, with_water=with_water) + extra_mass

        # Charge is calculated with the hardcoded mass of protons
        mi = (m + charge * proton) / charge
        return mi
