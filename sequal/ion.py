from sequal.sequence import Sequence
from sequal.mass import calculate_mass

modifier = {
    "b": -18-19,
}

class Ion(Sequence):
    def __init__(self, seq, charge=1, ion_type=None, fragment_number=None):
        super().__init__(seq)
        self.charge = charge
        self.ion_type = ion_type
        self.fragment_number = fragment_number

    def mass_calculate(self):
        m = calculate_mass(self.seq)
        ion = modifier.get(self.ion_type, 0)
        mi = (m + ion + self.charge)/self.charge
        return mi

