from sequence import Sequence


class Ion(Sequence):
    def __init__(self, seq, charge=None, ion_type=None, fragment_number=None):
        super().__init__(seq)
        self.charge = charge
        self.ion_type = ion_type
        self.fragment_number = fragment_number