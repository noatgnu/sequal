"""
Provides functionality for fragmenting sequences for mass spectrometry analysis.

This module contains a factory class for generating ion fragments from sequences,
along with functions for handling labile and non-labile modifications during
fragmentation.
"""

from sequal.ion import Ion

ax = "ax"
by = "by"
cz = "cz"


def fragment_non_labile(sequence, fragment_type):
    """
    Calculates non-labile modifications and yields associated transitions.

    For example, "by" would yield a tuple of "b" and "y" transitions.

    Args:
        sequence (sequal.sequence.Sequence): The sequence to be fragmented.
        fragment_type (str): The type of fragment transition (e.g., "by", "ax").

    Yields:
        tuple: A tuple containing the left and right ion fragments.
    """
    for i in range(1, sequence.seq_length, 1):
        left = Ion(sequence[:i], fragment_number=i, ion_type=fragment_type[0])
        right = Ion(
            sequence[i:],
            fragment_number=sequence.seq_length - i,
            ion_type=fragment_type[1],
        )
        yield left, right


def fragment_labile(sequence):
    """
    Calculates all labile modification variants for the sequence.

    Args:
        sequence (sequal.sequence.Sequence): The sequence to be fragmented.

    Returns:
        Ion: An Ion object representing the fragmented sequence with labile
             modifications.
    """
    fragment_number = 0
    for p in sequence.mods:
        for i in sequence.mods[p]:
            if i.labile:
                fragment_number += i.labile_number
    return Ion(sequence, fragment_number=fragment_number, ion_type="Y")


class FragmentFactory:
    """
    A factory class for generating ion fragments from sequences.

    Attributes:
        fragment_type (str): The type of fragment transition (e.g., "by", "ax").
        ignore (list): A list of modifications to ignore.
    """

    def __init__(self, fragment_type, ignore=None):
        """
        Initializes a FragmentFactory object.

        Args:
            fragment_type (str): The type of fragment transition (e.g., "by", "ax").
            ignore (list, optional): A list of modifications to ignore.
                Defaults to None.
        """
        self.fragment_type = fragment_type
        if ignore:
            self.ignore = ignore
        else:
            self.ignore = []

    def set_ignore(self, ignore):
        """
        Sets the list of modifications to ignore.

        Args:
            ignore (list): A list of modifications to ignore.
        """
        self.ignore = ignore
