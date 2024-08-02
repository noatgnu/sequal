"""
This module provides functionality for fragmenting sequences into ion fragments for mass spectrometry analysis.

Classes:
    FragmentFactory: A factory class for generating ion fragments from sequences.

Functions:
    fragment_non_labile(sequence, fragment_type): Calculate non-labile modifications and yield associated transitions.
    fragment_labile(sequence): Calculate all labile modification variants for the sequence and its associated labile modifications.
"""

from sequal.ion import Ion

ax = "ax"
by = "by"
cz = "cz"


def fragment_non_labile(sequence, fragment_type):
    """
    Calculate non-labile modifications and yield associated transitions.

    For example, "by" would yield a tuple of "b" and "y" transitions.

    :param sequence: sequal.sequence.Sequence
        The sequence to be fragmented.
    :param fragment_type: str
        The type of fragment transition (e.g., "by", "ax").

    :yield: tuple
        A tuple containing the left and right ion fragments.
    """
    for i in range(1, sequence.seq_length, 1):
        left = Ion(sequence[:i], fragment_number=i, ion_type=fragment_type[0])
        right = Ion(sequence[i:], fragment_number=sequence.seq_length - i, ion_type=fragment_type[1])
        yield left, right


def fragment_labile(sequence):
    """
    Calculate all labile modification variants for the sequence and its associated labile modifications.

    :param sequence: sequal.sequence.Sequence
        The sequence to be fragmented.

    :return: Ion
        An Ion object representing the fragmented sequence with labile modifications.
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

    :param fragment_type: str
        The type of fragment transition (e.g., "by", "ax").
    :param ignore: list, optional
        A list of modifications to ignore (default is None).
    """

    def __init__(self, fragment_type, ignore=None):
        """
        Initialize a FragmentFactory object.

        :param fragment_type: str
            The type of fragment transition (e.g., "by", "ax").
        :param ignore: list, optional
            A list of modifications to ignore (default is None).
        """
        self.fragment_type = fragment_type
        if ignore:
            self.ignore = ignore
        else:
            self.ignore = []

    def set_ignore(self, ignore):
        """
        Set the list of modifications to ignore.

        :param ignore: list
            A list of modifications to ignore.
        """
        self.ignore = ignore