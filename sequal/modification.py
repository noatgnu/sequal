"""
This module provides classes and functions for handling modifications in sequences.

Classes:
    Modification: Represents a modification block with various properties such as mass, type, and position.
    ModificationMap: Internal class that carries a dictionary of modifications, their positions, and names for quick reference.

Functions:
    None
"""
import re

from sequal.base_block import BaseBlock


# Base modification block.
class Modification(BaseBlock):
    """
    Represents a modification block with various properties such as mass, type, and position.

    :param value: str
        Short name of the modification.
    :param position: int, optional
        Position of the modification on the block it belongs to. Should be int and not None if it is assigned to a block. None when using as unassigned.
    :param regex_pattern: str, optional
        Regular expression pattern that can be used to represent where the modification is found.
    :param full_name: str, optional
        Full name of the modification.
    :param mod_type: str, optional
        Modification type. Either "static" or "variable" (default is "static").
    :param labile: bool, optional
        Whether or not the modification is labile. Important for mass spectrometry fragmentation (default is False).
    :param labile_number: int, optional
        The order of the fragment in a labile fragmentation event (default is 0).
    :param mass: float, optional
        Mass of the modification (default is 0).
    :param all_filled: bool, optional
        Whether or not the modification could be always found at all the expected sites (default is False).
    """
    def __init__(self, value, position=None, regex_pattern=None, full_name=None, mod_type="static", labile=False, labile_number=0, mass=0, all_filled=False):
        super().__init__(value, position=position, branch=True, mass=mass)
        if regex_pattern:
            self.regex = re.compile(regex_pattern)
        else:
            self.regex = None
        if mod_type in {"static", "variable"}:
            self.mod_type = mod_type
        else:
            print("Type can only be 'static' or 'variable'")
            raise ValueError

        assert(type(labile) == bool)
        self.labile = labile
        self.labile_number = labile_number
        self.full_name = full_name
        self.all_fill = all_filled

    def __repr__(self):
        """
        Return a string representation of the modification.

        :return: str
            The short name of the modification, with labile number if applicable.
        """
        if not self.labile:
            return self.value
        else:
            return self.value + str(self.labile_number)

    def __str__(self):
        """
        Return a string representation of the modification.

        :return: str
            The short name of the modification, with labile number if applicable.
        """
        if not self.labile:
            return self.value
        else:
            return self.value + str(self.labile_number)

    def find_positions(self, seq):
        """
        Find positions of the modification in the given sequence using the regex pattern.

        :param seq: str
            The sequence to search for the modification.

        :yield: tuple
            A tuple containing the start and end positions of the modification in the sequence.
        """
        for i in self.regex.finditer(seq):
            res = len(i.groups())
            if res > 0:
                for r in range(0, res + 1):
                    yield i.start(r), i.end(r)
            else:
                yield i.start(), i.end()


# Internal used object that carry dictionary of modification, modification position and their name for quick reference
class ModificationMap:
    """
    Internal class that carries a dictionary of modifications, their positions, and names for quick reference.

    :param seq: str
        The sequence to be analyzed.
    :param mods: list
        A list of Modification objects to be mapped.
    :param ignore_positions: list, optional
        A list of positions to ignore (default is None).
    :param parse_position: bool, optional
        Whether to parse positions of modifications in the sequence (default is True).
    :param mod_position_dict: dict, optional
        A dictionary of modification positions (default is None).
    """
    def __init__(self, seq, mods, ignore_positions=None, parse_position=True, mod_position_dict=None):
        self.ignore_positions = ignore_positions
        self.seq = seq
        self.mod_dict_by_name = {}
        if mod_position_dict:
            self.mod_position_dict = mod_position_dict
        else:
            self.mod_position_dict = {}

        for m in mods:
            self.mod_dict_by_name[str(m)] = m
            if parse_position:
                d = []
                for p_start, p_end in m.find_positions(self.seq):
                    if ignore_positions:
                        if p_start not in ignore_positions:
                            d.append(p_start)
                    else:
                        d.append(p_start)
                self.mod_position_dict[str(m)] = d
        # print(self.mod_position_dict)

    def get_mod_positions(self, mod_name):
        """
        Get the positions of a modification by its name.

        :param mod_name: str
            The name of the modification.

        :return: list
            A list of positions where the modification is found.
        """
        if mod_name in self.mod_position_dict:
            return self.mod_position_dict[mod_name]
        else:
            return None

    def get_mod(self, mod_name):
        """
        Get the Modification object by its name.

        :param mod_name: str
            The name of the modification.

        :return: Modification
            The Modification object.
        """
        if mod_name in self.mod_dict_by_name:
            return self.mod_dict_by_name[mod_name]
        else:
            return None


