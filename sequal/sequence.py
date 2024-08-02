
"""
This module provides the Sequence class for handling peptide or protein sequences and their fragments.

Classes:
    Sequence: Represents a sequence of amino acids or modifications.
    ModdedSequenceGenerator: Generates modified sequences based on static and variable modifications.

Functions:
    count_unique_elements(seq): Counts unique elements in a sequence.
    variable_position_placement_generator(positions): Generates different position combinations for modifications.
    ordered_serialize_position_dict(positions): Serializes a dictionary of positions in an ordered manner.
"""
import re
from typing import Set, Any, List

from sequal.amino_acid import AminoAcid
from sequal.modification import Modification, ModificationMap
from copy import deepcopy
import itertools
from json import dumps

mod_pattern = re.compile(r"[\(|\[]+([^\)]+)[\)|\]]+")
mod_enclosure_start = {"(", "[", "{"}
mod_enclosure_end = {")", "]", "}"}

# Base sequence object for peptide or protein sequences and their fragments
class Sequence:
    """
    Represents a sequence of amino acids or modifications.

    :param seq: iterable
        A string or array of strings or array of AminoAcid objects. The parser will recursively look over each string at
        the deepest level and identify individual modifications or amino acids for processing.
    :param encoder: BaseBlock, optional
        Class for encoding the sequence (default is AminoAcid).
    :param mods: dict, optional
        Dictionary whose keys are the positions within the sequence and values are arrays of modifications at those positions (default is None).
    :param parse: bool, optional
        Whether to parse the sequence (default is True).
    :param parser_ignore: list, optional
        List of items to ignore during parsing (default is None).
    :param mod_position: str, optional
        Indicates the position of the modifications relative to the base block it is supposed to modify (default is "right").
    """

    seq: List[Any]

    def __init__(self, seq, encoder=AminoAcid, mods=None, parse=True, parser_ignore=None, mod_position="right"):
        if type(seq) is not Sequence:
            if not mods:
                self.mods = {}
            else:
                self.mods = mods
            self.encoder = encoder
            if not parser_ignore:
                self.parser_ignore = []
            else:
                self.parser_ignore = parser_ignore
            self.seq = []
            current_mod = []
            current_position = 0
            if parse:
                self.sequence_parse(current_mod, current_position, mod_position, mods, seq)

        else:
            for k in seq.__dict__:
                if k != "mods":
                    setattr(self, k, deepcopy(seq.__dict__[k]))
        self.seq_length = len(self.seq)

    def __getitem__(self, key):
        return self.seq[key]

    def __len__(self):
        return self.seq_length

    def __repr__(self):
        a = ""
        for i in self.seq:
            a += str(i)
        return a

    def __str__(self):
        a = ""
        for i in self.seq:
            a += str(i)
        return a

    def sequence_parse(self, current_mod, current_position, mod_position, mods, seq):
        """
        Parse the sequence input.

        :param current_mod: list
            List of current modifications.
        :param current_position: int
            Current iterating amino acid position from the input sequence.
        :param mod_position: str
            Modification position relative to the modified residue.
        :param mods: dict
            External modification input.
        :param seq: iterable
            Sequence input.
        """
        for b, m in self.__load_sequence_iter(iter(seq)):
            if not m:
                if mod_position == "left":
                    if type(b) == AminoAcid:
                        current_unit = b
                        current_unit.position = current_position
                    else:
                        current_unit = self.encoder(b, current_position)

                    if current_mod and not mods:
                        for i in current_mod:
                            current_unit.set_modification(i)
                    elif current_position in self.mods and current_unit:
                        if type(self.mods[current_position]) == Modification:
                            current_unit.set_modification(self.mods[current_position])
                        else:
                            for mod in self.mods[current_position]:
                                current_unit.set_modification(mod)

                    self.seq.append(deepcopy(current_unit))

                    current_mod = []
                if mod_position == "right":

                    if current_mod and not mods:
                        for i in current_mod:
                            self.seq[current_position - 1].set_modification(i)
                    if type(b) == AminoAcid:
                        current_unit = b
                        current_unit.position = current_position
                    else:
                        current_unit = self.encoder(b, current_position)

                    if current_position in self.mods and current_unit:
                        if type(self.mods[current_position]) == Modification:
                            current_unit.set_modification(self.mods[current_position])

                        else:
                            for mod in self.mods[current_position]:
                                current_unit.set_modification(mod)

                    self.seq.append(deepcopy(current_unit))

                    current_mod = []
                current_position += 1
            else:
                if not mods:
                    # current_mod.append(Modification(b[1:-1]))
                    if mod_position == "right":
                        self.seq[current_position-1]\
                            .set_modification(Modification(b[1:-1]))
                    else:
                        current_mod.append(Modification(b[1:-1]))

    def __load_sequence_iter(self, seq=None, iter_seq=None):
        """
        Load the sequence iterator.

        :param seq: iterable, optional
            Sequence input.
        :param iter_seq: iterator, optional
            Sequence iterator.
        :yield: tuple
            A tuple containing the block and modification status.
        """
        mod_open = 0
        block = ""
        mod = False
        if not iter_seq:
            iter_seq = iter(seq)
        for i in iter_seq:
            if type(i) == str:
                if i in mod_enclosure_start:
                    mod = True
                    mod_open += 1
                elif i in mod_enclosure_end:
                    mod_open -= 1
                block += i
            elif type(i) == AminoAcid:
                block = i
            else:
                yield from self.__load_sequence_iter(iter_seq=iter_seq)
            if mod_open == 0:
                yield (block, mod)
                mod = False
                block = ""

    def __iter__(self):
        self.current_iter_count = 0
        return self

    def __next__(self):
        if self.current_iter_count == self.seq_length:
            raise StopIteration
        result = self.seq[self.current_iter_count]
        self.current_iter_count += 1
        return result

    def add_modifications(self, mod_dict):
        """
        Add modifications to the sequence.

        :param mod_dict: dict
            Dictionary of modifications to add.
        """
        for aa in self.seq:
            if aa.position in mod_dict:
                for mod in mod_dict[aa.position]:
                    aa.set_modification(mod)

    def to_stripped_string(self):
        """
        Return the sequence as a string without any modification annotations.

        :return: str
            The stripped sequence string.
        """
        seq = ""
        for i in self.seq:
            seq += i.value
        return seq

    def to_string_customize(self, data, annotation_placement="right", block_separator="", annotation_enclose_characters=("[", "]"),
                            individual_annotation_enclose=False, individual_annotation_enclose_characters=("[", "]"),
                            individual_annotation_separator=""):
        """
        Customize the sequence string with annotations.

        :param data: dict
            A dictionary where the key is the index position of the amino acid residue and the value is an iterable containing the items to be included in the sequence.
        :param annotation_placement: str, optional
            Whether the information should be included on the right or left of the residue (default is "right").
        :param block_separator: str, optional
            Separator between each block of annotation information (default is "").
        :param annotation_enclose_characters: tuple, optional
            Enclosure characters for each annotation cluster (default is ("[", "]")).
        :param individual_annotation_enclose: bool, optional
            Whether each individual annotation should be enclosed (default is False).
        :param individual_annotation_enclose_characters: tuple, optional
            Enclosure characters for each individual annotation (default is ("[", "]")).
        :param individual_annotation_separator: str, optional
            Separator for each individual annotation (default is "").
        :return: str
            The customized sequence string.
        """
        assert annotation_placement in {"left", "right"}
        seq = []
        for i in range(len(self.seq)):
            seq.append(self.seq[i].value)
            if i in data:
                annotation = []
                if individual_annotation_enclose:
                    for v in data[i]:
                        annotation.append("{}{}{}".format(individual_annotation_enclose_characters[0], v, individual_annotation_enclose_characters[1]))
                else:
                    annotation = data[i]
                if type(annotation) == str:
                    ann = annotation
                else:
                    ann = individual_annotation_separator.join(annotation)
                if annotation_enclose_characters:
                    seq.append("{}{}{}".format(annotation_enclose_characters[0], ann, annotation_enclose_characters[1]))
                else:
                    seq.append(individual_annotation_separator.join(ann))
        return block_separator.join(seq)

    def find_with_regex(self, motif, ignore=None):
        """
        Find positions in the sequence that match a given regex motif.

        :param motif: str
            The regex pattern to search for in the sequence.
        :param ignore: list of bool, optional
            A list indicating positions to ignore in the sequence. If provided, positions corresponding to True values will be ignored (default is None).

        :yield: slice
            A slice object representing the start and end positions of each match in the sequence.
        """
        pattern = re.compile(motif)
        new_str = ""
        if ignore is not None:
            for i in range(len(ignore)):
                if not ignore[i]:
                    new_str += self.seq[i].value
        else:
            new_str = self.to_stripped_string()

        for i in pattern.finditer(new_str):

            if not i.groups():
                yield slice(i.start(), i.end())
            else:
                for m in range(1, len(i.groups()) + 1):
                    yield slice(i.start(m), i.end(m))

    def gaps(self):
        """
        Identify gaps in the sequence.

        This method returns a list of boolean values indicating the presence of gaps in the sequence. A gap is represented by a '-' character.

        :return: list of bool
            A list where each element is True if the corresponding position in the sequence is a gap, and False otherwise.
        """
        s = [False for i in range(len(self.seq))]
        for i in range(len(s)):
            if self.seq[i].value == '-':
                s[i] = True

        return s

    def count(self, char, start, end):
        """
        Count the occurrences of a character in the sequence within a specified range.

        :param char: str
            The character to count in the sequence.
        :param start: int
            The starting index of the range to count the character.
        :param end: int
            The ending index of the range to count the character.

        :return: int
            The number of occurrences of the character in the specified range.
        """
        return self.to_stripped_string().count(char, start, end)

def count_unique_elements(seq):
    """
    Count unique elements in a sequence.

    This function iterates through the sequence and counts the occurrences of each unique element, including modifications.

    :param seq: iterable
        The sequence to count unique elements from. Each element should have a `value` attribute and optionally a `mods` attribute.
    :return: dict
        A dictionary where keys are unique element values and values are their counts.
    """
    elements = {}
    for i in seq:
        if i.value not in elements:
            elements[i.value] = 0
        elements[i.value] += 1
        if i.mods:
            for m in i.mods:
                if m.value not in elements:
                    elements[m.value] = 0
                elements[m.value] += 1
    return elements


def variable_position_placement_generator(positions):
    """
    Generate different position combinations for modifications.

    This function uses `itertools.product` to generate a list of tuples with different combinations of 0 and 1.
    The length of each tuple is the same as the length of the input positions. Using `itertools.compress`,
    for each output from `itertools.product` paired with input positions, it generates a list of positions
    where only those with the same index as 1 are yielded.

    :param positions: list
        A list of all identified positions for the modification on the sequence.
    :yield: list
        A list of positions for each combination.
    """
    for i in itertools.product([0, 1], repeat=len(positions)):
        yield list(itertools.compress(positions, i))


def ordered_serialize_position_dict(positions):
    """
    Serialize a dictionary of positions in an ordered manner.

    This function serializes the input dictionary of positions into a JSON string, ensuring the keys are sorted.

    :param positions: dict
        The dictionary of positions to serialize.
    :return: str
        The serialized JSON string of the positions dictionary.
    """
    return dumps(positions, sort_keys=True, default=str)


class ModdedSequenceGenerator:
    """
    Generator for creating modified sequences.

    This class generates all possible sequences with static and variable modifications applied to a base sequence.

    :param seq: str
        The base sequence to be modified.
    :param variable_mods: list of Modification, optional
        List of variable modifications to apply (default is None).
    :param static_mods: list of Modification, optional
        List of static modifications to apply (default is None).
    :param used_scenarios: set, optional
        Set of already used modification scenarios to avoid duplicates (default is None).
    :param parse_mod_position: bool, optional
        Whether to parse positions of modifications in the sequence (default is True).
    :param mod_position_dict: dict, optional
        Dictionary of modification positions (default is None).
    :param ignore_position: set, optional
        Set of positions to ignore when applying modifications (default is None).
    """
    used_scenarios_set: Set[str]

    def __init__(self, seq, variable_mods=None, static_mods=None, used_scenarios=None, parse_mod_position=True, mod_position_dict=None, ignore_position=None):
        self.seq = seq
        if static_mods:
            self.static_mods = static_mods

            self.static_map = ModificationMap(seq, static_mods, parse_position=parse_mod_position, mod_position_dict=mod_position_dict)
            self.static_mod_position_dict = self.static_mod_generate()
        else:
            self.static_mod_position_dict = {}
        if ignore_position:
            self.ignore_position = ignore_position
        else:
            self.ignore_position = set()

        for i in self.static_mod_position_dict:
            self.ignore_position.add(i)

        if variable_mods:
            self.variable_mods = variable_mods
            if self.static_mod_position_dict:
                self.variable_map = ModificationMap(seq, variable_mods, ignore_positions=self.ignore_position, parse_position=parse_mod_position, mod_position_dict=mod_position_dict)
            else:
                self.variable_map = ModificationMap(seq, variable_mods)
            self.variable_mod_number = len(variable_mods)
        else:
            self.variable_mods = None

        self.variable_map_scenarios = {}
        if used_scenarios:
            self.used_scenarios_set = used_scenarios
        else:
            self.used_scenarios_set = set()

    def generate(self):
        """
        Generate all possible modified sequences.

        This method yields dictionaries representing different modification scenarios applied to the base sequence.

        :yield: dict
            A dictionary where keys are positions and values are lists of modifications at those positions.
        """
        if self.variable_mods:
            self.variable_mod_generate_scenarios()
            for i in self.explore_scenarios():
                a = dict(self.static_mod_position_dict)
                a.update(i)
                serialized_a = ordered_serialize_position_dict(a)
                if serialized_a not in self.used_scenarios_set:
                    self.used_scenarios_set.add(serialized_a)
                    yield a
        else:
            serialized_a = ordered_serialize_position_dict(self.static_mod_position_dict)
            if serialized_a not in self.used_scenarios_set:
                yield self.static_mod_position_dict

    def static_mod_generate(self):
        """
        Generate positions for static modifications.

        This method creates a dictionary of positions for static modifications in the sequence.

        :return: dict
          A dictionary where keys are positions and values are lists of static modifications at those positions.
        """
        position_dict = {}

        for m in self.static_mods:

            for pm in self.static_map.get_mod_positions(str(m)):
                if pm not in position_dict:
                    position_dict[pm] = []
                position_dict[pm].append(m)
        return position_dict

    def variable_mod_generate_scenarios(self):
        """
        Generate all possible position combinations for variable modifications.

        This method populates the `variable_map_scenarios` dictionary with all possible position combinations for each variable modification.
        """
        for i in self.variable_mods:
            positions = self.variable_map.get_mod_positions(str(i))
            if i.value not in self.variable_map_scenarios:
                if not i.all_fill:
                    self.variable_map_scenarios[i.value] = list(
                        variable_position_placement_generator(positions))
                else:
                    self.variable_map_scenarios[i.value] = [[], positions]


    def explore_scenarios(self, current_mod=0, mod=None):
        """
        Recursively explore all modification scenarios.

        This method recursively generates all possible modification scenarios by exploring different combinations of variable modifications.

        :param current_mod: int, optional
            The current modification index being explored (default is 0).
        :param mod: dict, optional
            The current modification scenario being built (default is None).

        :yield: dict
            A dictionary representing a modification scenario.
        """
        if mod is None:
            mod = {}
        for pos in self.variable_map_scenarios[self.variable_mods[current_mod].value]:
            temp_dict = deepcopy(mod)
            if pos:
                for p in pos:
                    if p not in temp_dict:
                        temp_dict[p] = [self.variable_mods[current_mod]]
                    if current_mod != self.variable_mod_number - 1:
                        yield from self.explore_scenarios(current_mod + 1, temp_dict)
                    else:
                        yield temp_dict
            else:
                if current_mod != self.variable_mod_number - 1:
                    yield from self.explore_scenarios(current_mod + 1, temp_dict)
                else:
                    yield temp_dict


