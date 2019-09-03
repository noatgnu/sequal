import re
from amino_acid import AminoAcid
from modification import Modification, ModificationMap
from copy import deepcopy
import itertools
from json import dumps
mod_pattern = re.compile(r"[\(|\[]+([^\)]+)[\)|\]]+")
mod_enclosure_start = {"(", "[", "{"}
mod_enclosure_end = {")", "]", "}"}


class Sequence:
    def __init__(self, seq, encoder=AminoAcid, mods=None, parse=True, parser_ignore=None, mod_position="right"):
        """
        :param mod_position
        Indicate the position of the modifications relative to the base block it is supposed to modify
        :type mod_position: str
        :param mods
        Dictionary whose keys are the positions within the sequence and values are array of modifications at those
        positions
        :type mods: dict
        :param encoder
        Class for encoding of sequence.
        :type encoder: BaseBlock
        :param seq
        String or array of strings or array of AminoAcid objects. The parser will recursively look over each string at
        deepest level and identify individual modifications or amino acids for processing
        :type seq: iterable
        Python iterable where the deepest level is a string
            
        """
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
        :param seq: sequence input
        :param mods: external modification input
        :param mod_position: modification position relative to the modified residue
        :param current_position: current iterating amino acid position from the input sequence
        :type current_mod: Modification
        """
        for b, m in self.__load_sequence_iter(iter(seq)):
            if not m:
                if mod_position == "left":
                    if type(b) == AminoAcid:
                        current_unit = b
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
                    current_mod.append(Modification(b))

    def __load_sequence_iter(self, seq=None, iter_seq=None):
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
        for aa in self.seq:
            if aa.position in mod_dict:
                for mod in mod_dict[aa.position]:
                    aa.set_modification(mod)

    def to_stripped_string(self):
        seq = ""
        for i in self.seq:
            seq += i.value
        return seq


def count_unique_elements(seq):
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
    for i in itertools.product([0, 1], repeat=len(positions)):
        yield list(itertools.compress(positions, i))


def ordered_serialize_position_dict(positions):
    return dumps(positions, sort_keys=True, default=str)


class ModdedSequenceGenerator:
    def __init__(self, seq, variable_mods=None, static_mods=None, used_scenarios=None):
        self.seq = seq
        if static_mods:
            self.static_mods = static_mods
            self.static_map = ModificationMap(seq, static_mods)
            self.static_mod_position_dict = self.static_mod_generate()
        else:
            self.static_mod_position_dict = {}

        if variable_mods:
            self.variable_mods = variable_mods
            if self.static_mod_position_dict:
                self.variable_map = ModificationMap(seq, variable_mods, ignore_positions=set(self.static_mod_position_dict.keys()))
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
        position_dict = {}
        for m in self.static_mods:
            for pm in self.static_map.get_mod_positions(m.value):
                if pm not in position_dict:
                    position_dict[pm] = []
                position_dict[pm].append(m)
        return position_dict

    def variable_mod_generate_scenarios(self):
        for i in self.variable_mods:
            positions = self.variable_map.get_mod_positions(i.value)
            if i.value not in self.variable_map_scenarios:
                self.variable_map_scenarios[i.value] = list(
                    variable_position_placement_generator(positions))

    def explore_scenarios(self, current_mod=0, mod=None):
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


