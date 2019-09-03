import re

from base_block import BaseBlock


class Modification(BaseBlock):
    def __init__(self, value, position=None, regex_pattern=None, full_name=None, mod_type="static", labile=False, labile_number=0, mass=None):
        """
        :param position
        Position of the modification on the block it belongs to. Should be int and not None if it is assigned to a
        block. None when using as unassigned.
        :type position: int
        """
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

    def __repr__(self):
        return self.value

    def __str__(self):
        return self.value

    def find_positions(self, seq):
        for i in self.regex.finditer(seq):
            res = len(i.groups())
            if res > 0:
                for r in range(0, res + 1):
                    yield i.start(r), i.end(r)
            else:
                yield i.start(), i.end()


class ModificationMap:
    def __init__(self, seq, mods, ignore_positions=None):
        self.ignore_positions = ignore_positions
        self.seq = seq
        self.mod_dict_by_name = {}
        self.mod_position_dict = {}
        for m in mods:
            d = []
            self.mod_dict_by_name[m.value] = m
            for p_start, p_end in m.find_positions(self.seq):
                if ignore_positions:
                    if p_start not in ignore_positions:
                        d.append(p_start)
                else:
                    d.append(p_start)
            self.mod_position_dict[m.value] = d

    def get_mod_positions(self, mod_name):
        if mod_name in self.mod_position_dict:
            return self.mod_position_dict[mod_name]
        else:
            return None

    def get_mod(self, mod_name):
        if mod_name in self.mod_dict_by_name:
            return self.mod_dict_by_name[mod_name]
        else:
            return None


