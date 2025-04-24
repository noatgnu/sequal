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

import itertools
import json
import re
from collections import Counter, defaultdict
from copy import deepcopy
from typing import (
    Any,
    Dict,
    Generic,
    Iterable,
    Iterator,
    List,
    Optional,
    Set,
    Tuple,
    Type,
    TypeVar,
    Union,
)

from sequal.amino_acid import AminoAcid
from sequal.base_block import BaseBlock
from sequal.modification import GlobalModification, Modification, ModificationMap
from sequal.proforma import ProFormaParser, SequenceAmbiguity

mod_pattern = re.compile(r"[\(|\[]+([^\)]+)[\)|\]]+")
mod_enclosure_start = {"(", "[", "{"}
mod_enclosure_end = {")", "]", "}"}
T = TypeVar("T", bound="BaseBlock")


class Sequence:
    """
    Represents a sequence of amino acids or modifications.

    This class provides methods for building, modifying, and analyzing sequences
    of biochemical blocks such as amino acids with their modifications.

    Parameters
    ----------
    seq : Union[str, List, 'Sequence']
        A string, list of strings, or list of AminoAcid objects. The parser will
        recursively parse each element to identify amino acids and modifications.
    encoder : Type[BaseBlock], optional
        Class for encoding sequence elements, default is AminoAcid.
    mods : Optional[Dict[int, Union[Modification, List[Modification]]]], optional
        Dictionary mapping positions to modifications at those positions.
    parse : bool, optional
        Whether to parse the sequence, default is True.
    parser_ignore : Optional[List[str]], optional
        List of elements to ignore during parsing.
    mod_position : str, optional
        Indicates the position of modifications relative to their residue ("left" or "right").
    """

    # Regular expression patterns for parsing
    _MOD_PATTERN = re.compile(r"[\(|\[]+([^\)]+)[\)|\]]+")
    _MOD_ENCLOSURE_START = {"(", "[", "{"}
    _MOD_ENCLOSURE_END = {")", "]", "}"}

    def __init__(
        self,
        seq: Union[str, List, "Sequence"],
        encoder: Type[BaseBlock] = AminoAcid,
        mods: Optional[Dict[int, Union[Modification, List[Modification]]]] = None,
        parse: bool = True,
        parser_ignore: Optional[List[str]] = None,
        mod_position: str = "right",
        chains: Optional[List["Sequence"]] = None,
        global_mods: List[GlobalModification] = None,
        sequence_ambiguities: List[SequenceAmbiguity] = None,
        charge=None,
        ionic_species=None,
    ):
        self.encoder = encoder
        self.parser_ignore = parser_ignore or []
        self.seq: List[AminoAcid] = []
        self.chains = chains or [self]
        self.is_multi_chain = False
        self.mods = defaultdict(list)
        self.global_mods = global_mods or []
        self.sequence_ambiguities = sequence_ambiguities or []
        self.charge = charge
        self.ionic_species = ionic_species
        self.is_chimeric = False
        self.peptidoforms = []

        if isinstance(seq, Sequence):
            for attr_name, attr_value in seq.__dict__.items():
                if attr_name != "mods":
                    setattr(self, attr_name, deepcopy(attr_value))
            if hasattr(seq, "mods"):
                for pos, mod_list in seq.mods.items():
                    self.mods[pos] = deepcopy(mod_list)
        else:
            if mods:
                for pos, mod_items in mods.items():
                    if isinstance(mod_items, list):
                        self.mods[pos] = deepcopy(mod_items)
                    else:
                        self.mods[pos].append(mod_items)

            if parse:
                self._parse_sequence(seq, mod_position)
        self.seq_length = len(self.seq)

    @classmethod
    def from_proforma(self, proforma_str):
        """Create a Sequence object from a ProForma string with multi-chain support."""

        if "//" in proforma_str:
            chains = proforma_str.split("//")
            main_seq = self.from_proforma(chains[0])
            main_seq.is_multi_chain = True
            main_seq.chains = [main_seq]

            for chain_str in chains[1:]:
                chain = self.from_proforma(chain_str)
                main_seq.chains.append(chain)

            return main_seq

        peptidoforms = split_chimeric_proforma(proforma_str)

        if len(peptidoforms) > 1:
            # First peptidoform becomes the main sequence
            main_seq = self.from_proforma(peptidoforms[0])

            # Initialize chimeric properties if not already present
            main_seq.is_chimeric = True

            main_seq.peptidoforms = [main_seq]

            # Parse remaining peptidoforms
            for peptidoform_str in peptidoforms[1:]:
                peptidoform = self.from_proforma(peptidoform_str)
                peptidoform.is_chimeric = True
                main_seq.peptidoforms.append(peptidoform)

            return main_seq

        (
            base_sequence,
            modifications,
            global_mods,
            sequence_ambiquities,
            charge_info,
        ) = ProFormaParser.parse(proforma_str)
        charge = None
        ionic_species = None
        if charge_info:
            charge = charge_info[0]
            if len(charge_info) > 1:
                ionic_species = charge_info[1]
        seq = self(
            base_sequence,
            global_mods=global_mods,
            sequence_ambiguities=sequence_ambiquities,
            charge=charge,
            ionic_species=ionic_species,
        )
        for pos, mods in modifications.items():
            for mod in mods:
                if pos == -1:  # N-terminal
                    seq.mods[-1].append(mod)
                elif pos == -2:  # C-terminal
                    seq.mods[-2].append(mod)
                elif pos == -3:
                    seq.mods[-3].append(mod)
                elif pos == -4:
                    seq.mods[-4].append(mod)
                else:
                    seq.seq[pos].add_modification(mod)
        return seq

    def _chain_to_proforma(self, chain):
        result = ""
        for mod in self.global_mods:
            result += mod.to_proforma()

        ranges = []
        for i, aa in enumerate(self.seq):
            for mod in aa.mods:
                if (
                    mod.in_range
                    and mod.range_start is not None
                    and mod.range_end is not None
                ):
                    ranges.append((mod.range_start, mod.range_end, mod))
        ranges.sort(key=lambda x: x[0])
        if -4 in chain.mods:
            unknown_mods_by_value = defaultdict(int)
            for mod in chain.mods[-4]:
                unknown_mods_by_value[mod.to_proforma()] += 1
            for mod_value, count in unknown_mods_by_value.items():
                ambiguity_str = ""
                if count > 1:
                    ambiguity_str += f"[{mod_value}]"
                    ambiguity_str += f"^{count}?"
                else:
                    ambiguity_str = f"[{mod_value}]"
                    ambiguity_str += f"?"
                result += ambiguity_str
        if -3 in chain.mods:
            for mod in chain.mods[-3]:
                if mod.mod_type == "labile":
                    result += f"{{{mod.to_proforma()}}}"
        if -1 in chain.mods:
            n_mod_str = ""
            for mod in chain.mods[-1]:
                mod_str = f"[{mod.to_proforma()}]"
                n_mod_str += mod_str
            if n_mod_str:
                result += n_mod_str + "-"

        sorted_ambiguities = sorted(self.sequence_ambiguities, key=lambda a: a.position)
        ambiguity_index = 0

        for i, aa in enumerate(self.seq):
            if (
                ambiguity_index < len(sorted_ambiguities)
                and sorted_ambiguities[ambiguity_index].position == i
            ):
                ambiguity = sorted_ambiguities[ambiguity_index]
                result += f"(?{ambiguity.value})"
                ambiguity_index += 1

        range_start = False
        for i, aa in enumerate(chain.seq):
            if not range_start:
                for start, end, mod in ranges:
                    if i == start:
                        result += "("
                        range_start = True
                        break
            result += aa.value
            mod_str_data = ""
            if aa.mods:
                crosslink_refs_added = set()
                branch_refs_added = False
                for mod in aa.mods:
                    print(mod)
                    this_mod_str = mod.to_proforma()
                    print(this_mod_str)
                    if range_start and mod.in_range:
                        continue
                    if mod.mod_type == "ambiguous":
                        if not mod.has_ambiguity():
                            if mod.in_range:
                                this_mod_str = f"[{this_mod_str}]"
                            else:
                                this_mod_str = f"{{{this_mod_str}}}"
                        else:
                            this_mod_str = f"[{this_mod_str}]"
                    else:
                        this_mod_str = f"[{this_mod_str}]"

                    mod_str_data += this_mod_str
                result += mod_str_data

            if range_start:
                for start, end, mod in ranges:
                    if i == end:
                        result += ")"
                        range_start = False
                        break

                if not range_start:
                    result = self.get_mod_and_add_to_string(aa, result)

        if -2 in chain.mods:
            n_mod_str = ""
            for mod in chain.mods[-2]:
                mod_str = f"[{mod.to_proforma()}]"
                n_mod_str += mod_str
            if n_mod_str:
                result += "-" + n_mod_str
        if chain.charge:
            result += f"/{chain.charge}"
            if chain.ionic_species:
                result += f"[{chain.ionic_species}]"
        return result

    def _add_info_tags(self, result, mod):
        """
        Add info tags to the result string based on the modification.

        Parameters
        ----------
        result : str
            The current result string.
        mod : Modification
            The modification object.

        Returns
        -------
        str
            The updated result string with info tags.
        """
        info_str = ""
        added_info = set()
        stripped_result = result.lstrip("[")
        stripped_result = stripped_result.rstrip("]")
        added_info.add(stripped_result)
        for i in stripped_result.split(":", 1):
            i_splitted = i.split("#", 1)
            if len(i_splitted) > 1:
                for i2 in i_splitted:
                    if i2 and i2 not in added_info:
                        added_info.add(i2)
            added_info.add(i)
        if hasattr(mod, "mass") and mod.mass is not None:
            mass_str = f"+{mod.mass}" if mod.mass > 0 else f"{mod.mass}"
            if mass_str not in added_info:
                info_str += f"|{mass_str}"
                added_info.add(mass_str)

        if hasattr(mod, "synonyms") and mod.synonyms:
            for synonym in mod.synonyms:
                if synonym not in added_info:
                    info_str += f"|{synonym}"
                    added_info.add(synonym)

        if hasattr(mod, "observed_mass") and mod.observed_mass:
            obs_str = (
                f"Obs:+{mod.observed_mass}"
                if mod.observed_mass > 0
                else f"Obs:{mod.observed_mass}"
            )
            if obs_str not in added_info:
                info_str += f"|{obs_str}"
                added_info.add(obs_str)

        if hasattr(mod, "info_tags") and mod.info_tags:
            for tag in mod.info_tags:
                if tag not in added_info:
                    info_str += f"|{tag}"
                    added_info.add(tag)

        if info_str:
            if result.endswith("]"):
                result = result[:-1] + info_str + "]"
            else:
                result += info_str
        print(result)
        return result

    def get_mod_and_add_to_string(self, aa, result):
        mod_str = ""
        if aa.mods:
            for mod in aa.mods:
                this_mod_str = mod.to_proforma()
                if mod.mod_type == "ambiguous":
                    if not mod.has_ambiguity():
                        if mod.in_range:
                            this_mod_str = f"[{this_mod_str}]"
                        else:
                            this_mod_str = f"{{{this_mod_str}}}"
                    else:
                        this_mod_str = f"[{this_mod_str}]"
                else:
                    this_mod_str = f"[{this_mod_str}]"
                mod_str += this_mod_str
            result += mod_str
        return result

    def to_proforma(self):
        if hasattr(self, "is_multi_chain") and self.is_multi_chain:
            result = "//".join(self._chain_to_proforma(chain) for chain in self.chains)
            return result
        else:
            if hasattr(self, "is_chimeric") and self.is_chimeric:
                result = "+".join(
                    self._chain_to_proforma(peptide) for peptide in self.peptidoforms
                )
                return result
            else:
                return self._chain_to_proforma(self)

    def _parse_sequence(self, seq: Union[str, List], mod_position: str) -> None:
        """
        Parse the input sequence into a list of BaseBlock objects.

        Parameters
        ----------
        seq : Union[str, List]
            The input sequence to parse.
        mod_position : str
            Position of modifications relative to residues.
        """
        current_mod = []
        current_position = 0

        if mod_position not in {"left", "right"}:
            raise ValueError("mod_position must be either 'left' or 'right'")

        for block, is_mod in self._sequence_iterator(iter(seq)):
            if not is_mod:
                # Handle an amino acid/residue
                if mod_position == "left":
                    # Handle left-positioned modifications
                    if isinstance(block, self.encoder):
                        current_unit = block
                        current_unit.position = current_position
                    else:
                        current_unit = self.encoder(block, current_position)

                    # Apply pending modifications
                    self._apply_modifications(
                        current_unit, current_position, current_mod
                    )
                    self.seq.append(deepcopy(current_unit))
                    current_mod = []

                else:  # mod_position == "right"
                    # Apply modifications to previous residue
                    if current_mod and current_position > 0:
                        for mod in current_mod:
                            self.seq[current_position - 1].add_modification(mod)

                    # Create new residue
                    if isinstance(block, self.encoder):
                        current_unit = block
                        current_unit.position = current_position
                    else:
                        current_unit = self.encoder(block, current_position)

                    # Apply configured modifications
                    if current_position in self.mods:
                        mods_to_apply = self.mods[current_position]
                        if isinstance(mods_to_apply, list):
                            for mod in mods_to_apply:
                                current_unit.add_modification(mod)
                        else:
                            current_unit.add_modification(mods_to_apply)

                    self.seq.append(deepcopy(current_unit))
                    current_mod = []

                current_position += 1

            else:  # is_mod is True
                # Handle a modification
                if not self.mods:  # Only if not using predefined mods dict
                    # Extract mod string and create Modification object
                    mod_value = self._extract_mod_value(block)
                    mod_obj = Modification(mod_value)

                    if mod_position == "right" and current_position > 0:
                        # Apply directly to previous residue for right positioning
                        self.seq[current_position - 1].add_modification(mod_obj)
                    else:
                        # Store for later application with left positioning
                        current_mod.append(mod_obj)

    def _apply_modifications(
        self, block: BaseBlock, position: int, pending_mods: List[Modification]
    ) -> None:
        """
        Apply modifications to a block.

        Parameters
        ----------
        block : BaseBlock
            The block to modify.
        position : int
            Position of the block.
        pending_mods : List[Modification]
            List of pending modifications to apply.
        """
        # Apply pending modifications (from parsing)
        for mod in pending_mods:
            block.add_modification(mod)

        # Apply configured modifications (from mods dict)
        if position in self.mods:
            mods_to_apply = self.mods[position]
            if isinstance(mods_to_apply, list):
                for mod in mods_to_apply:
                    block.add_modification(mod)
            else:
                block.add_modification(mods_to_apply)

    def _extract_mod_value(self, mod_str: str) -> str:
        """
        Extract modification value from a string.

        Parameters
        ----------
        mod_str : str
            String containing modification.

        Returns
        -------
        str
            Extracted modification value.
        """
        # Find content between brackets/parentheses
        if (
            mod_str[0] in self._MOD_ENCLOSURE_START
            and mod_str[-1] in self._MOD_ENCLOSURE_END
        ):
            return mod_str[1:-1]
        return mod_str

    def _sequence_iterator(self, seq_iter: Iterator) -> Iterator[Tuple[Any, bool]]:
        """
        Iterate through sequence elements, identifying blocks and modifications.

        Parameters
        ----------
        seq_iter : Iterator
            Iterator over sequence elements.

        Yields
        ------
        Tuple[Any, bool]
            Tuple of (block, is_modification)
        """
        mod_open = 0
        block = ""
        is_mod = False

        for item in seq_iter:
            if isinstance(item, str):
                if item in self._MOD_ENCLOSURE_START:
                    is_mod = True
                    mod_open += 1
                elif item in self._MOD_ENCLOSURE_END:
                    mod_open -= 1
                block += item
            elif isinstance(item, self.encoder):
                block = item
            else:
                # Recursively handle nested iterables
                yield from self._sequence_iterator(iter(item))
                continue

            if mod_open == 0 and block:
                yield (block, is_mod)
                is_mod = False
                block = ""

    def __getitem__(self, key: Union[int, slice]) -> Union[BaseBlock, "Sequence"]:
        """
        Get item or slice from sequence.

        Parameters
        ----------
        key : Union[int, slice]
            Index or slice to retrieve.

        Returns
        -------
        Union[BaseBlock, Sequence]
            Single block or new Sequence containing the slice.
        """
        if isinstance(key, slice):
            new_seq = Sequence(self, parse=False)
            new_seq.seq = self.seq[key]
            new_seq.seq_length = len(new_seq.seq)
            return new_seq
        return self.seq[key]

    def __len__(self) -> int:
        """Return the length of the sequence."""
        return self.seq_length

    def __repr__(self) -> str:
        """Return a programmatic representation of the sequence."""
        return f"Sequence('{str(self)}')"

    def __str__(self) -> str:
        """Return a string representation of the sequence."""
        return "".join(str(block) for block in self.seq)

    def __iter__(self) -> "Sequence":
        """Initialize iteration over the sequence."""
        self.current_iter_count = 0
        return self

    def __next__(self) -> BaseBlock:
        """Get the next block in the sequence."""
        if self.current_iter_count == self.seq_length:
            raise StopIteration
        result = self.seq[self.current_iter_count]
        self.current_iter_count += 1
        return result

    def __eq__(self, other: object) -> bool:
        """Check if two sequences are equal."""
        if not isinstance(other, Sequence):
            return False
        if self.seq_length != other.seq_length:
            return False
        return all(a == b for a, b in zip(self.seq, other.seq))

    def add_modifications(self, mod_dict: Dict[int, List[Modification]]) -> None:
        """
        Add modifications to residues at specified positions.

        Parameters
        ----------
        mod_dict : Dict[int, List[Modification]]
            Dictionary mapping positions to lists of modifications.
        """
        for aa in self.seq:
            if aa.position in mod_dict:
                for mod in mod_dict[aa.position]:
                    aa.add_modification(mod)

    def to_stripped_string(self) -> str:
        """
        Return the sequence as a string without any modification annotations.

        Returns
        -------
        str
            The stripped sequence string.
        """
        return "".join(block.value for block in self.seq)

    def to_string_customize(
        self,
        data: Dict[int, Union[str, List[str]]],
        annotation_placement: str = "right",
        block_separator: str = "",
        annotation_enclose_characters: Tuple[str, str] = ("[", "]"),
        individual_annotation_enclose: bool = False,
        individual_annotation_enclose_characters: Tuple[str, str] = ("[", "]"),
        individual_annotation_separator: str = "",
    ) -> str:
        """
        Customize the sequence string with annotations.

        Parameters
        ----------
        data : Dict[int, Union[str, List[str]]]
            Dictionary mapping positions to annotations.
        annotation_placement : str, optional
            Placement of annotations ("left" or "right").
        block_separator : str, optional
            Separator between blocks.
        annotation_enclose_characters : Tuple[str, str], optional
            Characters to enclose annotation groups.
        individual_annotation_enclose : bool, optional
            Whether to enclose individual annotations.
        individual_annotation_enclose_characters : Tuple[str, str], optional
            Characters to enclose individual annotations.
        individual_annotation_separator : str, optional
            Separator between individual annotations.

        Returns
        -------
        str
            Customized sequence string with annotations.
        """
        if annotation_placement not in {"left", "right"}:
            raise ValueError("annotation_placement must be either 'left' or 'right'")

        elements = []

        for i in range(len(self.seq)):
            # Add annotation before residue if placement is left
            if annotation_placement == "left" and i in data:
                annotation = self._format_annotation(
                    data[i],
                    individual_annotation_enclose,
                    individual_annotation_enclose_characters,
                    individual_annotation_separator,
                    annotation_enclose_characters,
                )
                elements.append(annotation)

            # Add residue
            elements.append(self.seq[i].value)

            # Add annotation after residue if placement is right
            if annotation_placement == "right" and i in data:
                annotation = self._format_annotation(
                    data[i],
                    individual_annotation_enclose,
                    individual_annotation_enclose_characters,
                    individual_annotation_separator,
                    annotation_enclose_characters,
                )
                elements.append(annotation)

        return block_separator.join(elements)

    def _format_annotation(
        self,
        annotations: Union[str, List[str]],
        individual_enclose: bool,
        individual_enclose_chars: Tuple[str, str],
        separator: str,
        group_enclose_chars: Optional[Tuple[str, str]],
    ) -> str:
        """
        Format annotation strings.

        Parameters
        ----------
        annotations : Union[str, List[str]]
            Annotations to format.
        individual_enclose : bool
            Whether to enclose individual annotations.
        individual_enclose_chars : Tuple[str, str]
            Characters to enclose individual annotations.
        separator : str
            Separator between annotations.
        group_enclose_chars : Optional[Tuple[str, str]]
            Characters to enclose the entire annotation group.

        Returns
        -------
        str
            Formatted annotation string.
        """
        if isinstance(annotations, str):
            ann_text = annotations
        else:
            if individual_enclose:
                enclosed_annotations = [
                    f"{individual_enclose_chars[0]}{item}{individual_enclose_chars[1]}"
                    for item in annotations
                ]
                ann_text = separator.join(enclosed_annotations)
            else:
                ann_text = separator.join(annotations)

        if group_enclose_chars:
            return f"{group_enclose_chars[0]}{ann_text}{group_enclose_chars[1]}"
        return ann_text

    def find_with_regex(
        self, motif: str, ignore: Optional[List[bool]] = None
    ) -> Iterator[slice]:
        """
        Find positions in the sequence that match a given regex motif.

        Parameters
        ----------
        motif : str
            Regex pattern to search for.
        ignore : Optional[List[bool]], optional
            Positions to ignore (True = ignore).

        Yields
        ------
        slice
            Slice representing match position.
        """
        pattern = re.compile(motif)

        if ignore is not None:
            # Build string excluding ignored positions
            seq_str = "".join(
                self.seq[i].value for i in range(len(ignore)) if not ignore[i]
            )
        else:
            seq_str = self.to_stripped_string()

        for match in pattern.finditer(seq_str):
            if not match.groups():
                yield slice(match.start(), match.end())
            else:
                for group_idx in range(1, len(match.groups()) + 1):
                    yield slice(match.start(group_idx), match.end(group_idx))

    def gaps(self) -> List[bool]:
        """
        Identify gaps in the sequence.

        Returns
        -------
        List[bool]
            List where True indicates a gap at that position.
        """
        return [block.value == "-" for block in self.seq]

    def count(self, char: str, start: int, end: int) -> int:
        """
        Count occurrences of a character in a range.

        Parameters
        ----------
        char : str
            Character to count.
        start : int
            Start position.
        end : int
            End position.

        Returns
        -------
        int
            Number of occurrences.
        """
        return self.to_stripped_string().count(char, start, end)

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the sequence to a dictionary representation.

        Returns
        -------
        Dict[str, Any]
            Dictionary with sequence data.
        """
        # Collect all modifications by position
        mods_by_position = {}
        for i, aa in enumerate(self.seq):
            if hasattr(aa, "mods") and aa.mods:
                mods_by_position[i] = [mod.to_dict() for mod in aa.mods]

        return {
            "sequence": self.to_stripped_string(),
            "modifications": mods_by_position,
        }


def count_unique_elements(seq: Iterable[T]) -> Dict[str, int]:
    """
    Count unique elements in a sequence.

    Parameters
    ----------
    seq : Iterable[BaseBlock]
        The sequence of blocks to count elements from. Each element should
        have a `value` attribute and optionally a `mods` attribute.

    Returns
    -------
    Dict[str, int]
        Dictionary where keys are element values and values are their counts.

    Examples
    --------
    >>> from sequal.sequence import Sequence
    >>> seq = Sequence("ACDEFG")
    >>> counts = count_unique_elements(seq)
    >>> sorted(counts.items())
    [('A', 1), ('C', 1), ('D', 1), ('E', 1), ('F', 1), ('G', 1)]
    """
    if not seq:
        return {}

    elements = Counter()

    for item in seq:
        elements[item.value] += 1

        # Count modifications if present
        if hasattr(item, "mods") and item.mods:
            for mod in item.mods:
                elements[mod.value] += 1

    return dict(elements)


def variable_position_placement_generator(positions: List[int]) -> Iterator[List[int]]:
    """
    Generate all possible position combinations for modifications.

    This function creates all possible subsets of positions, including
    the empty set and the full set.

    Parameters
    ----------
    positions : List[int]
        List of positions where modifications could be applied.

    Yields
    ------
    List[int]
        Each possible combination of positions.

    Examples
    --------
    >>> list(variable_position_placement_generator([1, 2]))
    [[], [1], [2], [1, 2]]
    """
    if not positions:
        yield []
        return

    # Sort positions for consistent output
    sorted_positions = sorted(positions)

    # Generate all possible combinations (2^n possibilities)
    for mask in itertools.product([0, 1], repeat=len(sorted_positions)):
        yield [pos for pos, include in zip(sorted_positions, mask) if include]


def ordered_serialize_position_dict(positions: Dict[int, Any]) -> str:
    """
    Serialize a dictionary of positions with consistent ordering.

    Parameters
    ----------
    positions : Dict[int, Any]
        Dictionary mapping positions to modifications or other data.

    Returns
    -------
    str
        JSON string with sorted keys for consistent serialization.

    Raises
    ------
    TypeError
        If the dictionary contains values that cannot be serialized to JSON.
    """
    try:
        return json.dumps(positions, sort_keys=True, default=str)
    except TypeError as e:
        raise TypeError(f"Could not serialize positions dictionary: {e}") from e


class ModdedSequenceGenerator:
    """
    Generator for sequences with different modification combinations.

    This class creates all possible modified sequences by applying combinations
    of static and variable modifications to a base sequence.

    Parameters
    ----------
    seq : str
        The base sequence to modify.
    variable_mods : List[Modification], optional
        List of variable modifications to apply.
    static_mods : List[Modification], optional
        List of static modifications to apply.
    used_scenarios : Set[str], optional
        Set of serialized modification scenarios to avoid duplicates.
    parse_mod_position : bool, optional
        Whether to parse positions using modification regex patterns.
    mod_position_dict : Dict[str, List[int]], optional
        Pre-computed positions for modifications.
    ignore_position : Set[int], optional
        Set of positions to ignore when applying modifications.
    """

    def __init__(
        self,
        seq: str,
        variable_mods: Optional[List[Modification]] = None,
        static_mods: Optional[List[Modification]] = None,
        used_scenarios: Optional[Set[str]] = None,
        parse_mod_position: bool = True,
        mod_position_dict: Optional[Dict[str, List[int]]] = None,
        ignore_position: Optional[Set[int]] = None,
    ):
        self.seq = seq
        self.static_mods = static_mods or []
        self.variable_mods = variable_mods or []
        self.used_scenarios_set = used_scenarios or set()
        self.ignore_position = ignore_position or set()

        # Initialize modification maps and position dictionaries
        self.static_mod_position_dict = {}
        if self.static_mods:
            self.static_map = ModificationMap(
                seq,
                self.static_mods,
                parse_position=parse_mod_position,
                mod_position_dict=mod_position_dict,
            )
            self.static_mod_position_dict = self._generate_static_mod_positions()

            # Update ignore positions with static mod positions
            self.ignore_position.update(self.static_mod_position_dict.keys())

        self.variable_map_scenarios = {}
        if self.variable_mods:
            # Create variable modification map, considering ignored positions
            self.variable_map = ModificationMap(
                seq,
                self.variable_mods,
                ignore_positions=self.ignore_position,
                parse_position=parse_mod_position,
                mod_position_dict=mod_position_dict,
            )

    def generate(self) -> Iterator[Dict[int, List[Modification]]]:
        """
        Generate all possible modification combinations.

        Yields
        ------
        Dict[int, List[Modification]]
            Dictionary mapping positions to lists of modifications for each scenario.
        """
        if not self.variable_mods:
            # If only static mods, yield the single scenario if not already used
            serialized = ordered_serialize_position_dict(self.static_mod_position_dict)
            if serialized not in self.used_scenarios_set:
                self.used_scenarios_set.add(serialized)
                yield self.static_mod_position_dict
            return

        # Generate all variable mod scenarios
        self._generate_variable_mod_scenarios()

        # Explore all combinations
        for variable_scenario in self._explore_scenarios():
            # Combine static and variable modifications
            combined_scenario = dict(self.static_mod_position_dict)

            # Update with variable modifications
            for pos, mods in variable_scenario.items():
                if pos in combined_scenario:
                    combined_scenario[pos].extend(mods)
                else:
                    combined_scenario[pos] = mods

            # Serialize to check for duplicates
            serialized = ordered_serialize_position_dict(combined_scenario)
            if serialized not in self.used_scenarios_set:
                self.used_scenarios_set.add(serialized)
                yield combined_scenario

    def _generate_static_mod_positions(self) -> Dict[int, List[Modification]]:
        """
        Generate dictionary of positions for static modifications.

        Returns
        -------
        Dict[int, List[Modification]]
            Dictionary mapping positions to lists of static modifications.
        """
        position_dict = defaultdict(list)

        for mod in self.static_mods:
            positions = self.static_map.get_mod_positions(str(mod))
            if positions:
                for position in positions:
                    position_dict[position].append(mod)

        return dict(position_dict)

    def _generate_variable_mod_scenarios(self) -> None:
        """
        Generate all possible position combinations for variable modifications.

        Populates self.variable_map_scenarios with possible position combinations
        for each variable modification.
        """
        self.variable_map_scenarios = {}

        for mod in self.variable_mods:
            positions = self.variable_map.get_mod_positions(str(mod)) or []

            if not mod.all_filled:
                # Generate all possible subsets of positions
                self.variable_map_scenarios[mod.value] = list(
                    variable_position_placement_generator(positions)
                )
            else:
                # For all_filled mods, only empty list or all positions are valid
                self.variable_map_scenarios[mod.value] = [[], positions]

    def _explore_scenarios(
        self,
        current_mod_idx: int = 0,
        current_scenario: Optional[Dict[int, List[Modification]]] = None,
    ) -> Iterator[Dict[int, List[Modification]]]:
        """
        Recursively explore all possible modification scenarios.

        Parameters
        ----------
        current_mod_idx : int, optional
            Index of the current modification being processed.
        current_scenario : Dict[int, List[Modification]], optional
            Current scenario being built.

        Yields
        ------
        Dict[int, List[Modification]]
            Each possible scenario of variable modifications.
        """
        if current_scenario is None:
            current_scenario = {}

        # Base case: processed all modifications
        if current_mod_idx >= len(self.variable_mods):
            yield current_scenario
            return

        current_mod = self.variable_mods[current_mod_idx]
        position_combinations = self.variable_map_scenarios.get(current_mod.value, [[]])

        for positions in position_combinations:
            # Create a copy of the current scenario
            scenario_copy = deepcopy(current_scenario)

            # Add current modification to positions
            for pos in positions:
                if pos not in scenario_copy:
                    scenario_copy[pos] = []
                scenario_copy[pos].append(current_mod)

            # Recursively continue with next modification
            yield from self._explore_scenarios(current_mod_idx + 1, scenario_copy)


def split_chimeric_proforma(proforma_str: str) -> List[str]:
    """
    Splits a ProForma string representing multiple chimeric assignments
    into individual peptidoform strings.

    The split occurs at '+' characters that are not enclosed within any
    type of brackets ('[]', '()', '{}'). This ensures that '+' signs
    within modification tags or ionic species definitions are ignored.

    Args:
        proforma_str: The input ProForma string, potentially containing
                      multiple assignments separated by '+'.

    Returns:
        A list of strings, where each string is a ProForma representation
        of a single peptidoform assignment. Returns a list containing
        the original string if no valid '+' separator is found.
    """
    parts: List[str] = []
    current_part_start: int = 0
    bracket_level: int = 0

    for i, char in enumerate(proforma_str):
        if char in "[{(":
            bracket_level += 1
        elif char in "]})":
            bracket_level = max(0, bracket_level - 1)
        elif char == "+" and bracket_level == 0:
            parts.append(proforma_str[current_part_start:i].strip())
            current_part_start = i + 1

    parts.append(proforma_str[current_part_start:].strip())
    return [part for part in parts if part]
