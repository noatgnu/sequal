import re
from collections import defaultdict
from typing import Any, Dict, Iterator, List, Optional, Pattern, Set, Tuple

from sequal.base_block import BaseBlock
from sequal.resources import monosaccharides


class Modification(BaseBlock):
    """
    Represents a modification block with various properties for sequence analysis.

    This class extends BaseBlock to model biochemical modifications with additional
    properties like regex patterns, modification types, and fragmentation behavior.

    Parameters
    ----------
    value : str
        Short name of the modification.
    position : int, optional
        Position of the modification. Should be provided when assigned to a block.
    regex_pattern : str, optional
        Regular expression pattern for finding modification sites.
    full_name : str, optional
        Full descriptive name of the modification.
    mod_type : str, optional
        Type of modification: "static" or "variable". Default is "static".
    labile : bool, optional
        Whether the modification is labile (important for mass spectrometry).
    labile_number : int, optional
        Order of fragment in labile fragmentation events.
    mass : float, optional
        Mass delta of the modification in Daltons.
    all_filled : bool, optional
        Whether modification occurs at all expected sites.
    """

    KNOWN_SOURCES = {
        "Unimod",
        "U",
        "PSI-MOD",
        "M",
        "RESID",
        "R",
        "XL-MOD",
        "X",
        "XLMOD",
        "GNO",
        "G",
        "MOD",
        "Obs",
        "Formula",
        "Glycan",
    }

    def __init__(
        self,
        value: str,
        position: Optional[int] = None,
        regex_pattern: Optional[str] = None,
        full_name: Optional[str] = None,
        mod_type: str = "static",
        labile: bool = False,
        labile_number: int = 0,
        mass: float = 0.0,
        all_filled: bool = False,
        crosslink_id: Optional[str] = None,
        is_crosslink_ref: bool = False,
        is_branch_ref: bool = False,
        is_branch: bool = False,
        ambiguity_group: Optional[str] = None,
        is_ambiguity_ref: bool = False,
        in_range: bool = False,
        range_start: Optional[int] = None,
        range_end: Optional[int] = None,
        localization_score: Optional[float] = None,
        mod_value: Optional["ModificationValue"] = None,
        position_constraint: Optional[List[str]] = None,
        limit_per_position: int = 1,
        colocalize_known: bool = False,
        colocalize_unknown: bool = False,
        is_ion_type: bool = False,
    ):
        """
        Initializes a Modification object.

        Args:
            value (str): Short name of the modification.
            position (int, optional): Position of the modification. Should be
                provided when assigned to a block. Defaults to None.
            regex_pattern (str, optional): Regular expression pattern for
                finding modification sites. Defaults to None.
            full_name (str, optional): Full descriptive name of the
                modification. Defaults to None.
            mod_type (str, optional): Type of modification: "static" or
                "variable". Defaults to "static".
            labile (bool, optional): Whether the modification is labile.
                Defaults to False.
            labile_number (int, optional): Order of fragment in labile
                fragmentation events. Defaults to 0.
            mass (float, optional): Mass delta of the modification in Daltons.
                Defaults to 0.0.
            all_filled (bool, optional): Whether modification occurs at all
                expected sites. Defaults to False.
            crosslink_id (str, optional): The crosslink identifier. Defaults to None.
            is_crosslink_ref (bool, optional): Whether this modification is a
                crosslink reference. Defaults to False.
            is_branch_ref (bool, optional): Whether this modification is a
                branch reference. Defaults to False.
            is_branch (bool, optional): Whether this modification is a branch.
                Defaults to False.
            ambiguity_group (str, optional): The ambiguity group of the
                modification. Defaults to None.
            is_ambiguity_ref (bool, optional): Whether this modification is an
                ambiguity reference. Defaults to False.
            in_range (bool, optional): Whether the modification is in a range.
                Defaults to False.
            range_start (int, optional): The start of the range. Defaults to None.
            range_end (int, optional): The end of the range. Defaults to None.
            localization_score (float, optional): The localization score.
                Defaults to None.
            mod_value (ModificationValue, optional): The modification value
                object. Defaults to None.
            position_constraint (List[str], optional): A list of position
                constraints. Defaults to None.
            limit_per_position (int, optional): The limit per position.
                Defaults to 1.
            colocalize_known (bool, optional): Whether to colocalize with
                known modifications. Defaults to False.
            colocalize_unknown (bool, optional): Whether to colocalize with
                unknown modifications. Defaults to False.
            is_ion_type (bool, optional): Whether the modification is an ion
                type. Defaults to False.
        """
        self._source = None
        self._original_value = value
        self._crosslink_id = crosslink_id
        self._is_crosslink_ref = is_crosslink_ref
        self._is_branch_ref = is_branch_ref
        self._is_branch = is_branch
        self._is_ambiguity_ref = is_ambiguity_ref
        self._ambiguity_group = ambiguity_group
        self.in_range = in_range
        self.range_start = range_start
        self.range_end = range_end
        self.localization_score = localization_score
        self._mod_value = mod_value or ModificationValue(value, mass=mass)
        self.position_constraint = position_constraint
        self.limit_per_position = limit_per_position
        self.colocalize_known = colocalize_known
        self.colocalize_unknown = colocalize_unknown
        self.is_ion_type = is_ion_type

        if value.startswith("#") and is_crosslink_ref:
            self._crosslink_id = value[1:]
            value = "#" + self._crosslink_id

        super().__init__(value, position=position, branch=True, mass=mass)

        valid_mod_types = {
            "static",
            "variable",
            "terminal",
            "ambiguous",
            "crosslink",
            "branch",
            "gap",
            "labile",
            "unknown_position",
            "global",
        }

        if (crosslink_id or is_crosslink_ref) and mod_type not in {"crosslink"}:
            mod_type = "crosslink"
        if mod_type not in valid_mod_types:
            raise ValueError(f"mod_type must be one of: {', '.join(valid_mod_types)}")

        self._regex: Optional[Pattern] = (
            re.compile(regex_pattern) if regex_pattern else None
        )

        self._mod_type = mod_type
        self._labile = labile
        self._labile_number = labile_number
        self._full_name = full_name
        self._all_filled = all_filled
        if mod_type == "labile":
            self._labile = True
        if self.in_range:
            self._mod_type = "ambiguous"

    @property
    def value(self):
        """str: The modification value."""
        return self._mod_value.primary_value if self._mod_value else self._value

    @property
    def mass(self):
        """float: The mass of the modification."""
        return self._mod_value.mass if self._mod_value else self._mass

    @property
    def observed_mass(self):
        """float: The observed mass of the modification."""
        return self._mod_value.observed_mass if self._mod_value else 0

    @property
    def ambiguity_group(self):
        """str: The ambiguity group of the modification."""
        return self._mod_value.ambiguity_group if self._mod_value else None

    @property
    def is_ambiguity_ref(self):
        """bool: True if the modification is an ambiguity reference."""
        return (
            self._mod_value.is_ambiguity_ref
            if self._mod_value
            else self._is_ambiguity_ref
        )

    @property
    def synonyms(self):
        """List[str]: A list of synonyms for the modification."""
        return self._mod_value.synonyms

    @staticmethod
    def _validate_formula(formula: str) -> bool:
        """
        Validates a chemical formula according to ProForma specification.

        Args:
            formula (str): The chemical formula to validate.

        Returns:
            bool: True if the formula is valid, False otherwise.
        """
        # Empty formula is invalid
        if not formula.strip():
            return False

        if formula.count("[") != formula.count("]"):
            return False

        formula_no_spaces = formula.replace(" ", "")

        i = 0
        while i < len(formula_no_spaces):
            if formula_no_spaces[i] == "[":
                end_bracket = formula_no_spaces.find("]", i)
                if end_bracket == -1:
                    return False

                isotope_part = formula_no_spaces[i + 1 : end_bracket]
                if not re.match(r"\d+[A-Z][a-z]?(-?\d+)?", isotope_part):
                    return False
                i = end_bracket + 1

                if i < len(formula_no_spaces) and (
                    formula_no_spaces[i] == "-" or formula_no_spaces[i].isdigit()
                ):
                    start = i
                    if formula_no_spaces[i] == "-":
                        i += 1
                    while i < len(formula_no_spaces) and formula_no_spaces[i].isdigit():
                        i += 1
                    if int(formula_no_spaces[start:i]) == 0:
                        return False

            elif formula_no_spaces[i].isupper():
                if (
                    i + 1 < len(formula_no_spaces)
                    and formula_no_spaces[i + 1].islower()
                ):
                    i += 2
                else:
                    i += 1

                if i < len(formula_no_spaces) and (
                    formula_no_spaces[i] == "-" or formula_no_spaces[i].isdigit()
                ):
                    start = i
                    if formula_no_spaces[i] == "-":
                        i += 1
                    while i < len(formula_no_spaces) and formula_no_spaces[i].isdigit():
                        i += 1
                    if int(formula_no_spaces[start:i]) == 0:
                        return False
            else:
                return False

        return True

    @staticmethod
    def _validate_glycan(glycan: str) -> bool:
        """
        Validates a glycan string per ProForma specification.

        Args:
            glycan (str): The glycan string to validate.

        Returns:
            bool: True if the glycan is valid, False otherwise.
        """
        glycan_clean = glycan.replace(" ", "")
        monos = list(monosaccharides)
        monos.sort(key=len, reverse=True)
        mono_pattern = (
            r"^("
            + "|".join(re.escape(m) for m in monos)
            + r")((\(([1-9]\d*)\))|[1-9]\d*)?"
        )
        i = 0
        while i < len(glycan_clean):
            if glycan_clean[i] == "{":
                close_brace = glycan_clean.find("}", i)
                if close_brace == -1:
                    return False

                i = close_brace + 1
                is_at_end = i == len(glycan_clean)

                if i < len(glycan_clean) and glycan_clean[i] == "(":
                    close_paren = glycan_clean.find(")", i)
                    if close_paren == -1:
                        return False
                    count_str = glycan_clean[i + 1 : close_paren]
                    if not re.match(r"^[1-9]\d*$", count_str):
                        return False
                    i = close_paren + 1
                elif i < len(glycan_clean) and glycan_clean[i].isdigit():
                    start = i
                    while i < len(glycan_clean) and glycan_clean[i].isdigit():
                        i += 1
                    count_str = glycan_clean[start:i]
                    if not re.match(r"^[1-9]\d*$", count_str):
                        return False
                elif not is_at_end:
                    return False
            else:
                match = re.match(mono_pattern, glycan_clean[i:])
                if not match:
                    return False

                mono_length = len(match.group(0))
                has_count = match.group(2) is not None and match.group(2) != ""
                i += mono_length

                is_at_end = i == len(glycan_clean)
                if not has_count and not is_at_end:
                    return False

        return i == len(glycan_clean)

    @property
    def mod_value(self) -> Optional["ModificationValue"]:
        """Optional[ModificationValue]: The modification value object."""
        return self._mod_value

    @property
    def info_tags(self) -> List[str]:
        """List[str]: A list of information tags associated with the modification."""
        return self._mod_value.info_tags

    @property
    def crosslink_id(self) -> Optional[str]:
        """Optional[str]: The crosslink identifier."""
        return self.mod_value.crosslink_id if self._mod_value else self._crosslink_id

    @property
    def is_crosslink_ref(self) -> bool:
        """bool: True if this modification is a crosslink reference."""
        return (
            self._mod_value.is_crosslink_ref
            if self._mod_value
            else self._is_crosslink_ref
        )

    @property
    def source(self) -> Optional[str]:
        """Optional[str]: The modification database source."""
        return self._mod_value.source if self._mod_value else self._source

    @property
    def original_value(self) -> str:
        """str: The original value including any source prefix."""
        return self._original_value

    @property
    def regex(self) -> Optional[Pattern]:
        """Optional[Pattern]: The compiled regex pattern for finding modification sites."""
        return self._regex

    @property
    def mod_type(self) -> str:
        """str: The modification type (e.g., 'static', 'variable')."""
        return self._mod_type

    @property
    def labile(self) -> bool:
        """bool: True if the modification is labile."""
        return self._labile

    @property
    def labile_number(self) -> int:
        """int: The labile fragmentation order number."""
        return self._labile_number

    @property
    def full_name(self) -> Optional[str]:
        """Optional[str]: The full descriptive name of the modification."""
        return self._full_name

    @property
    def all_filled(self) -> bool:
        """bool: True if the modification occurs at all expected sites."""
        return self._all_filled

    def find_positions(self, seq: str) -> Iterator[Tuple[int, int]]:
        """
        Find positions of the modification in the given sequence.

        Parameters
        ----------
        seq : str
            The sequence to search for modification sites.

        Yields
        ------
        Tuple[int, int]
            Start and end positions of each match in the sequence.

        Raises
        ------
        ValueError
            If no regex pattern was defined for this modification.
        """
        if not self._regex:
            raise ValueError(
                f"No regex pattern defined for modification '{self.value}'"
            )

        for match in self._regex.finditer(seq):
            groups = match.groups()
            if groups:
                for group_idx in range(len(groups) + 1):
                    yield match.start(group_idx), match.end(group_idx)
            else:
                yield match.start(), match.end()

    def to_dict(self) -> Dict[str, Any]:
        """
        Converts the modification to a dictionary representation.

        Returns:
            Dict[str, Any]: A dictionary containing the modification's attributes.
        """
        result = super().to_dict()
        result.update(
            {
                "source": self._source,
                "original_value": self._original_value,
                "regex_pattern": self._regex.pattern if self._regex else None,
                "full_name": self._full_name,
                "mod_type": self._mod_type,
                "labile": self._labile,
                "labile_number": self._labile_number,
                "all_filled": self._all_filled,
                "crosslink_id": self._crosslink_id,
                "is_crosslink_ref": self._is_crosslink_ref,
            }
        )
        return result

    def __eq__(self, other) -> bool:
        """
        Checks if two modifications are equal.

        Args:
            other (Any): The object to compare with.

        Returns:
            bool: True if the modifications are equal, False otherwise.
        """
        if not super().__eq__(other):
            return False
        if not isinstance(other, Modification):
            return False
        return (
            self._mod_type == other.mod_type
            and self._labile == other.labile
            and self._labile_number == other.labile_number
        )

    def __hash__(self) -> int:
        """
        Generates a hash for the modification.

        Returns:
            int: The hash of the modification.
        """
        base_hash = super().__hash__()
        return hash((base_hash, self._mod_type, self._labile, self._labile_number))

    def __str__(self) -> str:
        """
        Returns a string representation of the modification.

        Returns:
            str: The string representation of the modification.
        """
        if self._is_crosslink_ref and self._crosslink_id:
            return f"#{self._crosslink_id}"
        if self._is_branch_ref:
            return "#BRANCH"

        result = self._mod_value.to_string()
        if self._crosslink_id and not self._is_crosslink_ref:
            result += f"#{self._crosslink_id}"
        if self._is_branch and not self._is_branch_ref:
            result += "#BRANCH"
        if self._labile:
            result += f"{self._labile_number}"

        return result

    def __repr__(self) -> str:
        """
        Returns a detailed string representation for debugging.

        Returns:
            str: The detailed string representation of the modification.
        """
        return f"Modification(value='{self.value}', position={self.position},mod_type='{self._mod_type}', labile={self._labile}, crosslink_id={self._crosslink_id!r}, is_crosslink_ref={self._is_crosslink_ref}, is_branch={self._is_branch}, is_branch_ref={self._is_branch_ref})"

    def has_ambiguity(self) -> bool:
        """
        Checks if the modification has ambiguity.

        Returns:
            bool: True if the modification has ambiguity, False otherwise.
        """
        return any([v.type == PipeValue.AMBIGUITY for v in self.mod_value._pipe_values])

    def has_crosslink(self) -> bool:
        """
        Checks if the modification has a crosslink.

        Returns:
            bool: True if the modification has a crosslink, False otherwise.
        """
        return any([v.type == PipeValue.CROSSLINK for v in self.mod_value._pipe_values])

    def has_branch(self) -> bool:
        """
        Checks if the modification has a branch.

        Returns:
            bool: True if the modification has a branch, False otherwise.
        """
        return any([v.type == PipeValue.BRANCH for v in self.mod_value._pipe_values])

    def to_proforma(self) -> str:
        """
        Convert the modification to ProForma notation string.

        Returns
        -------
        str
            The ProForma string representation of this modification.
        """
        parts = []
        if self.mod_value:
            seen = set()
            for pv in self.mod_value.pipe_values:
                mod_part = ""
                if pv.source:
                    mod_part = f"{pv.source}:"
                    if pv.mass:
                        if pv.mass > 0:
                            mod_part += f"+{pv.mass}"
                            seen.add(f"+{pv.mass}")
                        elif pv.mass < 0:
                            mod_part += f"-{pv.mass}"
                            seen.add(f"{pv.mass}")
                    else:
                        mod_part += f"{pv.value}"
                        if pv.charge is not None and pv.source.upper() == "FORMULA":
                            mod_part += f":z{pv.charge}"
                else:
                    if pv.mass:
                        if pv.mass > 0:
                            mod_part = f"+{pv.mass}"
                        elif pv.mass < 0:
                            mod_part = f"{pv.mass}"
                    elif pv.type == PipeValue.SYNONYM:
                        mod_part = f"{pv.value}"
                    else:
                        if "#" not in pv.value:
                            mod_part = f"{pv.value}"

                if pv.type == PipeValue.CROSSLINK and pv.crosslink_id:
                    mod_part += f"#{pv.crosslink_id}"
                elif pv.type == PipeValue.BRANCH and pv.is_branch:
                    mod_part += f"#BRANCH"
                elif pv.type == PipeValue.AMBIGUITY and pv.ambiguity_group:
                    score_str = (
                        f"({pv.localization_score:.2f})"
                        if pv.localization_score
                        else ""
                    )
                    mod_part += f"#{pv.ambiguity_group}{score_str}"

                if mod_part in seen:
                    continue
                parts.append(mod_part)
                seen.add(mod_part)

            result = "|".join(parts)

            if self.position_constraint:
                result += f"|Position:{','.join(self.position_constraint)}"
            if self.limit_per_position > 1:
                result += f"|Limit:{self.limit_per_position}"
            if self.colocalize_known:
                result += "|CoMKP"
            if self.colocalize_unknown:
                result += "|CoMUP"

            return result
        else:
            if self.mass is not None and self.value.startswith(("+", "-")):
                return str(self.mass)
            return self.value


class ModificationMap:
    """
    Maps modifications to their positions in a sequence for quick lookup.

    This class provides efficient access to modification positions and
    modification objects by name or position.

    Parameters
    ----------
    seq : str
        The sequence to be analyzed.
    mods : List[Modification]
        A list of Modification objects to be mapped.
    ignore_positions : Set[int], optional
        A set of positions to ignore when mapping.
    parse_position : bool, optional
        Whether to parse positions of modifications. Default is True.
    mod_position_dict : Dict[str, List[int]], optional
        Pre-computed dict of modification positions.
    """

    def __init__(
        self,
        seq: str,
        mods: List["Modification"],
        ignore_positions: Optional[Set[int]] = None,
        parse_position: bool = True,
        mod_position_dict: Optional[Dict[str, List[int]]] = None,
    ):
        self.seq = seq
        self.ignore_positions = ignore_positions or set()
        self.mod_dict_by_name: Dict[str, "Modification"] = {}
        self.mod_position_dict: Dict[str, List[int]] = mod_position_dict or {}
        self.position_to_mods: Dict[int, List["Modification"]] = defaultdict(list)

        self._build_mappings(mods, parse_position)

    def _build_mappings(self, mods: List["Modification"], parse_position: bool) -> None:
        """
        Build internal mappings between modifications and positions.

        Parameters
        ----------
        mods : List[Modification]
            List of modifications to map
        parse_position : bool
            Whether to use regex to find positions
        """
        for mod in mods:
            mod_name = str(mod)
            self.mod_dict_by_name[mod_name] = mod

            if parse_position:
                positions = []
                try:
                    for p_start, _ in mod.find_positions(self.seq):
                        if p_start not in self.ignore_positions:
                            positions.append(p_start)
                            self.position_to_mods[p_start].append(mod)
                except ValueError:
                    # No regex pattern defined, skip position parsing
                    pass

                self.mod_position_dict[mod_name] = positions

    def get_mod_positions(self, mod_name: str) -> Optional[List[int]]:
        """
        Get the positions of a modification by its name.

        Parameters
        ----------
        mod_name : str
            The name of the modification.

        Returns
        -------
        List[int] or None
            List of positions where the modification is found, or None if not found.
        """
        return self.mod_position_dict.get(mod_name)

    def get_mod(self, mod_name: str) -> Optional["Modification"]:
        """
        Get the Modification object by its name.

        Parameters
        ----------
        mod_name : str
            The name of the modification.

        Returns
        -------
        Modification or None
            The Modification object, or None if not found.
        """
        return self.mod_dict_by_name.get(mod_name)

    def get_mods_at_position(self, position: int) -> List["Modification"]:
        """
        Get all modifications at a specific position.

        Parameters
        ----------
        position : int
            The position to check for modifications.

        Returns
        -------
        List[Modification]
            List of modifications at the specified position.
        """
        return self.position_to_mods.get(position, [])

    def has_mod_at_position(
        self, position: int, mod_name: Optional[str] = None
    ) -> bool:
        """
        Check if a position has any modification or a specific modification.

        Parameters
        ----------
        position : int
            The position to check
        mod_name : str, optional
            The name of a specific modification to check for

        Returns
        -------
        bool
            True if position has the specified modification(s)
        """
        mods = self.get_mods_at_position(position)
        if not mods:
            return False
        if mod_name is None:
            return True
        return any(str(mod) == mod_name for mod in mods)

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the modification map to a dictionary representation.

        Returns
        -------
        Dict[str, Any]
            Dictionary containing the map's data.
        """
        return {
            "sequence": self.seq,
            "modifications": {
                name: [pos for pos in positions]
                for name, positions in self.mod_position_dict.items()
                if positions
            },
        }


class GlobalModification(Modification):
    """
    Represents a global modification that applies to specified residues.

    Attributes:
        target_residues (Optional[List[str]]): For fixed modifications, the
            target residue types (e.g., ["C", "M"]).
        global_mod_type (str): Type of global modification: "isotope" or "fixed".
    """

    def __init__(
        self,
        value: str,
        target_residues: Optional[List[str]] = None,
        mod_type: str = "isotope",
    ):
        """
        Initializes a global modification.

        Args:
            value (str): The modification value (name, accession, or formula).
            target_residues (Optional[List[str]]): For fixed modifications, the
                target residue types (e.g., ["C", "M"]). Defaults to None.
            mod_type (str): Type of global modification: "isotope" or "fixed".
                Defaults to "isotope".

        Raises:
            ValueError: If mod_type is not 'isotope' or 'fixed'.
        """
        if mod_type not in ["isotope", "fixed"]:
            raise ValueError("Global modification type must be 'isotope' or 'fixed'")
        super().__init__(
            value=value,
            position=None,
            mod_type="global",
        )
        mod_value = ModificationValue(value)
        self._mod_value = mod_value
        self.target_residues = target_residues
        self.global_mod_type = mod_type

    def to_proforma(self) -> str:
        """
        Converts to ProForma notation.

        Returns:
            str: The ProForma notation string.
        """
        if self.global_mod_type == "isotope":
            return f"<{super().to_proforma()}>"
        else:
            mod_value = super().to_proforma()
            if not mod_value.startswith("["):
                mod_str = f"[{mod_value}]"
            else:
                mod_str = mod_value

            target_strings = []
            for target in self.target_residues:
                if isinstance(target, dict):
                    if target.get("type") == "terminal_specific":
                        target_strings.append(
                            f"{target['terminal']}:{target['amino_acid']}"
                        )
                    elif target.get("type") == "terminal":
                        target_strings.append(target["terminal"])
                    elif target.get("type") == "amino_acid":
                        target_strings.append(target["residue"])
                else:
                    target_strings.append(target)

            targets = ",".join(target_strings)
            return f"<{mod_str}@{targets}>"

    def __repr__(self) -> str:
        """
        Returns a detailed string representation for debugging.

        Returns:
            str: The detailed string representation.
        """
        base = f"GlobalModification(value='{self.value}'"
        if self.source:
            base += f", source='{self.source}'"
        if self.target_residues:
            base += f", target_residues={self.target_residues}"
        base += f", mod_type='{self.global_mod_type}')"
        return base


class PipeValue:
    """
    Represents a single pipe-separated value in a modification.
    """

    SYNONYM = "synonym"
    INFO_TAG = "info_tag"
    MASS = "mass"
    OBSERVED_MASS = "observed_mass"
    CROSSLINK = "crosslink"
    BRANCH = "branch"
    AMBIGUITY = "ambiguity"
    GLYCAN = "glycan"
    GAP = "gap"
    FORMULA = "formula"

    def __init__(self, value: str, value_type: str, original_value: str = None):
        """
        Initializes a PipeValue object.

        Args:
            value (str): The value of the pipe-separated component.
            value_type (str): The type of the pipe-separated component.
            original_value (str, optional): The original value of the
                component. Defaults to None.
        """
        self.value = value
        self._type = value_type
        self.crosslink_id = None
        self.is_branch = False
        self.is_branch_ref = False
        self.is_crosslink_ref = False
        self.ambiguity_group = None
        self.is_ambiguity_ref = False
        self.localization_score = None
        self.source = None
        self.original_value = original_value
        self.mass = None
        self.observed_mass = None
        self.is_valid_glycan = False
        self.is_valid_formula = False
        self.charge = None
        self.charge_value = 0
        self._extract_properties()
        self.assigned_types: List[str] = []

    def _extract_properties(self):
        """Extracts special properties from the value based on type."""
        if self.type == self.CROSSLINK and "#" in self.value:
            parts = self.value.split("#", 1)
            if parts[1] == "BRANCH":
                self.is_branch = True
            else:
                self.crosslink_id = parts[1]

        elif self.type == self.AMBIGUITY and "#" in self.value:
            parts = self.value.split("#", 1)
            self.ambiguity_group = parts[1]
            if "(" in self.ambiguity_group and ")" in self.ambiguity_group:
                score_match = re.search(r"\(([\d.]+)\)", self.ambiguity_group)
                if score_match:
                    try:
                        self.localization_score = float(score_match.group(1))
                    except ValueError:
                        pass

    def __str__(self) -> str:
        """
        Returns a string representation of the PipeValue.

        Returns:
            str: The string representation of the PipeValue.
        """
        return self.value

    @property
    def type(self) -> str:
        """str: The type of the pipe-separated component."""
        return self._type

    @type.setter
    def type(self, value: str):
        self._type = value
        if len(self.assigned_types) > 0:
            self.assigned_types[0] = value
        else:
            self.assign_type(value)

    def assign_type(self, value: str):
        """
        Assigns a type to the PipeValue.

        Args:
            value (str): The type to assign.
        """
        if value not in self.assigned_types:
            self.assigned_types.append(value)


class ModificationValue:
    """
    Represents a modification value with unified pipe value handling.
    """

    KNOWN_SOURCES = {
        "Unimod",
        "U",
        "PSI-MOD",
        "M",
        "RESID",
        "R",
        "XL-MOD",
        "X",
        "XLMOD",
        "GNO",
        "G",
        "MOD",
        "Obs",
        "Formula",
        "FORMULA",
        "GLYCAN",
        "Glycan",
        "Info",
        "INFO",
        "OBS",
        "INFO",
        "XL",
    }

    def __init__(self, value: str, mass: Optional[float] = None):
        """
        Initializes a ModificationValue object.

        Args:
            value (str): The modification value string.
            mass (Optional[float]): The mass of the modification. Defaults to None.
        """
        self._primary_value = ""
        self._source = None
        self._mass = mass

        self._pipe_values = []
        self._parse_value(value)

    def _parse_value(self, value: str):
        """
        Parses the modification value string.

        Args:
            value (str): The modification value string.
        """
        if "|" in value:
            components = value.split("|")
            self._process_primary_value(components[0])
            for component in components[1:]:
                self._process_pipe_component(component)
        else:
            self._process_primary_value(value)

    def _process_primary_value(self, value: str):
        """
        Processes the primary value component of the modification value.

        Args:
            value (str): The primary value component.
        """
        if value == "#BRANCH":
            self._primary_value = ""
            pipe_val = PipeValue(value, PipeValue.BRANCH, value)
            pipe_val.is_branch_ref = True
            pipe_val.is_branch = True
            self._pipe_values.append(pipe_val)
            return
        elif value.startswith("#"):
            self._primary_value = ""
            pipe_val = PipeValue(
                value,
                (
                    PipeValue.CROSSLINK
                    if value[1:].startswith("XL")
                    else PipeValue.AMBIGUITY
                ),
                value,
            )
            pipe_val.is_crosslink_ref = pipe_val.type == PipeValue.CROSSLINK
            pipe_val.is_ambiguity_ref = pipe_val.type == PipeValue.AMBIGUITY
            pipe_val.crosslink_id = value[1:] if pipe_val.is_crosslink_ref else None
            pipe_val.ambiguity_group = value[1:] if pipe_val.is_ambiguity_ref else None
            if pipe_val.ambiguity_group:
                if "(" in pipe_val.ambiguity_group and ")" in pipe_val.ambiguity_group:
                    score_match = re.search(r"\(([\d.]+)\)", pipe_val.ambiguity_group)
                    if score_match:
                        try:
                            pipe_val.localization_score = float(score_match.group(1))
                            pipe_val.ambiguity_group = pipe_val.ambiguity_group.replace(
                                score_match.group(0), ""
                            )
                        except ValueError:
                            pass
            self._pipe_values.append(pipe_val)
            return

        # Handle source prefix
        if ":" in value:
            parts = value.split(":", 1)
            if parts[0] in self.KNOWN_SOURCES:
                self._source = parts[0]
                self._primary_value = parts[1]
                is_valid_glycan = False
                is_valid_formula = False
                if self._source.upper() == "FORMULA":
                    is_valid_formula = self._validate_formula(self._primary_value)
                elif self._source.upper() == "GLYCAN":
                    is_valid_glycan = self._validate_glycan(self._primary_value)

                if "#" in self._primary_value:
                    pv_parts = self._primary_value.split("#", 1)
                    self._primary_value = pv_parts[0]
                    if self._source in ["XL", "XLMOD", "XL-MOD", "X"]:
                        pipe_val = PipeValue(
                            f"{self._primary_value}", PipeValue.CROSSLINK
                        )
                        pipe_val.source = self._source
                        pipe_val.crosslink_id = pv_parts[1]
                    elif pv_parts[1] == "BRANCH":
                        pipe_val = PipeValue(pv_parts[0], PipeValue.BRANCH)
                        pipe_val.source = self._source
                        pipe_val.is_branch = True
                    else:
                        pipe_val = PipeValue(
                            f"{self._primary_value}", PipeValue.AMBIGUITY, value
                        )
                        if is_valid_glycan:
                            pipe_val.is_valid_glycan = is_valid_glycan
                            pipe_val.assign_type("glycan")
                        elif is_valid_formula:
                            pipe_val.is_valid_formula = is_valid_formula
                            pipe_val.assign_type("formula")
                        if self._source.upper() == "GNO" or self._source.upper() == "G":
                            pipe_val.is_valid_glycan = True
                            pipe_val.assign_type("glycan")
                        pipe_val.source = self._source
                        pipe_val.ambiguity_group = pv_parts[1]
                        if (
                            "(" in pipe_val.ambiguity_group
                            and ")" in pipe_val.ambiguity_group
                        ):
                            score_match = re.search(
                                r"\(([\d.]+)\)", pipe_val.ambiguity_group
                            )
                            if score_match:
                                try:
                                    pipe_val.localization_score = float(
                                        score_match.group(1)
                                    )
                                    pipe_val.ambiguity_group = (
                                        pipe_val.ambiguity_group.replace(
                                            score_match.group(0), ""
                                        )
                                    )
                                except ValueError:
                                    pass
                        pipe_val.source = self._source
                    self._pipe_values.append(pipe_val)
                else:
                    if self._source.upper() == "INFO":
                        pipe_val = PipeValue(parts[1], PipeValue.INFO_TAG, value)
                        pipe_val.source = self._source
                    elif self._source.upper() == "OBS":
                        pipe_val = PipeValue(parts[1], PipeValue.OBSERVED_MASS, value)
                        pipe_val.source = self._source
                        pipe_val.observed_mass = float(parts[1])
                    elif self._source.upper() == "GLYCAN":
                        pipe_val = PipeValue(parts[1], PipeValue.GLYCAN, value)
                        pipe_val.source = self._source
                        pipe_val.is_valid_glycan = is_valid_glycan
                    elif self._source.upper() == "GNO" or self._source.upper() == "G":
                        pipe_val = PipeValue(parts[1], PipeValue.GAP, value)
                        pipe_val.source = self._source
                        pipe_val.is_valid_glycan = True
                    elif self._source.upper() == "FORMULA":
                        formula_str = parts[1]
                        charge_str = None
                        charge_val = 0

                        if ":z" in formula_str:
                            formula_parts = formula_str.rsplit(":z", 1)
                            formula_str = formula_parts[0]
                            charge_str = formula_parts[1]
                            try:
                                charge_val = int(charge_str)
                            except ValueError:
                                pass

                        pipe_val = PipeValue(formula_str, PipeValue.FORMULA, value)
                        pipe_val.source = self._source
                        pipe_val.is_valid_formula = (
                            is_valid_formula or self._validate_formula(formula_str)
                        )
                        pipe_val.charge = charge_str
                        pipe_val.charge_value = charge_val
                    else:
                        pipe_val = PipeValue(parts[1], PipeValue.SYNONYM, value)
                        pipe_val.source = self._source
                    self._pipe_values.append(pipe_val)

            elif parts[0].upper() == "MASS":
                self._primary_value = value
                try:
                    self._mass = float(parts[1])
                    pipe_val = PipeValue(parts[1], PipeValue.MASS, value)
                    pipe_val.mass = self._mass
                    if "#" in parts[1]:
                        pv_parts = self._primary_value.split("#", 1)
                        self._primary_value = pv_parts[0]
                        pipe_val.value = pv_parts[0]
                        if pv_parts[1] == "BRANCH":
                            pipe_val.is_branch = True
                            pipe_val.type = PipeValue.BRANCH

                        elif pv_parts[1].startswith("XL"):
                            pipe_val.crosslink_id = pv_parts[1]
                            pipe_val.type = PipeValue.CROSSLINK
                        else:
                            pipe_val.ambiguity_group = pv_parts[1]
                            pipe_val.type = PipeValue.AMBIGUITY
                            if (
                                "(" in pipe_val.ambiguity_group
                                and ")" in pipe_val.ambiguity_group
                            ):
                                score_match = re.search(
                                    r"\(([\d.]+)\)", pipe_val.ambiguity_group
                                )
                                if score_match:
                                    try:
                                        pipe_val.localization_score = float(
                                            score_match.group(1)
                                        )
                                        pipe_val.ambiguity_group = (
                                            pipe_val.ambiguity_group.replace(
                                                score_match.group(0), ""
                                            )
                                        )
                                    except ValueError:
                                        pass
                        pipe_val.assign_type(PipeValue.MASS)
                    self._pipe_values.append(pipe_val)

                except ValueError:
                    pass
            else:
                self._primary_value = value
                pipe_val = PipeValue(value, PipeValue.SYNONYM, value)
        else:
            if "#" in value:
                parts = value.split("#", 1)
                self._primary_value = parts[0]
                if parts[1] == "BRANCH":
                    pipe_val = PipeValue(f"{parts[0]}", PipeValue.BRANCH, value)
                    pipe_val.is_branch = True

                elif parts[1].startswith("XL"):
                    pipe_val = PipeValue(f"{parts[0]}", PipeValue.CROSSLINK, value)
                    pipe_val.crosslink_id = parts[1]
                else:
                    pipe_val = PipeValue(f"{parts[0]}", PipeValue.AMBIGUITY, value)
                    pipe_val.ambiguity_group = parts[1]
                    if (
                        "(" in pipe_val.ambiguity_group
                        and ")" in pipe_val.ambiguity_group
                    ):
                        score_match = re.search(
                            r"\(([\d.]+)\)", pipe_val.ambiguity_group
                        )
                        if score_match:
                            try:
                                pipe_val.localization_score = float(
                                    score_match.group(1)
                                )
                                pipe_val.ambiguity_group = (
                                    pipe_val.ambiguity_group.replace(
                                        score_match.group(0), ""
                                    )
                                )
                            except ValueError:
                                pass
                if parts[0].startswith("+") or parts[0].startswith("-"):
                    try:
                        self._mass = float(parts[0])
                        pipe_val.mass = self._mass
                        pipe_val.assign_type(PipeValue.MASS)
                    except ValueError:
                        pass
                else:
                    pipe_val.assign_type(PipeValue.SYNONYM)
                self._pipe_values.append(pipe_val)
            else:
                self._primary_value = value

                if (
                    self._primary_value.startswith("+")
                    or self._primary_value.startswith("-")
                ) and any(c.isdigit() for c in self._primary_value):
                    try:
                        self._mass = float(self._primary_value)
                        pipe_val = PipeValue(self._primary_value, PipeValue.MASS, value)
                        pipe_val.mass = self._mass
                        self._pipe_values.append(pipe_val)
                    except ValueError:
                        pipe_val = PipeValue(value, PipeValue.SYNONYM, value)
                        self._pipe_values.append(pipe_val)
                else:
                    pipe_val = PipeValue(value, PipeValue.SYNONYM, value)
                    self._pipe_values.append(pipe_val)

    def _process_pipe_component(self, component: str):
        """
        Processes a single pipe-separated component.

        Args:
            component (str): The pipe-separated component.
        """
        if component == "#BRANCH":
            pipe_val = PipeValue(component, PipeValue.BRANCH, component)
            pipe_val.is_branch_ref = True
            pipe_val.is_branch = True
            self._pipe_values.append(pipe_val)
            return
        elif component.startswith("#"):
            pipe_val = PipeValue(
                component,
                (
                    PipeValue.CROSSLINK
                    if component[1:].startswith("XL")
                    else PipeValue.AMBIGUITY
                ),
                component,
            )
            pipe_val.is_crosslink_ref = pipe_val.type == PipeValue.CROSSLINK
            pipe_val.is_ambiguity_ref = pipe_val.type == PipeValue.AMBIGUITY
            pipe_val.crosslink_id = component[1:] if pipe_val.is_crosslink_ref else None
            pipe_val.ambiguity_group = (
                component[1:] if pipe_val.is_ambiguity_ref else None
            )
            if pipe_val.ambiguity_group:
                if "(" in pipe_val.ambiguity_group and ")" in pipe_val.ambiguity_group:
                    score_match = re.search(r"\(([\d.]+)\)", pipe_val.ambiguity_group)
                    if score_match:
                        try:
                            pipe_val.localization_score = float(score_match.group(1))
                            pipe_val.ambiguity_group = pipe_val.ambiguity_group.replace(
                                score_match.group(0), ""
                            )
                        except ValueError:
                            pass
            self._pipe_values.append(pipe_val)
            return

        # Handle source prefix
        if ":" in component:
            parts = component.split(":", 1)
            if parts[0] in self.KNOWN_SOURCES:
                source = parts[0]
                value = parts[1]
                is_valid_glycan = False
                is_valid_formula = False
                if source.upper() == "FORMULA":
                    is_valid_formula = self._validate_formula(value)
                elif source.upper() == "GLYCAN":
                    is_valid_glycan = self._validate_glycan(value)

                # Handle crosslinks or ambiguity in value
                if "#" in value:
                    pv_parts = value.split("#", 1)
                    value = pv_parts[0]
                    # Recalculate if value changed due to # split
                    if source.upper() == "FORMULA":
                        is_valid_formula = self._validate_formula(value)
                    elif source.upper() == "GLYCAN":
                        is_valid_glycan = self._validate_glycan(value)

                    if source in ["XL", "XLMOD", "XL-MOD", "X"]:
                        pipe_val = PipeValue(value, PipeValue.CROSSLINK, component)
                        pipe_val.source = source
                        pipe_val.crosslink_id = pv_parts[1]
                    elif pv_parts[1] == "BRANCH":
                        pipe_val = PipeValue(value, PipeValue.BRANCH, component)
                        pipe_val.source = source
                        pipe_val.is_branch = True
                    elif source.upper() == "GLYCAN":
                        pipe_val = PipeValue(value, PipeValue.GLYCAN, component)
                        pipe_val.source = source
                        pipe_val.is_valid_glycan = is_valid_glycan
                    elif source.upper() == "GNO" or source.upper() == "G":
                        pipe_val = PipeValue(value, PipeValue.GLYCAN, component)
                        pipe_val.source = source
                        pipe_val.is_valid_glycan = True
                    elif source.upper() == "FORMULA":
                        pipe_val = PipeValue(value, PipeValue.FORMULA, component)
                        pipe_val.source = source
                        pipe_val.is_valid_formula = is_valid_formula
                    else:
                        pipe_val = PipeValue(value, PipeValue.AMBIGUITY, component)
                        pipe_val.source = source
                        if is_valid_glycan:
                            pipe_val.is_valid_glycan = is_valid_glycan
                            pipe_val.assign_type("glycan")
                        elif is_valid_formula:
                            pipe_val.is_valid_formula = is_valid_formula
                            pipe_val.assign_type("formula")
                        if source.upper() == "GNO" or source.upper() == "G":
                            pipe_val.is_valid_glycan = True
                            pipe_val.assign_type("glycan")
                        pipe_val.ambiguity_group = pv_parts[1]
                        if (
                            "(" in pipe_val.ambiguity_group
                            and ")" in pipe_val.ambiguity_group
                        ):
                            score_match = re.search(
                                r"\(([\d.]+)\)", pipe_val.ambiguity_group
                            )
                            if score_match:
                                try:
                                    pipe_val.localization_score = float(
                                        score_match.group(1)
                                    )
                                    pipe_val.ambiguity_group = (
                                        pipe_val.ambiguity_group.replace(
                                            score_match.group(0), ""
                                        )
                                    )
                                except ValueError:
                                    pass
                    self._pipe_values.append(pipe_val)
                else:
                    if source.upper() == "INFO":
                        pipe_val = PipeValue(value, PipeValue.INFO_TAG, component)
                    elif source.upper() == "OBS":
                        pipe_val = PipeValue(value, PipeValue.OBSERVED_MASS, component)
                        pipe_val.observed_mass = float(value)
                    elif source.upper() == "GLYCAN":
                        pipe_val = PipeValue(parts[1], PipeValue.GLYCAN, value)
                        pipe_val.source = parts[0]
                        pipe_val.is_valid_glycan = is_valid_glycan
                    elif source.upper() == "GNO" or source.upper() == "G":
                        pipe_val = PipeValue(value, PipeValue.GLYCAN, component)
                        pipe_val.source = source
                        pipe_val.is_valid_glycan = True
                    elif source.upper() == "FORMULA":
                        pipe_val = PipeValue(parts[1], PipeValue.FORMULA, value)
                        pipe_val.source = parts[0]
                        pipe_val.is_valid_formula = is_valid_formula

                    else:
                        pipe_val = PipeValue(value, PipeValue.SYNONYM, component)

                    pipe_val.source = source
                    self._pipe_values.append(pipe_val)

            elif parts[0].upper() == "MASS":
                try:
                    mass = float(parts[1])
                    pipe_val = PipeValue(parts[1], PipeValue.MASS, component)
                    pipe_val.mass = mass

                    if "#" in parts[1]:
                        pv_parts = parts[1].split("#", 1)
                        pipe_val.value = pv_parts[0]
                        if pv_parts[1] == "BRANCH":
                            pipe_val.is_branch = True
                            pipe_val.type = PipeValue.BRANCH
                        elif pv_parts[1].startswith("XL"):
                            pipe_val.crosslink_id = pv_parts[1]
                            pipe_val.type = PipeValue.CROSSLINK
                        else:
                            pipe_val.ambiguity_group = pv_parts[1]
                            pipe_val.type = PipeValue.AMBIGUITY
                            if (
                                "(" in pipe_val.ambiguity_group
                                and ")" in pipe_val.ambiguity_group
                            ):
                                score_match = re.search(
                                    r"\(([\d.]+)\)", pipe_val.ambiguity_group
                                )
                                if score_match:
                                    try:
                                        pipe_val.localization_score = float(
                                            score_match.group(1)
                                        )
                                        pipe_val.ambiguity_group = (
                                            pipe_val.ambiguity_group.replace(
                                                score_match.group(0), ""
                                            )
                                        )
                                    except ValueError:
                                        pass
                        pipe_val.assign_type(PipeValue.MASS)
                    self._pipe_values.append(pipe_val)
                except ValueError:
                    self._pipe_values.append(
                        PipeValue(component, PipeValue.SYNONYM, component)
                    )
            else:
                self._pipe_values.append(
                    PipeValue(component, PipeValue.SYNONYM, component)
                )
        else:
            # Handle crosslink ID or ambiguity for values without source prefix
            if "#" in component:
                parts = component.split("#", 1)
                value = parts[0]

                if parts[1] == "BRANCH":
                    pipe_val = PipeValue(value, PipeValue.BRANCH, component)
                    pipe_val.is_branch = True
                elif parts[1].startswith("XL"):
                    pipe_val = PipeValue(value, PipeValue.CROSSLINK, component)
                    pipe_val.crosslink_id = parts[1]
                else:
                    pipe_val = PipeValue(value, PipeValue.AMBIGUITY, component)
                    pipe_val.ambiguity_group = parts[1]
                    if (
                        "(" in pipe_val.ambiguity_group
                        and ")" in pipe_val.ambiguity_group
                    ):
                        score_match = re.search(
                            r"\(([\d.]+)\)", pipe_val.ambiguity_group
                        )
                        if score_match:
                            try:
                                pipe_val.localization_score = float(
                                    score_match.group(1)
                                )
                                pipe_val.ambiguity_group = (
                                    pipe_val.ambiguity_group.replace(
                                        score_match.group(0), ""
                                    )
                                )
                            except ValueError:
                                pass

                if (value.startswith("+") or value.startswith("-")) and any(
                    c.isdigit() for c in value
                ):
                    try:
                        pipe_val.mass = float(value)
                        pipe_val.assign_type(PipeValue.MASS)
                    except ValueError:
                        pass
                else:
                    pipe_val.assign_type(PipeValue.SYNONYM)
                self._pipe_values.append(pipe_val)
            else:
                if (component.startswith("+") or component.startswith("-")) and any(
                    c.isdigit() for c in component
                ):
                    try:
                        mass = float(component)
                        pipe_val = PipeValue(component, PipeValue.MASS, component)
                        pipe_val.mass = mass
                        self._pipe_values.append(pipe_val)
                    except ValueError:
                        self._pipe_values.append(
                            PipeValue(component, PipeValue.SYNONYM, component)
                        )
                else:
                    self._pipe_values.append(
                        PipeValue(component, PipeValue.SYNONYM, component)
                    )

    @property
    def source(self) -> Optional[str]:
        """Optional[str]: The source of the modification."""
        return self._source

    @property
    def primary_value(self) -> str:
        """str: The primary value of the modification."""
        return self._primary_value

    @property
    def mass(self) -> Optional[float]:
        """Optional[float]: The mass of the modification."""
        return self._mass

    @property
    def synonyms(self) -> List[str]:
        """List[str]: A list of synonyms for the modification."""
        return [pv.value for pv in self._pipe_values if pv.type == PipeValue.SYNONYM]

    @property
    def observed_mass(self) -> Optional[str]:
        """Optional[str]: The observed mass of the modification."""
        for pv in self._pipe_values:
            if pv.type == PipeValue.OBSERVED_MASS:
                parts = pv.value.split(":", 1)
                if len(parts) > 1:
                    return parts[1]
        return None

    @property
    def pipe_values(self) -> List[PipeValue]:
        """List[PipeValue]: A list of PipeValue objects."""
        return self._pipe_values

    @property
    def info_tags(self) -> List[str]:
        """List[str]: A list of information tags."""
        return [pv.value for pv in self._pipe_values if pv.type == PipeValue.INFO_TAG]

    @property
    def crosslink_id(self) -> Optional[str]:
        """Optional[str]: The crosslink identifier."""
        for pv in self._pipe_values:
            if pv.type == PipeValue.CROSSLINK and pv.crosslink_id:
                return pv.crosslink_id
        return None

    @property
    def is_branch(self) -> bool:
        """bool: True if the modification is a branch."""
        return any(pv.is_branch for pv in self._pipe_values)

    @property
    def is_branch_ref(self) -> bool:
        """bool: True if the modification is a branch reference."""
        return any(pv.is_branch_ref for pv in self._pipe_values)

    @property
    def is_crosslink_ref(self) -> bool:
        """bool: True if the modification is a crosslink reference."""
        return any(pv.is_crosslink_ref for pv in self._pipe_values)

    @property
    def ambiguity_group(self) -> Optional[str]:
        """Optional[str]: The ambiguity group of the modification."""
        for pv in self._pipe_values:
            if pv.type == PipeValue.AMBIGUITY and pv.ambiguity_group:
                return pv.ambiguity_group
        return None

    @property
    def is_ambiguity_ref(self) -> bool:
        """bool: True if the modification is an ambiguity reference."""
        return any(pv.is_ambiguity_ref for pv in self._pipe_values)

    def to_string(self) -> str:
        """
        Converts to string representation for ProForma output.

        Returns:
            str: The ProForma string representation.
        """
        parts = []

        # Start with primary value and source
        if self._source:
            parts.append(f"{self._source}:{self._primary_value}")
        else:
            parts.append(self._primary_value)

        # Add all pipe values without duplicates
        seen = set()
        for pv in self._pipe_values:
            if pv.value and pv.value not in seen:
                seen.add(pv.value)
                parts.append(pv.value)

        return "|".join(parts)

    @staticmethod
    def _validate_formula(formula: str) -> bool:
        """
        Validates a chemical formula according to ProForma specification.

        Args:
            formula (str): The chemical formula to validate.

        Returns:
            bool: True if the formula is valid, False otherwise.
        """
        # Empty formula is invalid
        if not formula.strip():
            return False

        # Check for balanced brackets
        if formula.count("[") != formula.count("]"):
            return False

        # Remove spaces for processing (allowed by spec)
        formula_no_spaces = formula.replace(" ", "")

        # Process through the formula
        i = 0
        while i < len(formula_no_spaces):
            # Handle isotopes [13C2]
            if formula_no_spaces[i] == "[":
                end_bracket = formula_no_spaces.find("]", i)
                if end_bracket == -1:
                    return False
                # Extract isotope content
                isotope_part = formula_no_spaces[i + 1 : end_bracket]
                # Must start with digits followed by element
                if not re.match(r"\d+[A-Z][a-z]?(-?\d+)?", isotope_part):
                    return False

                i = end_bracket + 1

                # Check for cardinality after bracket
                if i < len(formula_no_spaces) and (
                    formula_no_spaces[i] == "-" or formula_no_spaces[i].isdigit()
                ):
                    start = i
                    if formula_no_spaces[i] == "-":
                        i += 1
                    while i < len(formula_no_spaces) and formula_no_spaces[i].isdigit():
                        i += 1
                    if int(formula_no_spaces[start:i]) == 0:
                        return False

            # Handle regular elements (C12, H, Na+)
            elif formula_no_spaces[i].isupper():
                # Element symbol (1-2 chars)
                if (
                    i + 1 < len(formula_no_spaces)
                    and formula_no_spaces[i + 1].islower()
                ):
                    i += 2
                else:
                    i += 1

                # Check for cardinality
                if i < len(formula_no_spaces) and (
                    formula_no_spaces[i] == "-" or formula_no_spaces[i].isdigit()
                ):
                    start = i
                    if formula_no_spaces[i] == "-":
                        i += 1
                    while i < len(formula_no_spaces) and formula_no_spaces[i].isdigit():
                        i += 1
                    if int(formula_no_spaces[start:i]) == 0:
                        return False
            else:
                # Unexpected character
                return False

        return True

    @staticmethod
    def _validate_glycan(glycan: str) -> bool:
        """
        Validates a glycan string per ProForma specification.

        Args:
            glycan (str): The glycan string to validate.

        Returns:
            bool: True if the glycan is valid, False otherwise.
        """
        glycan_clean = glycan.replace(" ", "")

        monos = list(monosaccharides)
        monos.sort(key=len, reverse=True)
        mono_pattern = (
            r"^("
            + "|".join(re.escape(m) for m in monos)
            + r")((\(([1-9]\d*)\))|[1-9]\d*)?"
        )

        i = 0
        while i < len(glycan_clean):
            if glycan_clean[i] == "{":
                close_brace = glycan_clean.find("}", i)
                if close_brace == -1:
                    return False

                i = close_brace + 1
                is_at_end = i == len(glycan_clean)

                if i < len(glycan_clean) and glycan_clean[i] == "(":
                    close_paren = glycan_clean.find(")", i)
                    if close_paren == -1:
                        return False
                    count_str = glycan_clean[i + 1 : close_paren]
                    if not re.match(r"^[1-9]\d*$", count_str):
                        return False
                    i = close_paren + 1
                elif i < len(glycan_clean) and glycan_clean[i].isdigit():
                    start = i
                    while i < len(glycan_clean) and glycan_clean[i].isdigit():
                        i += 1
                    count_str = glycan_clean[start:i]
                    if not re.match(r"^[1-9]\d*$", count_str):
                        return False
                elif not is_at_end:
                    return False
            else:
                match = re.match(mono_pattern, glycan_clean[i:])
                if not match:
                    return False

                mono_length = len(match.group(0))
                has_count = match.group(2) is not None and match.group(2) != ""
                i += mono_length

                is_at_end = i == len(glycan_clean)
                if not has_count and not is_at_end:
                    return False

        return i == len(glycan_clean)

    def __getitem__(self, index):
        """
        Accesses pipe values by index or slice.

        Args:
            index (int or slice): Index or slice to access pipe values.

        Returns:
            PipeValue or list: Single PipeValue for integer index, list of
                PipeValues for slice.

        Raises:
            IndexError: If index is out of range.
        """
        return self._pipe_values[index]

    def __len__(self):
        """
        Gets the number of pipe values.

        Returns:
            int: The number of pipe values.
        """
        return len(self._pipe_values)

    def __iter__(self):
        """
        Allows iteration through pipe values.

        Returns:
            iterator: An iterator over pipe values.
        """
        return iter(self._pipe_values)
