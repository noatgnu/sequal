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
    ):
        self._source = None
        self._original_value = value
        self._crosslink_id = crosslink_id
        self._is_crosslink_ref = is_crosslink_ref
        self._is_branch_ref = is_branch_ref
        self._is_branch = is_branch

        if ":" in value:
            parts = value.split(":", 1)
            if parts[0] in self.KNOWN_SOURCES:
                self._source = parts[0]
                value = parts[1]
                if "#" in value:
                    value_parts = value.split("#", 1)
                    value = value_parts[0]
                    self._crosslink_id = value_parts[1]
                if self._source == "Formula":
                    if not self._validate_formula(value):
                        raise ValueError(f"Invalid formula: {value}")
                elif self._source == "Glycan":
                    if not self._validate_glycan(value):
                        raise ValueError(f"Invalid glycan: {value}")

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

    @staticmethod
    def _validate_formula(formula: str) -> bool:
        """
        Validate a chemical formula according to the specified rules.

        Validates:
        1. Element symbols followed by optional numbers (C12, H20, O)
        2. Isotopes in brackets ([13C2])
        3. Spaces between elements
        4. Negative cardinalities (C-2)
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
        """Validate a glycan string per ProForma specification."""
        # List of supported monosaccharides

        # Remove spaces for processing
        glycan_clean = glycan.replace(" ", "")

        # Build pattern to match monosaccharide with optional number
        monos = list(monosaccharides)
        monos.sort(key=len, reverse=True)
        mono_pattern = r"^(" + "|".join(re.escape(m) for m in monos) + r")(\d+)?"
        # Check if entire string matches consecutive monosaccharide patterns
        i = 0
        while i < len(glycan_clean):
            match = re.match(mono_pattern, glycan_clean[i:])
            if not match:
                return False
            i += len(match.group(0))

        return i == len(glycan_clean)  # Ensure we consumed the entire string

    @property
    def crosslink_id(self) -> Optional[str]:
        """Get the crosslink identifier."""
        return self._crosslink_id

    @property
    def is_crosslink_ref(self) -> bool:
        """Check if this modification is a crosslink reference."""
        return self._is_crosslink_ref

    @property
    def source(self) -> Optional[str]:
        """Get the modification database source."""
        return self._source

    @property
    def original_value(self) -> str:
        """Get the original value including any source prefix."""
        return self._original_value

    @property
    def regex(self) -> Optional[Pattern]:
        """Get the compiled regex pattern for finding modification sites."""
        return self._regex

    @property
    def mod_type(self) -> str:
        """Get the modification type (static or variable)."""
        return self._mod_type

    @property
    def labile(self) -> bool:
        """Check if the modification is labile."""
        return self._labile

    @property
    def labile_number(self) -> int:
        """Get the labile fragmentation order number."""
        return self._labile_number

    @property
    def full_name(self) -> Optional[str]:
        """Get the full descriptive name of the modification."""
        return self._full_name

    @property
    def all_filled(self) -> bool:
        """Check if the modification occurs at all expected sites."""
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
        """Convert the modification to a dictionary representation."""
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
        """Check if two modifications are equal."""
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
        """Generate a hash for the modification."""
        base_hash = super().__hash__()
        return hash((base_hash, self._mod_type, self._labile, self._labile_number))

    def __str__(self) -> str:
        """Return a string representation of the modification."""
        if self._is_crosslink_ref and self._crosslink_id:
            return f"#{self._crosslink_id}"
        if self._is_branch_ref:
            return "#BRANCH"
        result = ""
        if self._source:
            result = f"{self._source}:{self.value}"
        else:
            result = self.value

        if self._crosslink_id and not self._is_crosslink_ref:
            result += f"#{self._crosslink_id}"
        if self._is_branch and not self._is_branch_ref:
            result += "#BRANCH"
        if self._labile:
            result += f"{self._labile_number}"

        return result

    def __repr__(self) -> str:
        """Return a detailed string representation for debugging."""
        return (
            f"Modification(value='{self.value}', position={self.position}, "
            f"mod_type='{self._mod_type}', labile={self._labile}, "
            f"crosslink_id={self._crosslink_id!r}, is_crosslink_ref={self._is_crosslink_ref}, "
            f"is_branch={self._is_branch}, is_branch_ref={self._is_branch_ref})"
        )


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

        # Maps mod name to Modification object
        self.mod_dict_by_name: Dict[str, "Modification"] = {}

        # Maps mod name to list of positions
        self.mod_position_dict: Dict[str, List[int]] = mod_position_dict or {}

        # Maps position to list of modifications at that position
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
