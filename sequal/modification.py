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
    ):
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
        """Get the modification value."""
        return self._mod_value.primary_value if self._mod_value else self._value

    @property
    def mass(self):
        """Get the mass of the modification."""
        return self._mod_value.mass if self._mod_value else self._mass

    @property
    def observed_mass(self):
        """Get the observed mass of the modification."""
        return self._mod_value.observed_mass if self._mod_value else 0

    @property
    def ambiguity_group(self):
        """Get the ambiguity group of the modification."""
        return self._mod_value.ambiguity_group if self._mod_value else None

    @property
    def is_ambiguity_ref(self):
        """Check if the modification is an ambiguity reference."""
        return (
            self._mod_value.is_ambiguity_ref
            if self._mod_value
            else self._is_ambiguity_ref
        )

    @property
    def synonyms(self):
        """Get the synonyms of the modification."""
        return self._mod_value.synonyms

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
    def mod_value(self) -> Optional["ModificationValue"]:
        """Get the modification value object."""
        return self._mod_value

    @property
    def info_tags(self) -> List[str]:
        """Get the list of information tags associated with the modification."""

        return self._mod_value.info_tags

    @property
    def crosslink_id(self) -> Optional[str]:
        """Get the crosslink identifier."""
        return self.mod_value.crosslink_id if self._mod_value else self._crosslink_id

    @property
    def is_crosslink_ref(self) -> bool:
        """Check if this modification is a crosslink reference."""
        return (
            self._mod_value.is_crosslink_ref
            if self._mod_value
            else self._is_crosslink_ref
        )

    @property
    def source(self) -> Optional[str]:
        """Get the modification database source."""
        return self._mod_value.source if self._mod_value else self._source

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

        result = self._mod_value.to_string()

        # Add special notations after the main value parts
        if self._crosslink_id and not self._is_crosslink_ref:
            result += f"#{self._crosslink_id}"
        if self._is_branch and not self._is_branch_ref:
            result += "#BRANCH"
        if self._labile:
            result += f"{self._labile_number}"

        return result

    def __repr__(self) -> str:
        """Return a detailed string representation for debugging."""
        return f"Modification(value='{self.value}', position={self.position},mod_type='{self._mod_type}', labile={self._labile}, crosslink_id={self._crosslink_id!r}, is_crosslink_ref={self._is_crosslink_ref}, is_branch={self._is_branch}, is_branch_ref={self._is_branch_ref})"

    def has_ambiguity(self) -> bool:
        """Check if the modification has ambiguity."""
        return any([v.type == PipeValue.AMBIGUITY for v in self.mod_value._pipe_values])

    def has_crosslink(self) -> bool:
        """Check if the modification has crosslink."""
        return any([v.type == PipeValue.CROSSLINK for v in self.mod_value._pipe_values])

    def has_branch(self) -> bool:
        """Check if the modification has branch."""
        return any([v.type == PipeValue.BRANCH for v in self.mod_value._pipe_values])

    def to_proforma(self) -> str:
        """
        Convert the modification to ProForma notation string.

        Returns
        -------
        str
            The ProForma string representation of this modification.
        """

        # Start with the base value
        parts = []

        # Get primary value with source prefix if present
        if self.mod_value:
            seen = set()
            for pv in self.mod_value._pipe_values:
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

            return "|".join(parts)
        # If mod_value isn't available, use basic properties
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


class GlobalModification(Modification):
    """Represents a global modification that applies to specified residues."""

    def __init__(
        self,
        value: str,
        target_residues: Optional[List[str]] = None,
        mod_type: str = "isotope",
    ):
        """Initialize a global modification.

        Parameters
        ----------
        value : str
            The modification value (name, accession, or formula)
        target_residues : List[str], optional
            For fixed modifications, the target residue types (e.g., ["C", "M"])
        mod_type : str
            Type of global modification: "isotope" or "fixed"
        """
        if mod_type not in ["isotope", "fixed"]:
            raise ValueError("Global modification type must be 'isotope' or 'fixed'")

        # Initialize base Modification class
        # This will handle source parsing (MOD:, Unimod:, etc.)
        super().__init__(
            value=value,
            position=None,
            mod_type="global",  # Use a specific mod_type to identify global mods
        )
        mod_value = ModificationValue(value)
        self._mod_value = mod_value

        # Global mod specific attributes
        self.target_residues = target_residues
        self.global_mod_type = mod_type

    def to_proforma(self) -> str:
        """Convert to ProForma notation."""
        if self.global_mod_type == "isotope":
            return f"<{super().to_proforma()}>"
        else:
            # For fixed modifications with target residues
            mod_value = super().to_proforma()
            # Format with brackets if not already present
            if not mod_value.startswith("["):
                mod_str = f"[{mod_value}]"
            else:
                mod_str = mod_value
            # Join target residues with commas
            targets = ",".join(self.target_residues)
            return f"<{mod_str}@{targets}>"

    def __repr__(self) -> str:
        """Return a detailed string representation for debugging."""
        base = f"GlobalModification(value='{self.value}'"
        if self.source:
            base += f", source='{self.source}'"
        if self.target_residues:
            base += f", target_residues={self.target_residues}"
        base += f", mod_type='{self.global_mod_type}')"
        return base


class PipeValue:
    """Represents a single pipe-separated value in a modification."""

    # Value types
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
        self.value = value
        self._type = value_type
        # Additional properties for special types
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
        # Extract special properties based on type
        self._extract_properties()
        self.assigned_types: List[str] = []

    def _extract_properties(self):
        """Extract special properties from the value based on type."""
        if self.type == self.CROSSLINK and "#" in self.value:
            parts = self.value.split("#", 1)
            if parts[1] == "BRANCH":
                self.is_branch = True
            else:
                self.crosslink_id = parts[1]

        elif self.type == self.AMBIGUITY and "#" in self.value:
            parts = self.value.split("#", 1)
            self.ambiguity_group = parts[1]
            # Check for localization score
            if "(" in self.ambiguity_group and ")" in self.ambiguity_group:
                score_match = re.search(r"\(([\d.]+)\)", self.ambiguity_group)
                if score_match:
                    try:
                        self.localization_score = float(score_match.group(1))
                    except ValueError:
                        pass

    def __str__(self) -> str:
        return self.value

    @property
    def type(self) -> str:
        return self._type

    @type.setter
    def type(self, value: str):
        self._type = value
        if len(self.assigned_types) > 0:
            self.assigned_types[0] = value
        else:
            self.assign_type(value)

    def assign_type(self, value: str):
        """Assign a type to the PipeValue."""
        if value not in self.assigned_types:
            self.assigned_types.append(value)


class ModificationValue:
    """Represents a modification value with unified pipe value handling."""

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
        # Core properties
        self._primary_value = ""
        self._source = None
        self._mass = mass

        # Use a single collection for all pipe values
        self._pipe_values = []

        # Parse the input value
        self._parse_value(value)

    def _parse_value(self, value: str):
        """Parse modification value with unified pipe value handling."""
        if "|" in value:
            components = value.split("|")
            # Process primary component
            self._process_primary_value(components[0])

            # Process additional pipe components
            for component in components[1:]:
                self._process_pipe_component(component)
        else:
            self._process_primary_value(value)

    def _process_primary_value(self, value: str):
        """Process the primary value component."""
        # Handle special cases first
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
                elif self._source.upper() == "GlYCAN":
                    is_valid_glycan = self._validate_glycan(self._primary_value)
                # Handle crosslinks or ambiguity in primary value
                if "#" in self._primary_value:
                    pv_parts = self._primary_value.split("#", 1)
                    self._primary_value = pv_parts[0]

                    # Create appropriate pipe value
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
                        # parse localization score if present
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
                        pipe_val = PipeValue(parts[1], PipeValue.FORMULA, value)
                        pipe_val.source = self._source
                        pipe_val.is_valid_formula = is_valid_formula
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
            # Handle crosslink ID or ambiguity for values without source prefix
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
        """Process a single pipe-separated component."""
        # Handle special cases first
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

                # Handle crosslinks or ambiguity in value
                if "#" in value:
                    pv_parts = value.split("#", 1)
                    value = pv_parts[0]
                    is_valid_glycan = False
                    is_valid_formula = False
                    if source.upper() == "FORMULA":
                        is_valid_formula = self._validate_formula(value)
                    elif source.upper() == "GlYCAN":
                        is_valid_glycan = self._validate_glycan(value)
                    # Create appropriate pipe value
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
                    elif source.upper() == "GNO" or self._source.upper() == "G":
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
                        if source.upper() == "GNO" or self._source.upper() == "G":
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
                    elif self._source.upper() == "GLYCAN":
                        pipe_val = PipeValue(parts[1], PipeValue.GLYCAN, value)
                        pipe_val.source = parts[0]
                        pipe_val.is_valid_glycan = is_valid_glycan
                    elif source.upper() == "GNO" or self._source.upper() == "G":
                        pipe_val = PipeValue(value, PipeValue.GLYCAN, component)
                        pipe_val.source = source
                        pipe_val.is_valid_glycan = True
                    elif self._source.upper() == "FORMULA":
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

                    # Handle special notation in mass
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

                # Check if base value is a mass
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
                # Handle mass shifts
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
        return self._source

    @property
    def primary_value(self) -> str:
        return self._primary_value

    @property
    def mass(self) -> Optional[float]:
        return self._mass

    @property
    def synonyms(self) -> List[str]:
        return [pv.value for pv in self._pipe_values if pv.type == PipeValue.SYNONYM]

    @property
    def observed_mass(self) -> Optional[str]:
        for pv in self._pipe_values:
            if pv.type == PipeValue.OBSERVED_MASS:
                parts = pv.value.split(":", 1)
                if len(parts) > 1:
                    return parts[1]
        return None

    @property
    def pipe_values(self) -> List[PipeValue]:
        return self._pipe_values

    @property
    def info_tags(self) -> List[str]:
        return [pv.value for pv in self._pipe_values if pv.type == PipeValue.INFO_TAG]

    @property
    def crosslink_id(self) -> Optional[str]:
        for pv in self._pipe_values:
            if pv.type == PipeValue.CROSSLINK and pv.crosslink_id:
                return pv.crosslink_id
        return None

    @property
    def is_branch(self) -> bool:
        return any(pv.is_branch for pv in self._pipe_values)

    @property
    def is_branch_ref(self) -> bool:
        return any(pv.is_branch_ref for pv in self._pipe_values)

    @property
    def is_crosslink_ref(self) -> bool:
        return any(pv.is_crosslink_ref for pv in self._pipe_values)

    @property
    def ambiguity_group(self) -> Optional[str]:
        for pv in self._pipe_values:
            if pv.type == PipeValue.AMBIGUITY and pv.ambiguity_group:
                return pv.ambiguity_group
        return None

    @property
    def is_ambiguity_ref(self) -> bool:
        return any(pv.is_ambiguity_ref for pv in self._pipe_values)

    def to_string(self) -> str:
        """Convert to string representation for ProForma output."""
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

        return i == len(glycan_clean)

    def __getitem__(self, index):
        """
        Access pipe values by index or slice.

        Parameters
        ----------
        index : int or slice
            Index or slice to access pipe values

        Returns
        -------
        PipeValue or list
            Single PipeValue for integer index, list of PipeValues for slice

        Raises
        ------
        IndexError
            If index is out of range
        """
        return self._pipe_values[index]

    def __len__(self):
        """
        Get the number of pipe values.

        Returns
        -------
        int
            Number of pipe values
        """
        return len(self._pipe_values)

    def __iter__(self):
        """
        Allow iteration through pipe values.

        Returns
        -------
        iterator
            Iterator over pipe values
        """
        return iter(self._pipe_values)
