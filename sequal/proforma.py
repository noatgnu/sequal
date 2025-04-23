import re
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from sequal.modification import GlobalModification, Modification, ModificationValue


class ProFormaParser:
    """Parser for the ProForma peptide notation format."""

    # Regex patterns for ProForma notation components
    MASS_SHIFT_PATTERN = re.compile(r"^[+-]\d+(\.\d+)?$")
    TERMINAL_PATTERN = re.compile(r"^\[([^\]]+)\]-(.+)-\[([^\]]+)\]$")
    N_TERMINAL_PATTERN = re.compile(r"^\[([^\]]+)\]-(.+)$")
    C_TERMINAL_PATTERN = re.compile(r"^(.+)-\[([^\]]+)\]$")
    CROSSLINK_PATTERN = re.compile(r"^([^#]+)#(XL[A-Za-z0-9]+)$")
    CROSSLINK_REF_PATTERN = re.compile(r"^#(XL[A-Za-z0-9]+)$")
    BRANCH_PATTERN = re.compile(r"^([^#]+)#BRANCH$")
    BRANCH_REF_PATTERN = re.compile(r"^#BRANCH$")
    UNKNOWN_POSITION_PATTERN = re.compile(r"(\[([^\]]+)\])(\^(\d+))?(\?)")

    @staticmethod
    def parse(proforma_str: str) -> Tuple[str, Dict[int, List[Modification]]]:
        """
        Parse a ProForma string into a base sequence and modifications.

        Parameters
        ----------
        proforma_str : str
            ProForma formatted peptide string

        Returns
        -------
        Tuple[str, Dict[int, List[Modification]]]
            Base sequence and modifications dictionary
        """
        base_sequence = ""
        modifications = defaultdict(list)
        global_mods = []
        sequence_ambiguities = []

        while proforma_str.startswith("<"):
            end_bracket = proforma_str.find(">")
            if end_bracket == -1:
                raise ValueError("Unclosed global modification angle bracket")

            global_mod_str = proforma_str[1:end_bracket]
            proforma_str = proforma_str[end_bracket + 1 :]  # Remove processed part

            if "@" in global_mod_str:
                # Fixed protein modification
                mod_part, targets = global_mod_str.split("@")
                if mod_part.startswith("[") and mod_part.endswith("]"):
                    mod_value = mod_part[1:-1]  # Remove brackets
                else:
                    mod_value = mod_part

                target_residues = targets.split(",")
                global_mods.append(
                    GlobalModification(mod_value, target_residues, "fixed")
                )
            else:
                # Isotope labeling
                global_mods.append(GlobalModification(global_mod_str, None, "isotope"))

        if "?" in proforma_str:
            i = 0
            unknown_pos_mods = []
            while i < len(proforma_str):
                # If not a bracket, check if we've collected unknown position mods
                if proforma_str[i] != "[":
                    if (
                        unknown_pos_mods
                        and i < len(proforma_str)
                        and proforma_str[i] == "?"
                    ):
                        print(unknown_pos_mods)
                        # Add all collected mods to position -4
                        for mod_str in unknown_pos_mods:
                            mod = ProFormaParser._create_modification(
                                mod_str, is_unknown_position=True
                            )
                            modifications[-4].append(mod)
                        i += 1  # Move past the question mark
                    unknown_pos_mods = []  # Reset collection
                    break  # Done with unknown position section

                # Find matching closing bracket
                bracket_count = 1
                j = i + 1
                while j < len(proforma_str) and bracket_count > 0:
                    if proforma_str[j] == "[":
                        bracket_count += 1
                    elif proforma_str[j] == "]":
                        bracket_count -= 1
                    j += 1

                if bracket_count > 0:
                    raise ValueError(f"Unclosed bracket at position {i}")

                mod_str = proforma_str[i + 1 : j - 1]  # Extract modification string

                # Check if followed by caret notation
                count = 1
                if j < len(proforma_str) and proforma_str[j] == "^":
                    j += 1
                    num_start = j
                    while j < len(proforma_str) and proforma_str[j].isdigit():
                        j += 1
                    if j > num_start:
                        count = int(proforma_str[num_start:j])

                # Add to our collection of mods that might be unknown position
                for _ in range(count):
                    unknown_pos_mods.append(mod_str)
                i = j  # Move to next character (might be ? or another [)
            proforma_str = proforma_str[i:]

        i = 0
        while i < len(proforma_str) and proforma_str[i] == "{":
            j = proforma_str.find("}", i)
            if j == -1:
                raise ValueError(f"Unclosed curly brace at position {i}")

            mod_str = proforma_str[i + 1 : j]
            if not mod_str.startswith("Glycan:"):
                raise ValueError(
                    f"Labile modification must start with 'Glycan:', found: {mod_str}"
                )

            mod = ProFormaParser._create_modification(mod_str, is_labile=True)
            modifications[-3].append(mod)  # Store labile modifications at position -3
            i = j + 1

        proforma_str = proforma_str[i:]

        if proforma_str.startswith("["):
            bracket_level = 0
            terminator_pos = -1

            # Find the terminal hyphen that's outside all brackets
            for i in range(len(proforma_str)):
                if proforma_str[i] == "[":
                    bracket_level += 1
                elif proforma_str[i] == "]":
                    bracket_level -= 1
                elif proforma_str[i] == "-" and bracket_level == 0:
                    terminator_pos = i
                    break

            if terminator_pos != -1:
                n_terminal_part = proforma_str[:terminator_pos]
                proforma_str = proforma_str[terminator_pos + 1 :]

                # Parse N-terminal modifications
                current_pos = 0
                while current_pos < len(n_terminal_part):
                    if n_terminal_part[current_pos] == "[":
                        bracket_depth = 1
                        end_pos = current_pos + 1

                        # Find matching closing bracket
                        while end_pos < len(n_terminal_part) and bracket_depth > 0:
                            if n_terminal_part[end_pos] == "[":
                                bracket_depth += 1
                            if n_terminal_part[end_pos] == "]":
                                bracket_depth -= 1
                            end_pos += 1

                        if bracket_depth == 0:
                            mod_string = n_terminal_part[current_pos + 1 : end_pos - 1]
                            n_term_mod = ProFormaParser._create_modification(
                                mod_string, is_terminal=True
                            )
                            modifications[-1].append(n_term_mod)

                        current_pos = end_pos
                    else:
                        current_pos += 1

        # Check for C-terminal modifications
        if "-" in proforma_str:
            bracket_level = 0
            terminator_pos = -1

            # Find the terminal hyphen that's outside all brackets, scanning from right to left
            for i in range(len(proforma_str) - 1, -1, -1):
                if proforma_str[i] == "]":
                    bracket_level += 1
                elif proforma_str[i] == "[":
                    bracket_level -= 1
                elif proforma_str[i] == "-" and bracket_level == 0:
                    terminator_pos = i
                    break

            if terminator_pos != -1:
                c_terminal_part = proforma_str[terminator_pos + 1 :]
                proforma_str = proforma_str[:terminator_pos]

                # Parse C-terminal modifications
                current_pos = 0
                while current_pos < len(c_terminal_part):
                    if c_terminal_part[current_pos] == "[":
                        bracket_depth = 1
                        end_pos = current_pos + 1

                        # Find matching closing bracket
                        while end_pos < len(c_terminal_part) and bracket_depth > 0:
                            if c_terminal_part[end_pos] == "[":
                                bracket_depth += 1
                            if c_terminal_part[end_pos] == "]":
                                bracket_depth -= 1
                            end_pos += 1

                        if bracket_depth == 0:
                            mod_string = c_terminal_part[current_pos + 1 : end_pos - 1]
                            c_term_mod = ProFormaParser._create_modification(
                                mod_string, is_terminal=True
                            )
                            modifications[-2].append(c_term_mod)

                        current_pos = end_pos
                    else:
                        current_pos += 1

        i = 0
        next_mod_is_gap = False
        range_stack = []  # Stack to keep track of ranges
        current_position = 0
        while i < len(proforma_str):
            char = proforma_str[i]
            if i + 1 < len(proforma_str) and proforma_str[i : i + 2] == "(?":
                closing_paren = proforma_str.find(")", i + 2)
                if closing_paren == -1:
                    raise ValueError("Unclosed sequence ambiguity parenthesis")

                # Extract ambiguous sequence
                ambiguous_seq = proforma_str[i + 2 : closing_paren]
                # Add to ambiguities list with current position
                sequence_ambiguities.append(
                    SequenceAmbiguity(ambiguous_seq, current_position)
                )

                # Skip past the ambiguity notation
                i = closing_paren + 1
                continue
            if char == "(":
                # Start of a range
                range_stack.append(len(base_sequence))
                i += 1
                continue

            elif char == ")":
                # End of a range
                if not range_stack:
                    raise ValueError("Unmatched closing parenthesis")

                range_start = range_stack.pop()
                range_end = len(base_sequence) - 1

                # Look for modification after the range
                j = i + 1
                while j < len(proforma_str) and proforma_str[j] == "[":
                    # Extract the modification that applies to the range
                    mod_start = j
                    bracket_count = 1
                    j += 1

                    while j < len(proforma_str) and bracket_count > 0:
                        if proforma_str[j] == "[":
                            bracket_count += 1
                        elif proforma_str[j] == "]":
                            bracket_count -= 1
                        j += 1

                    if bracket_count == 0:
                        # Get modification string and create modification
                        mod_str = proforma_str[mod_start + 1 : j - 1]
                        mod = ProFormaParser._create_modification(
                            mod_str,
                            in_range=True,
                            range_start=range_start,
                            range_end=range_end,
                        )

                        # Apply to all positions in range
                        for pos in range(mod.range_start, mod.range_end + 1):
                            modifications[pos].append(mod)

                i = j  # Skip past the modification
                continue

            elif char == "[":
                # Parse modification in square brackets
                bracket_count = 1
                j = i + 1
                while j < len(proforma_str) and bracket_count > 0:
                    if proforma_str[j] == "[":
                        bracket_count += 1
                    elif proforma_str[j] == "]":
                        bracket_count -= 1
                    j += 1

                if bracket_count > 0:
                    raise ValueError(f"Unclosed square bracket at position {i}")
                j -= 1
                if j == -1:
                    raise ValueError(f"Unclosed square bracket at position {i}")

                mod_str = proforma_str[i + 1 : j]
                if next_mod_is_gap:
                    mod = ProFormaParser._create_modification(mod_str, is_gap=True)
                    next_mod_is_gap = False
                # Check if this is a crosslink reference
                elif ProFormaParser.CROSSLINK_REF_PATTERN.match(mod_str):
                    mod = ProFormaParser._create_modification(
                        mod_str, is_crosslink_ref=True
                    )
                elif ProFormaParser.BRANCH_REF_PATTERN.match(mod_str):
                    mod = ProFormaParser._create_modification(
                        mod_str, is_branch_ref=True
                    )
                else:
                    # Check for crosslink or branch notation within the modification
                    crosslink_match = ProFormaParser.CROSSLINK_PATTERN.match(mod_str)
                    branch_match = ProFormaParser.BRANCH_PATTERN.match(mod_str)

                    if crosslink_match:
                        mod_base, crosslink_id = crosslink_match.groups()
                        mod = ProFormaParser._create_modification(
                            mod_str, crosslink_id=crosslink_id
                        )
                    elif branch_match:
                        mod_base = branch_match.group(1)
                        mod = ProFormaParser._create_modification(
                            mod_str, is_branch=True
                        )
                    else:
                        mod = ProFormaParser._create_modification(mod_str)

                # Add modification to the last amino acid
                if base_sequence:
                    modifications[len(base_sequence) - 1].append(mod)

                i = j + 1
            elif char == "{":
                # Parse ambiguous modification in curly braces
                j = proforma_str.find("}", i)
                if j == -1:
                    raise ValueError(f"Unclosed curly brace at position {i}")

                mod_str = proforma_str[i + 1 : j]
                mod = ProFormaParser._create_modification(mod_str, is_ambiguous=True)

                # Add ambiguous modification to the last amino acid
                if base_sequence:
                    modifications[len(base_sequence) - 1].append(mod)

                i = j + 1
            else:
                # Regular amino acid
                base_sequence += char
                is_gap = (
                    char == "X"
                    and i + 1 < len(proforma_str)
                    and proforma_str[i + 1] == "["
                )
                if is_gap:
                    next_mod_is_gap = True
                i += 1

        return base_sequence, modifications, global_mods, sequence_ambiguities

    @staticmethod
    def _create_modification(
        mod_str: str,
        is_terminal: bool = False,
        is_ambiguous: bool = False,
        is_labile: bool = False,
        is_unknown_position: bool = False,
        crosslink_id: Optional[str] = None,
        is_crosslink_ref: bool = False,
        is_branch: bool = False,
        is_branch_ref: bool = False,
        is_gap: bool = False,
        in_range: bool = False,
        range_start: Optional[int] = None,
        range_end: Optional[int] = None,
    ) -> Modification:
        """
        Create a Modification object from a ProForma modification string.

        Parameters
        ----------
        mod_str : str
            Modification string from ProForma notation
        is_terminal : bool
            Whether this is a terminal modification
        is_ambiguous : bool
            Whether this is an ambiguous modification
        crosslink_id : str, optional
            The crosslink identifier, if applicable
        is_crosslink_ref : bool
            Whether this is a crosslink reference
        is_branch : bool
            Whether this is a branch modification
        is_branch_ref : bool
            Whether this is a branch reference
        is_gap : bool
            Whether this is a gap modification


        Returns
        -------
        Modification
            Created modification object
        """
        mod_value = ModificationValue(mod_str)
        mod_type = "static"
        if is_terminal:
            mod_type = "terminal"
        elif is_ambiguous:
            mod_type = "ambiguous"
        elif is_labile:
            mod_type = "labile"
        elif is_unknown_position:
            mod_type = "unknown_position"
        elif crosslink_id or is_crosslink_ref:
            mod_type = "crosslink"
        elif is_branch or is_branch_ref:
            mod_type = "branch"
        elif is_gap:
            mod_type = "gap"

        ambiguity_match = re.match(r"(.+?)#([A-Za-z0-9]+)(?:\(([0-9.]+)\))?$", mod_str)
        ambiguity_ref_match = re.match(r"#([A-Za-z0-9]+)(?:\(([0-9.]+)\))?$", mod_str)
        ambiguity_group = None
        localization_score = None
        is_ambiguity_ref = False

        if ProFormaParser.MASS_SHIFT_PATTERN.match(mod_str) and "#" not in mod_str:
            mass_value = float(mod_str)
            if is_gap:
                return Modification(
                    mod_str, mass=mass_value, mod_type="gap", mod_value=mod_value
                )
            elif in_range:
                return Modification(
                    mod_str,
                    mass=mass_value,
                    mod_type="variable",
                    in_range=True,
                    range_start=range_start,
                    range_end=range_end,
                    mod_value=mod_value,
                )
            return Modification(
                f"Mass:{mod_str}",
                mass=mass_value,
                in_range=in_range,
                range_start=range_start,
                range_end=range_end,
                mod_value=mod_value,
            )

        if (
            "#" in mod_str
            and not is_crosslink_ref
            and not is_branch
            and not is_branch_ref
            and not crosslink_id
        ):
            if ambiguity_match and not ambiguity_match.group(2).startswith("XL"):
                mod_str = ambiguity_match.group(1)
                ambiguity_group = ambiguity_match.group(2)
                if ambiguity_match.group(3):  # Score is present
                    localization_score = float(ambiguity_match.group(3))
                mod = Modification(
                    mod_str,
                    mod_type="ambiguous",
                    ambiguity_group=ambiguity_group,
                    is_ambiguity_ref=False,
                    in_range=in_range,
                    range_start=range_start,
                    range_end=range_end,
                    localization_score=localization_score,
                    mod_value=mod_value,
                )
                return mod
            elif ambiguity_ref_match and not ambiguity_ref_match.group(1).startswith(
                "XL"
            ):
                ambiguity_group = ambiguity_ref_match.group(1)
                if ambiguity_ref_match.group(2):  # Score is present
                    localization_score = float(ambiguity_ref_match.group(2))
                mod = Modification(
                    "",
                    mod_type="ambiguous",
                    ambiguity_group=ambiguity_group,
                    is_ambiguity_ref=True,
                    in_range=in_range,
                    range_start=range_start,
                    range_end=range_end,
                    localization_score=localization_score,
                    mod_value=mod_value,
                )
                return mod

        # Create the modification with appropriate attributes
        return Modification(
            mod_str,
            mod_type=mod_type,
            crosslink_id=crosslink_id,
            is_crosslink_ref=is_crosslink_ref,
            is_branch=is_branch,
            is_branch_ref=is_branch_ref,
            in_range=in_range,
            range_start=range_start,
            range_end=range_end,
            mod_value=mod_value,
        )


class SequenceAmbiguity:
    """Represents ambiguity in the amino acid sequence."""

    def __init__(self, value: str, position: int):
        """Initialize a sequence ambiguity.

        Parameters
        ----------
        value : str
            The ambiguous sequence possibilities
        position : int
            The position in the sequence where ambiguity occurs
        """
        self.value = value
        self.position = position

    def __repr__(self) -> str:
        return f"SequenceAmbiguity(value='{self.value}', position={self.position})"
