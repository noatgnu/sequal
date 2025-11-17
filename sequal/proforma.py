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
    CHARGE_PATTERN = re.compile(r"/(-?\d+)(?:\[(.*?)\])?(?=-|\Z)")

    @staticmethod
    def _find_balanced_paren(s: str, start: int) -> int:
        """
        Finds the closing parenthesis matching the opening one at start.

        Handles nested parentheses correctly for ProForma 2.1 named entities.

        Args:
            s (str): The string to search.
            start (int): The position of the opening parenthesis.

        Returns:
            int: Position of the matching closing parenthesis, or -1 if not found.
        """
        if start >= len(s) or s[start] != "(":
            return -1

        count = 1
        i = start + 1
        while i < len(s) and count > 0:
            if s[i] == "(":
                count += 1
            elif s[i] == ")":
                count -= 1
            i += 1

        return i - 1 if count == 0 else -1

    @staticmethod
    def parse(
        proforma_str: str,
    ) -> Tuple[
        str,
        Dict[int, List[Modification]],
        List[GlobalModification],
        List["SequenceAmbiguity"],
        Optional[Tuple[int, str]],
        Dict[str, Optional[str]],
    ]:
        """
        Parses a ProForma string into a base sequence and modifications.

        Args:
            proforma_str (str): ProForma formatted peptide string.

        Returns:
            Tuple[str, Dict[int, List[Modification]], List[GlobalModification], List["SequenceAmbiguity"], Optional[Tuple[int, str]], Dict[str, Optional[str]]]:
                A tuple containing:
                    - The base sequence.
                    - A dictionary of modifications keyed by position.
                    - A list of global modifications.
                    - A list of sequence ambiguities.
                    - A tuple of charge and ionic species.
                    - A dictionary of names.
        """
        base_sequence = ""
        modifications = defaultdict(list)
        global_mods = []
        sequence_ambiguities = []

        names = {
            "peptidoform_name": None,
            "peptidoform_ion_name": None,
            "compound_ion_name": None,
        }

        if proforma_str.startswith("(>>>"):
            close_idx = ProFormaParser._find_balanced_paren(proforma_str, 0)
            if close_idx != -1:
                names["compound_ion_name"] = proforma_str[4:close_idx]
                proforma_str = proforma_str[close_idx + 1 :]

        if proforma_str.startswith("(>>"):
            close_idx = ProFormaParser._find_balanced_paren(proforma_str, 0)
            if close_idx != -1:
                names["peptidoform_ion_name"] = proforma_str[3:close_idx]
                proforma_str = proforma_str[close_idx + 1 :]

        if proforma_str.startswith("(>"):
            close_idx = ProFormaParser._find_balanced_paren(proforma_str, 0)
            if close_idx != -1:
                names["peptidoform_name"] = proforma_str[2:close_idx]
                proforma_str = proforma_str[close_idx + 1 :]

        while proforma_str.startswith("<"):
            if proforma_str.startswith("<["):
                closing_bracket = proforma_str.find("]", 2)
                if closing_bracket == -1:
                    raise ValueError("Unclosed square bracket in global modification")
                end_bracket = proforma_str.find(">", closing_bracket)
            else:
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

                target_list = []
                has_terminal_targets = False
                for target in targets.split(","):
                    target = target.strip()

                    if ":" in target:
                        has_terminal_targets = True
                        term_part, aa_part = target.split(":", 1)
                        target_list.append(
                            {
                                "type": "terminal_specific",
                                "terminal": term_part,
                                "amino_acid": aa_part,
                            }
                        )
                    elif target in ["N-term", "C-term"]:
                        has_terminal_targets = True
                        target_list.append({"type": "terminal", "terminal": target})
                    else:
                        if has_terminal_targets:
                            target_list.append(
                                {"type": "amino_acid", "residue": target}
                            )
                        else:
                            target_list.append(target)

                global_mods.append(GlobalModification(mod_value, target_list, "fixed"))
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

        proforma_str, charge, species = ProFormaParser.parse_charge_info(proforma_str)

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

        return (
            base_sequence,
            modifications,
            global_mods,
            sequence_ambiguities,
            (charge, species),
            names,
        )

    @staticmethod
    def parse_charge_info(proforma_str):
        """
        Parses charge information from a ProForma string.

        Args:
            proforma_str (str): The ProForma string to parse.

        Returns:
            tuple: A tuple containing:
                - The updated ProForma string.
                - The charge value.
                - The ionic species.
        """
        if "/" not in proforma_str:
            return proforma_str, None, None

        charge_pos = -1
        bracket_level = 0

        for i, char in enumerate(proforma_str):
            if char in "[(":
                bracket_level += 1
            elif char in "])":
                bracket_level -= 1
            elif char == "/" and bracket_level == 0:
                charge_pos = i
                break

        if charge_pos == -1:
            return proforma_str, None, None

        before_charge = proforma_str[:charge_pos]
        after_charge = proforma_str[charge_pos + 1 :]

        i = 0
        sign = 1

        if i < len(after_charge) and after_charge[i] == "-":
            sign = -1
            i += 1

        start_digit = i
        while i < len(after_charge) and after_charge[i].isdigit():
            i += 1

        if start_digit == i:
            return proforma_str, None, None

        charge_value = int(after_charge[start_digit:i]) * sign

        remaining = after_charge[i:]
        ionic_species = None

        if remaining and remaining[0] == "[":
            # Find the matching closing bracket
            bracket_level = 1
            end_pos = 0

            for j, char in enumerate(remaining[1:], 1):
                if char == "[":
                    bracket_level += 1
                elif char == "]":
                    bracket_level -= 1

                if bracket_level == 0:
                    end_pos = j
                    break

            if end_pos > 0:
                ionic_species = remaining[1:end_pos]
                remaining = remaining[end_pos + 1 :]

        result_str = before_charge
        if remaining:
            result_str += remaining

        return result_str, charge_value, ionic_species

    @staticmethod
    def _parse_placement_controls(mod_str: str) -> Tuple[str, Dict]:
        """
        Parses placement control tags from a modification string.

        Args:
            mod_str (str): Modification string that may contain placement controls.

        Returns:
            Tuple[str, Dict]: A tuple containing:
                - The cleaned modification string with controls removed.
                - A dictionary of placement controls.
        """
        controls = {
            "position_constraint": None,
            "limit_per_position": 1,
            "colocalize_known": False,
            "colocalize_unknown": False,
        }

        if "|" not in mod_str:
            return mod_str, controls

        parts = mod_str.split("|")
        mod_parts = []

        for i, part in enumerate(parts):
            part = part.strip()

            if part.startswith("Position:"):
                positions_str = part[9:]
                controls["position_constraint"] = [
                    p.strip() for p in positions_str.split(",")
                ]
            elif part.startswith("Limit:"):
                try:
                    controls["limit_per_position"] = int(part[6:])
                except ValueError:
                    mod_parts.append(part)
            elif part in ("CoMKP", "ColocaliseModificationsOfKnownPosition"):
                controls["colocalize_known"] = True
            elif part in ("CoMUP", "ColocaliseModificationsOfUnknownPosition"):
                controls["colocalize_unknown"] = True
            else:
                mod_parts.append(part)

        cleaned_mod_str = "|".join(mod_parts) if mod_parts else parts[0]
        return cleaned_mod_str, controls

    @staticmethod
    def _is_ion_type_modification(mod_str: str) -> bool:
        """
        Checks if a modification represents an ion type.

        Args:
            mod_str (str): Modification string to check.

        Returns:
            bool: True if the modification is an ion type, False otherwise.
        """
        mod_str_lower = mod_str.lower()

        # Check for -type-ion suffix
        if mod_str_lower.endswith("-type-ion"):
            return True

        # Check for known Unimod ion type IDs
        # UNIMOD:140 = a-type-ion, UNIMOD:2132 = b-type-ion, etc.
        ion_type_unimod_ids = {
            "140",  # a-type-ion
            "2132",  # b-type-ion
            "4",  # c-type-ion
            "24",  # x-type-ion
            "2133",  # y-type-ion
            "23",  # z-type-ion
        }

        if mod_str.startswith(("UNIMOD:", "U:")):
            unimod_id = mod_str.split(":")[1]
            if unimod_id in ion_type_unimod_ids:
                return True

        return False

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
        Creates a Modification object from a ProForma modification string.

        Args:
            mod_str (str): Modification string from ProForma notation.
            is_terminal (bool, optional): Whether this is a terminal modification.
                Defaults to False.
            is_ambiguous (bool, optional): Whether this is an ambiguous
                modification. Defaults to False.
            is_labile (bool, optional): Whether this is a labile modification.
                Defaults to False.
            is_unknown_position (bool, optional): Whether this is a modification
                of unknown position. Defaults to False.
            crosslink_id (str, optional): The crosslink identifier, if applicable.
                Defaults to None.
            is_crosslink_ref (bool, optional): Whether this is a crosslink
                reference. Defaults to False.
            is_branch (bool, optional): Whether this is a branch modification.
                Defaults to False.
            is_branch_ref (bool, optional): Whether this is a branch reference.
                Defaults to False.
            is_gap (bool, optional): Whether this is a gap modification.
                Defaults to False.
            in_range (bool, optional): Whether this modification is in a range.
                Defaults to False.
            range_start (int, optional): The start of the range. Defaults to None.
            range_end (int, optional): The end of the range. Defaults to None.

        Returns:
            Modification: The created modification object.
        """
        cleaned_mod_str, placement_controls = ProFormaParser._parse_placement_controls(
            mod_str
        )
        mod_value = ModificationValue(cleaned_mod_str)
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

        ambiguity_match = re.match(
            r"(.+?)#([A-Za-z0-9]+)(?:\(([0-9.]+)\))?$", cleaned_mod_str
        )
        ambiguity_ref_match = re.match(
            r"#([A-Za-z0-9]+)(?:\(([0-9.]+)\))?$", cleaned_mod_str
        )
        ambiguity_group = None
        localization_score = None
        is_ambiguity_ref = False

        if (
            ProFormaParser.MASS_SHIFT_PATTERN.match(cleaned_mod_str)
            and "#" not in cleaned_mod_str
        ):
            mass_value = float(cleaned_mod_str)
            if is_gap:
                return Modification(
                    cleaned_mod_str,
                    mass=mass_value,
                    mod_type="gap",
                    mod_value=mod_value,
                    **placement_controls,
                )
            elif in_range:
                return Modification(
                    cleaned_mod_str,
                    mass=mass_value,
                    mod_type="variable",
                    in_range=True,
                    range_start=range_start,
                    range_end=range_end,
                    mod_value=mod_value,
                    **placement_controls,
                )
            return Modification(
                f"Mass:{cleaned_mod_str}",
                mass=mass_value,
                in_range=in_range,
                range_start=range_start,
                range_end=range_end,
                mod_value=mod_value,
                **placement_controls,
            )

        if (
            "#" in cleaned_mod_str
            and not is_crosslink_ref
            and not is_branch
            and not is_branch_ref
            and not crosslink_id
        ):
            if ambiguity_match and not ambiguity_match.group(2).startswith("XL"):
                cleaned_mod_str = ambiguity_match.group(1)
                ambiguity_group = ambiguity_match.group(2)
                if ambiguity_match.group(3):  # Score is present
                    localization_score = float(ambiguity_match.group(3))
                mod = Modification(
                    cleaned_mod_str,
                    mod_type="ambiguous",
                    ambiguity_group=ambiguity_group,
                    is_ambiguity_ref=False,
                    in_range=in_range,
                    range_start=range_start,
                    range_end=range_end,
                    localization_score=localization_score,
                    mod_value=mod_value,
                    **placement_controls,
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
                    **placement_controls,
                )
                return mod

        # Check if this is an ion type modification (ProForma 2.1 Section 11.6)
        is_ion = ProFormaParser._is_ion_type_modification(cleaned_mod_str)

        # Create the modification with appropriate attributes
        return Modification(
            cleaned_mod_str,
            mod_type=mod_type,
            crosslink_id=crosslink_id,
            is_crosslink_ref=is_crosslink_ref,
            is_branch=is_branch,
            is_branch_ref=is_branch_ref,
            in_range=in_range,
            range_start=range_start,
            range_end=range_end,
            mod_value=mod_value,
            is_ion_type=is_ion,
            **placement_controls,
        )


class SequenceAmbiguity:
    """
    Represents ambiguity in the amino acid sequence.

    Attributes:
        value (str): The ambiguous sequence possibilities.
        position (int): The position in the sequence where ambiguity occurs.
    """

    def __init__(self, value: str, position: int):
        """
        Initializes a sequence ambiguity.

        Args:
            value (str): The ambiguous sequence possibilities.
            position (int): The position in the sequence where ambiguity occurs.
        """
        self.value = value
        self.position = position

    def __repr__(self) -> str:
        return f"SequenceAmbiguity(value='{self.value}', position={self.position})"
