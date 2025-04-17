import re
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from sequal.modification import Modification


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

        # Check for terminal modifications
        n_term_mod_str = None
        c_term_mod_str = None

        # Parse terminal modifications if present
        if match := ProFormaParser.TERMINAL_PATTERN.match(proforma_str):
            n_term_mod_str, proforma_str, c_term_mod_str = match.groups()
        elif match := ProFormaParser.N_TERMINAL_PATTERN.match(proforma_str):
            n_term_mod_str, proforma_str = match.groups()
        elif match := ProFormaParser.C_TERMINAL_PATTERN.match(proforma_str):
            proforma_str, c_term_mod_str = match.groups()

        # Add terminal modifications
        if n_term_mod_str:
            n_term_mod = ProFormaParser._create_modification(
                n_term_mod_str, is_terminal=True
            )
            modifications[-1].append(n_term_mod)

        if c_term_mod_str:
            c_term_mod = ProFormaParser._create_modification(
                c_term_mod_str, is_terminal=True
            )
            modifications[-2].append(c_term_mod)

        # Process the main sequence
        i = 0
        next_mod_is_gap = False
        while i < len(proforma_str):
            char = proforma_str[i]

            if char == "[":
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
                            mod_base, crosslink_id=crosslink_id
                        )
                    elif branch_match:
                        mod_base = branch_match.group(1)
                        mod = ProFormaParser._create_modification(
                            mod_base, is_branch=True
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

        return base_sequence, modifications

    @staticmethod
    def _create_modification(
        mod_str: str,
        is_terminal: bool = False,
        is_ambiguous: bool = False,
        crosslink_id: Optional[str] = None,
        is_crosslink_ref: bool = False,
        is_branch: bool = False,
        is_branch_ref: bool = False,
        is_gap: bool = False,
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
        # Handle mass shift notation

        # Determine modification type
        mod_type = "static"  # Default
        if is_terminal:
            mod_type = "terminal"
        elif is_ambiguous:
            mod_type = "ambiguous"
        elif crosslink_id or is_crosslink_ref:
            mod_type = "crosslink"
        elif is_branch or is_branch_ref:
            mod_type = "branch"
        elif is_gap:
            mod_type = "gap"

        if ProFormaParser.MASS_SHIFT_PATTERN.match(mod_str):
            mass_value = float(mod_str)
            if is_gap:
                return Modification(mod_str, mass=mass_value, mod_type="gap")
            return Modification(f"Mass:{mod_str}", mass=mass_value)

        # Create the modification with appropriate attributes
        return Modification(
            mod_str,
            mod_type=mod_type,
            crosslink_id=crosslink_id,
            is_crosslink_ref=is_crosslink_ref,
            is_branch=is_branch,
            is_branch_ref=is_branch_ref,
        )
