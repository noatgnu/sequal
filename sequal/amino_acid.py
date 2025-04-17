from typing import List, Optional, Dict, Any, Union
from sequal.base_block import BaseBlock
from sequal.modification import Modification
from sequal.resources import AA_mass


class AminoAcid(BaseBlock):
    """
    Represents an amino acid block that can carry position, modifications, and amino acid value.

    Inherits from the BaseBlock class and adds functionality specific to amino acids, such as
    handling modifications and inferring mass from a predefined dictionary.
    """

    def __init__(self,
                 value: str,
                 position: Optional[int] = None,
                 mass: Optional[float] = None):
        """
        Initialize an AminoAcid object.

        Parameters
        ----------
        value : str
            The amino acid one letter or three letter code.
        position : int, optional
            The position of this amino acid in a sequence.
        mass : float, optional
            The mass of the amino acid. If not provided, inferred from AA_mass dictionary.
        """
        if value not in AA_mass and not mass:
            raise ValueError(f"Unknown amino acid '{value}' and no mass provided")

        inferred_mass = AA_mass.get(value) if not mass else mass
        super().__init__(value, position, branch=False, mass=inferred_mass)
        self._mods: List[Modification] = []

    @property
    def mods(self) -> List[Modification]:
        """Get the list of modifications applied to this amino acid."""
        return self._mods.copy()  # Return a copy to prevent direct modification

    def add_modification(self, mod: Modification) -> None:
        """
        Add a modification to this amino acid.

        Parameters
        ----------
        mod : Modification
            The modification to add.
        """
        self._mods.append(mod)

    def set_modification(self, mod: Modification) -> None:
        """
        Add a modification to this amino acid (legacy method).

        Parameters
        ----------
        mod : Modification
            The modification to add.
        """
        self.add_modification(mod)

    def remove_modification(self, mod: Union[Modification, str]) -> bool:
        """
        Remove a modification from this amino acid.

        Parameters
        ----------
        mod : Modification or str
            The modification or modification value to remove.

        Returns
        -------
        bool
            True if modification was removed, False if not found.
        """
        if isinstance(mod, str):
            for i, existing_mod in enumerate(self._mods):
                if existing_mod.value == mod:
                    self._mods.pop(i)
                    return True
        else:
            if mod in self._mods:
                self._mods.remove(mod)
                return True
        return False

    def has_modification(self, mod: Union[Modification, str]) -> bool:
        """
        Check if this amino acid has a specific modification.

        Parameters
        ----------
        mod : Modification or str
            The modification or modification value to check for.

        Returns
        -------
        bool
            True if the modification exists, False otherwise.
        """
        if isinstance(mod, str):
            return any(m.value == mod for m in self._mods)
        return mod in self._mods

    def get_total_mass(self) -> float:
        """
        Calculate the total mass including all modifications.

        Returns
        -------
        float
            The total mass of the amino acid with all modifications.
        """
        total = self.mass or 0
        for mod in self._mods:
            if mod.mass:
                total += mod.mass
        return total

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the amino acid to a dictionary representation.

        Returns
        -------
        Dict[str, Any]
            Dictionary containing the amino acid's attributes including modifications.
        """
        result = super().to_dict()
        result['mods'] = [mod.to_dict() for mod in self._mods]
        result['total_mass'] = self.get_total_mass()
        return result

    def __eq__(self, other) -> bool:
        """Check if two amino acids are equal including their modifications."""
        if not super().__eq__(other):
            return False
        if not isinstance(other, AminoAcid):
            return False
        if len(self._mods) != len(other._mods):
            return False
        # Sort mods by value for comparison
        self_mods = sorted(self._mods, key=lambda m: m.value)
        other_mods = sorted(other._mods, key=lambda m: m.value)
        return all(a == b for a, b in zip(self_mods, other_mods))

    def __hash__(self) -> int:
        """Generate a hash for the amino acid including modifications."""
        mod_hash = hash(tuple(sorted(m.value for m in self._mods)))
        return hash((super().__hash__(), mod_hash))

    def __str__(self) -> str:
        """Return a string representation with modifications."""
        s = self.value
        for i in self._mods:
            s += f"[{i.value}]"
        return s

    def __repr__(self) -> str:
        """Return a detailed string representation for debugging."""
        mod_str = ", ".join(repr(m) for m in self._mods)
        return f"AminoAcid(value='{self.value}', position={self.position}, mods=[{mod_str}])"