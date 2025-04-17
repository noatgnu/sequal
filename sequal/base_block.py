from abc import ABC, abstractmethod
from typing import Any, Dict, Optional


class BaseBlock(ABC):
    """
    Base class for biochemical building blocks with position and mass properties.

    This abstract base class provides core functionality for blocks such as
    amino acids, modifications, and other biochemical components.
    """

    def __init__(
        self,
        value: str,
        position: Optional[int] = None,
        branch: bool = False,
        mass: Optional[float] = None,
    ):
        """
        Initialize a BaseBlock object.

        Parameters
        ----------
        value : str
            The identifier of the block.
        position : int, optional
            The position of the block within a chain.
        branch : bool, optional
            Indicates whether this block is a branch of another block.
        mass : float, optional
            The mass of the block in Daltons.
        """
        self._value = value
        self._position = position
        self._branch = branch
        self._mass = mass
        self._extra = None

    @property
    def value(self) -> str:
        """Get the identifier of the block."""
        return self._value

    @property
    def position(self) -> Optional[int]:
        """Get the position of the block."""
        return self._position

    @position.setter
    def position(self, position: Optional[int]) -> None:
        """Set the position of the block."""
        self._position = position

    @property
    def branch(self) -> bool:
        """Check if the block is a branch."""
        return self._branch

    @property
    def mass(self) -> Optional[float]:
        """Get the mass of the block."""
        return self._mass

    @mass.setter
    def mass(self, mass: Optional[float]) -> None:
        """Set the mass of the block."""
        self._mass = mass

    @property
    def extra(self) -> Any:
        """Get extra information associated with the block."""
        return self._extra

    @extra.setter
    def extra(self, value: Any) -> None:
        """Set extra information for the block."""
        self._extra = value

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the block to a dictionary representation.

        Returns
        -------
        Dict[str, Any]
            Dictionary containing the block's attributes.
        """
        return {
            "value": self._value,
            "position": self._position,
            "branch": self._branch,
            "mass": self._mass,
            "extra": self._extra,
        }

    def __eq__(self, other) -> bool:
        """Check if two blocks are equal."""
        if not isinstance(other, BaseBlock):
            return False
        return (
            self._value == other.value
            and self._position == other.position
            and self._branch == other.branch
        )

    def __hash__(self) -> int:
        """Generate a hash for the block."""
        return hash((self._value, self._position, self._branch))

    def __str__(self) -> str:
        """Return a string representation of the block."""
        return self._value

    def __repr__(self) -> str:
        """Return a detailed string representation of the block."""
        return f"{self.__class__.__name__}(value='{self._value}', position={self._position})"
