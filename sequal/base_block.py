# Base block object to be used for blocks that can have position and mass properties
class BaseBlock:
    def __init__(self, value, position, branch=False, mass=None):
        """
        Initialize a BaseBlock object.

        :param value: str
            The name of the block.
        :param position: int
            The position of the block within a chain.
        :param branch: bool, optional
            Indicates whether this block is a branch of another block (default is False).
        :param mass: float, optional
            The mass of the block (default is None).
        """
        self.value = value
        self.position = position
        self.branch = branch
        self.mass = mass
        self.extra = None

    def __str__(self):
        """
        Return a string representation of the block.

        :return: str
            The name of the block.
        """
        return self.value

    def __repr__(self):
        """
        Return a string representation of the block for debugging.

        :return: str
            The name of the block.
        """
        return self.value
