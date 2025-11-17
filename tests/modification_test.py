import unittest

from sequal.modification import Modification


class ModificationTestCase(unittest.TestCase):
    """Test cases for the Modification class."""

    def test_find_positions(self):
        """Test finding modification positions with regex."""
        mod = Modification("HexNAc", regex_pattern="N[^P][S|T]")
        for ps, pe in mod.find_positions("TESNEST"):
            self.assertEqual(
                ps, 3, "HexNAc is at expected index position {}".format(ps)
            )


class ModificationMapTestCase(unittest.TestCase):
    """Test cases for the ModificationMap class."""

    def test_map_creation(self):
        """Test the creation of a modification map."""
        pass


if __name__ == "__main__":
    unittest.main()
