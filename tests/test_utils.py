import unittest

from higgs_dna.utils import setup_logger


class TestUtils(unittest.TestCase):
    """
    Test main utilities provided in the utils submodule
    """
    def setUp(self):
        pass

    def test_logger(self):
        """
        Test main aspects of setup_logger function
        """
        null_info_level = "NOT_ALLOWED"
        self.assertRaises(ValueError, setup_logger, level=null_info_level)


if __name__ == '__main__':
    unittest.main()
