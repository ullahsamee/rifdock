import unittest
import NEST_wrap
import util_tests

class TestSequenceFunctions(unittest.TestCase):
    def test_sanity(self):
        assert 1 == 1

class TestSchemeUtilities(unittest.TestCase):
	def test_dilation(self):
		assert util_tests.test_zorder()


