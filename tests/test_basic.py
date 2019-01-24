# -*- coding: utf-8 -*-

from .context import sirpy2

import unittest


class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_fn2dt(self):
        fname = "msfa-a-NAm07-181-185.sir"
        result = sirpy2.fn2dt(fname)
        self.assertEqual(result.year, 2007)
        self.assertEqual(result.month, 7)
        self.assertEqual(result.day, 2)


if __name__ == '__main__':
    unittest.main()
