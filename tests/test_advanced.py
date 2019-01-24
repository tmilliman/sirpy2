# -*- coding: utf-8 -*-

from .context import sirpy2

import unittest

class AdvancedTestSuite(unittest.TestCase):
    """Advanced test cases."""

    def test_extractsirhead(self):
        file = "tests/msfa-a-NAm07-181-185.sir"
        sir = sirpy2.loadsir(file)
        sh = sirpy2.extractsirhead(sir[1])
        self.assertEqual(sh['region'], 'NAm')
        self.assertEqual(sh['year'], 2007)
        self.assertEqual(sh['nsx'], 1890)
        self.assertEqual(sh['nsy'], 1150)        


if __name__ == '__main__':
    unittest.main()
