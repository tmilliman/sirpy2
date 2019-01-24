# -*- coding: utf-8 -*-

import sirpy2

def test_extractsirhead():
        file = "tests/msfa-a-NAm07-181-185.sir"
        sir = sirpy2.loadsir(file)
        sh = sirpy2.extractsirhead(sir[1])
        assert sh['region'] == 'NAm'
        assert sh['year'] == 2007
        assert sh['nsx'] == 1890
        assert sh['nsy'] == 1150
