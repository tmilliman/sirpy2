# -*- coding: utf-8 -*-

import sirpy2


def test_loadsir():
    fname = "tests/msfa-V-NAm07-181-185.sir"
    image, head, descrip, iaopt = sirpy2.loadsir(fname)
    assert image.shape == (1150, 1890)
    assert image.dtype == "float64"
    assert head[2] == -92.5
    assert head[3] == 45.0
    assert head[7] == -4200.0
    assert head[8] == -2300.0
