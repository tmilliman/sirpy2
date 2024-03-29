# -*- coding: utf-8 -*-

import sirpy2


def test_fn2dt():
    fname = "msfa-V-NAm07-181-185.sir"
    result = sirpy2.fn2dt(fname, date_flag="center")
    assert result.year == 2007
    assert result.month == 7
    assert result.day == 2
    result = sirpy2.fn2dt(fname)
    assert result.year == 2007
    assert result.month == 6
    assert result.day == 30


def test_parseFilename():
    fname = "msfa-V-NAm07-181-185.sir"
    fparts = sirpy2.parseFilename(fname)
    assert fparts["product"] == "msfa"
    assert fparts["itype"] == "V"
    assert fparts["region"] == "NAm"
    assert fparts["year"] == 2007
