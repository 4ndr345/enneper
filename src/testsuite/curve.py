#!/usr/bin/env python
# -*- coding: utf-8 -*-

#***************************************************************************
#*   Copyright (C) 2012 by Andreas Kührmann [andreas.kuehrmann@gmail.com]  *
#*                                                                         *
#*   This program is free software; you can redistribute it and/or modify  *
#*   it under the terms of the GNU General Public License as published by  *
#*   the Free Software Foundation; either version 3 of the License, or     *
#*   (at your option) any later version.                                   *
#*                                                                         *
#*   This program is distributed in the hope that it will be useful,       *
#*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
#*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
#*   GNU General Public License for more details.                          *
#*                                                                         *
#*   You should have received a copy of the GNU General Public License     *
#*   along with this program; if not, write to the                         *
#*   Free Software Foundation, Inc.,                                       *
#*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
#***************************************************************************


import unittest

from enneper.curve import Curve


CTRL_PNTS = [[0, 0, 1], [4, 4, 4], [3, 2, 1], [4, 1, 1], [5, -1, 1]]
KNOTS = [0, 0, 0, 1, 2, 3, 3, 3]
DEG = 2


class TestCurve(unittest.TestCase):

    def test_resize(self):
        curve = Curve()
        curve.resize(n=17, dim=4, deg=3)
        self.assertEqual(curve.ctrl_pnts.shape, (17, 4))
        self.assertEqual(curve.knots.size, 17+3+1)
        self.assertEqual(curve.deg, 3)


if __name__ == '__main__':
    unittest.main()
