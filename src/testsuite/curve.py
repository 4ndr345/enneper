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

import numpy as np

from enneper.curve import Curve


# The NURBS Book 2nd edition: example 4.1 (page 122 ff.)
CTRL_PNTS = [[0, 0, 1], [4, 4, 4], [3, 2, 1], [4, 1, 1], [5, -1, 1]]
KNOTS = [0, 0, 0, 1, 2, 3, 3, 3]
DIM = len(CTRL_PNTS[0])
N = len(CTRL_PNTS)
DEG = 2


class TestCurve(unittest.TestCase):

    def test_resize(self):
        curve = Curve()
        curve.resize(N, DIM, DEG)
        self.assertEqual(curve.ctrl_pnts.shape, (N, DIM))
        self.assertEqual(curve.knots.size, N + DEG + 1)
        self.assertEqual(curve.deg, DEG)

    def test_from_curve(self):
        original = Curve()
        original.resize(N, DIM, DEG)
        original.ctrl_pnts[:] = CTRL_PNTS
        original.knots[:] = KNOTS
        copy = Curve.from_curve(original)
        np.testing.assert_equal(original.ctrl_pnts, copy.ctrl_pnts)
        np.testing.assert_equal(original.knots, copy.knots)
        self.assertEqual(original.deg, copy.deg)
        self.assertIsNot(original.ctrl_pnts, copy.ctrl_pnts)
        self.assertIsNot(original.knots, copy.knots)

    def test_from_parameters(self):
        curve = Curve.from_ctrl_pnts_knots_and_deg(CTRL_PNTS, KNOTS, DEG)
        np.testing.assert_equal(curve.ctrl_pnts, CTRL_PNTS)
        np.testing.assert_equal(curve.knots, KNOTS)
        self.assertEqual(curve.deg, DEG)

    def test_evaluate_at(self):
        curve = Curve.from_ctrl_pnts_knots_and_deg(CTRL_PNTS, KNOTS, DEG)
        actual = curve.evaluate_at(1)
        desired = np.asarray([3.5, 3, 2.5])
        np.testing.assert_equal(actual, desired)


if __name__ == '__main__':
    unittest.main()
