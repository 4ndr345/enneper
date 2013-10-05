#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ***************************************************************************
# *   Copyright (C) 2013 by Andreas Kührmann [andreas.kuehrmann@gmail.com]  *
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU General Public License as published by  *
# *   the Free Software Foundation; either version 3 of the License, or     *
# *   (at your option) any later version.                                   *
# *                                                                         *
# *   This program is distributed in the hope that it will be useful,       *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU General Public License for more details.                          *
# *                                                                         *
# *   You should have received a copy of the GNU General Public License     *
# *   along with this program; if not, write to the                         *
# *   Free Software Foundation, Inc.,                                       *
# *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
# ***************************************************************************


import os
import unittest

import numpy as np

import enneper.curve as enneper


DIR = os.path.join(os.path.dirname(__file__))

# The NURBS Book 2nd edition: example 4.1 (page 122 ff.)
CTRL_PNTS = [[0, 0, 1], [4, 4, 4], [3, 2, 1], [4, 1, 1], [5, -1, 1]]
KNOTS = [0, 0, 0, 1, 2, 3, 3, 3]
DEG = 2


class TestCurve(unittest.TestCase):

    def test_designated_initializer(self):
        curve = enneper.Curve(CTRL_PNTS, KNOTS, DEG)
        np.testing.assert_equal(curve.ctrl_pnts, CTRL_PNTS)
        np.testing.assert_equal(curve.knots, KNOTS)
        self.assertEqual(curve.deg, DEG)

    def test_from_curve_constructor(self):
        original = enneper.Curve(CTRL_PNTS, KNOTS, DEG)
        copy = enneper.Curve.from_curve(original)
        np.testing.assert_equal(original.ctrl_pnts, copy.ctrl_pnts)
        np.testing.assert_equal(original.knots, copy.knots)
        self.assertEqual(original.deg, copy.deg)
        self.assertIsNot(original.ctrl_pnts, copy.ctrl_pnts)
        self.assertIsNot(original.knots, copy.knots)

    def test_from_json_constructor(self):
        curve = enneper.Curve.from_json(os.path.join(DIR, 'curve.json'))
        np.testing.assert_equal(curve.ctrl_pnts, CTRL_PNTS)
        np.testing.assert_equal(curve.knots, KNOTS)
        self.assertEqual(curve.deg, DEG)

    def test_evaluate_at(self):
        curve = enneper.Curve(CTRL_PNTS, KNOTS, DEG)
        np.testing.assert_equal(curve.evaluate_at(1), [3.5, 3., 2.5])


if __name__ == '__main__':
    unittest.main()
