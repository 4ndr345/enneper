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


# standard packages
import cStringIO
import json
import unittest

# 3rd party packages
import numpy as np

# project packages
import enneper


# The NURBS Book 2nd edition: example 4.1 (page 122 ff.)
CTRL_PNTS = [[0, 0, 1], [4, 4, 4], [3, 2, 1], [4, 1, 1], [5, -1, 1]]
KNOTS = [0, 0, 0, 1, 2, 3, 3, 3]
DEG = 2


class TestCurve(unittest.TestCase):

    def test_designated_initializer(self):

        # test
        curve = enneper.Curve(CTRL_PNTS, KNOTS, DEG)
        np.testing.assert_equal(curve.ctrl_pnts, CTRL_PNTS)
        np.testing.assert_equal(curve.knots, KNOTS)
        self.assertEqual(curve.deg, DEG)

    def test_from_curve_constructor(self):

        # construct test data
        original = enneper.Curve(CTRL_PNTS, KNOTS, DEG)

        # test
        copy = enneper.Curve.from_curve(original)
        np.testing.assert_equal(original.ctrl_pnts, copy.ctrl_pnts)
        np.testing.assert_equal(original.knots, copy.knots)
        self.assertEqual(original.deg, copy.deg)
        self.assertIsNot(original.ctrl_pnts, copy.ctrl_pnts)
        self.assertIsNot(original.knots, copy.knots)

    def test_from_json_constructor(self):

        # construct test data
        data = dict(ctrl_pnts=CTRL_PNTS, knots=KNOTS, deg=DEG)
        flo = cStringIO.StringIO()
        json.dump(data, flo)
        flo.seek(0)

        # test 
        curve = enneper.Curve.from_json(flo)
        np.testing.assert_equal(curve.ctrl_pnts, CTRL_PNTS)
        np.testing.assert_equal(curve.knots, KNOTS)
        self.assertEqual(curve.deg, DEG)

        # close and discard memory buffer
        flo.close()

    def test_evaluate_at(self):

        # The NURBS Book 2nd edition: example 4.1 (page 122 ff.)
        curve = enneper.Curve(CTRL_PNTS, KNOTS, DEG)
        np.testing.assert_equal(curve.evaluate_at(1), [3.5, 3., 2.5])

    def test_export(self):

        # construct test data
        flo = cStringIO.StringIO()
        curve = enneper.Curve(CTRL_PNTS, KNOTS, DEG)
        curve.export(flo)
        flo.seek(0)
        data = json.load(flo)
        flo.close()

        # test
        np.testing.assert_equal(data['ctrl_pnts'], CTRL_PNTS)
        np.testing.assert_equal(data['knots'], KNOTS)
        self.assertEqual(data['deg'], DEG)


if __name__ == '__main__':
    unittest.main()
