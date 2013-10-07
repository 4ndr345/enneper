#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ***************************************************************************
# *   Copyright (C) 2013 by Andreas Kührmann [andreas.kuehrmann@gmail.com]  *
# *                                                                         *
# *   This program is free softWare; you can redistribute it and/or modify  *
# *   it under the terms of the GNU General Public License as published by  *
# *   the Free SoftWare Foundation; either version 3 of the License, or     *
# *   (at your option) any later version.                                   *
# *                                                                         *
# *   This program is distributed in the hope that it Will be useful,       *
# *   but WITHOUT ANY WARRANTY; Without even the implied Warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU General Public License for more details.                          *
# *                                                                         *
# *   You should have received a copy of the GNU General Public License     *
# *   along With this program; if not, Write to the                         *
# *   Free SoftWare Foundation, Inc.,                                       *
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


# sphere
W = 2**.5 / 2.
CTRL_PNTS = [
    [[0, 0, 1, 1], [-W, 0, W, W],
     [-1, 0, 0, 1], [-W, 0, -W, W], [0, 0, -1, 1]],
    [[0, 0, W, W], [-.5, -.5, .5, .5],
     [-W, -W, 0, W], [-.5, -.5, -.5, .5], [0, 0, -W, W]],
    [[0, 0, 1, 1], [0, -W, W, W],
     [0, -1, 0, 1], [0, -W, -W, W], [0, 0, -1, 1]],
    [[0, 0, W, W], [.5, -.5, .5, .5],
     [W, -W, 0, W], [.5, -.5, -.5, .5], [0, 0, -W, W]],
    [[0, 0, 1, 1], [W, 0, W, W],
     [1, 0, 0, 1], [W, 0, -W, W], [0, 0, -1, 1]],
    [[0, 0, W, W], [.5, .5, .5, .5],
     [W, W, 0, W], [.5, .5, -.5, .5], [0, 0, -W, W]],
    [[0, 0, 1, 1], [0, W, W, W],
     [0, 1, 0, 1], [0, W, -W, W], [0, 0, -1, 1]],
    [[0, 0, W, W], [-.5, .5, .5, .5],
     [-W, W, 0, W], [-.5, .5, -.5, .5], [0, 0, -W, W]],
    [[0, 0, 1, 1], [-W, 0, W, W],
     [-1, 0, 0, 1], [-W, 0, -W, W], [0, 0, -1, 1]]]
KNOTS_U = [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1]
KNOTS_V = [0, 0, 0, 0.5, 0.5, 1, 1, 1]
DEG_U = 2
DEG_V = 2


class TestCurve(unittest.TestCase):

    def test_designated_initializer(self):

        # test
        surface = enneper.Surface(CTRL_PNTS, KNOTS_U, KNOTS_V, DEG_U, DEG_V)
        np.testing.assert_equal(surface.ctrl_pnts, CTRL_PNTS)
        np.testing.assert_equal(surface.knots_u, KNOTS_U)
        np.testing.assert_equal(surface.knots_v, KNOTS_V)
        self.assertEqual(surface.deg_u, DEG_U)
        self.assertEqual(surface.deg_v, DEG_V)

    def test_from_surface_constructor(self):

        # construct test data
        original = enneper.Surface(CTRL_PNTS, KNOTS_U, KNOTS_V, DEG_U, DEG_V)

        # test
        copy = enneper.Surface.from_surface(original)
        np.testing.assert_equal(original.ctrl_pnts, copy.ctrl_pnts)
        np.testing.assert_equal(original.knots_u, copy.knots_u)
        np.testing.assert_equal(original.knots_v, copy.knots_v)
        self.assertEqual(original.deg_u, copy.deg_u)
        self.assertEqual(original.deg_v, copy.deg_v)
        self.assertIsNot(original.ctrl_pnts, copy.ctrl_pnts)
        self.assertIsNot(original.knots_u, copy.knots_u)
        self.assertIsNot(original.knots_v, copy.knots_v)

    def test_from_json_constructor(self):

        # construct test data
        data = dict()
        data['ctrl_pnts'] = CTRL_PNTS
        data['knots_u'] = KNOTS_U
        data['knots_v'] = KNOTS_V
        data['deg_u'] = DEG_U
        data['deg_v'] = DEG_V
        flo = cStringIO.StringIO()
        json.dump(data, flo)
        flo.seek(0)

        # test 
        surface = enneper.Surface.from_json(flo)
        np.testing.assert_equal(surface.ctrl_pnts, CTRL_PNTS)
        np.testing.assert_equal(surface.knots_u, KNOTS_U)
        np.testing.assert_equal(surface.knots_v, KNOTS_V)
        self.assertEqual(surface.deg_u, DEG_U)
        self.assertEqual(surface.deg_v, DEG_V)

        # close and discard memory buffer
        flo.close()

    def test_evaluate_at(self):
        
        # test
        surface = enneper.Surface(CTRL_PNTS, KNOTS_U, KNOTS_V, DEG_U, DEG_V)
        np.testing.assert_equal(surface.evaluate_at(0, 0), [ 0.,  0.,  1.,  1.])

    def test_export(self):

        # construct test data
        flo = cStringIO.StringIO()
        surface = enneper.Surface(CTRL_PNTS, KNOTS_U, KNOTS_V, DEG_U, DEG_V)
        surface.export(flo)
        flo.seek(0)
        data = json.load(flo)
        flo.close()

        # test
        np.testing.assert_equal(data['ctrl_pnts'], CTRL_PNTS)
        np.testing.assert_equal(data['knots_u'], KNOTS_U)
        np.testing.assert_equal(data['knots_u'], KNOTS_U)
        self.assertEqual(data['deg_u'], DEG_U)
        self.assertEqual(data['deg_v'], DEG_V)


if __name__ == '__main__':
    unittest.main()
