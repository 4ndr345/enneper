#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ***************************************************************************
# *   Copyright (C) 2012 by Andreas Kührmann [andreas.kuehrmann@gmail.com]  *
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

import enneper.surface as enneper


DIR = os.path.join(os.path.dirname(__file__), 'sphere')

CTRL_PNTS = np.fromfile(os.path.join(DIR, 'ctrl_pnts')).reshape(-1, 5, 4)
KNOTS_U = np.fromfile(os.path.join(DIR, 'knots_u'))
KNOTS_V = np.fromfile(os.path.join(DIR, 'knots_v'))
DEG_U = 2
DEG_V = 2


class TestCurve(unittest.TestCase):

    def test_designated_initializer(self):
        surface = enneper.Surface(CTRL_PNTS, KNOTS_U, KNOTS_V, DEG_U, DEG_V)
        np.testing.assert_array_equal(surface.ctrl_pnts, CTRL_PNTS)
        np.testing.assert_array_equal(surface.knots_u, KNOTS_U)
        np.testing.assert_array_equal(surface.knots_v, KNOTS_V)
        self.assertEqual(surface.deg_u, DEG_U)
        self.assertEqual(surface.deg_v, DEG_V)

    def test_from_surface_constructor(self):
        original = enneper.Surface(CTRL_PNTS, KNOTS_U, KNOTS_V, DEG_U, DEG_V)
        copy = enneper.Surface.from_surface(original)
        np.testing.assert_array_equal(original.ctrl_pnts, copy.ctrl_pnts)
        np.testing.assert_array_equal(original.knots_u, copy.knots_u)
        np.testing.assert_array_equal(original.knots_v, copy.knots_v)
        self.assertEqual(original.deg_u, copy.deg_u)
        self.assertEqual(original.deg_v, copy.deg_v)
        self.assertIsNot(original.ctrl_pnts, copy.ctrl_pnts)
        self.assertIsNot(original.knots_u, copy.knots_u)
        self.assertIsNot(original.knots_v, copy.knots_v)

    def test_evaluate_at(self):
        surface = enneper.Surface(CTRL_PNTS, KNOTS_U, KNOTS_V, DEG_U, DEG_V)
        np.testing.assert_array_equal(surface.evaluate_at(0, 0), [ 0.,  0.,  1.,  1.])


if __name__ == '__main__':
    unittest.main()
