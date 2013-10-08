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


import unittest

import numpy as np

import enneper


class TestPrimitives(unittest.TestCase):

    def test_line(self):
        line = enneper.primitives.Line([0, 0, 0], [1, 1, 1])
        np.testing.assert_equal(line.evaluate_at(0), [0, 0, 0, 1])
        np.testing.assert_equal(line.evaluate_at(1), [1, 1, 1, 1])

    def test_circular_arc (self):
        circle = enneper.primitives.CircularArc(2 * np.pi)
        np.testing.assert_almost_equal(circle.evaluate_at(0), [1, 0, 1])
        np.testing.assert_almost_equal(circle.evaluate_at(.25), [0, 1, 1])
        np.testing.assert_almost_equal(circle.evaluate_at(.5), [-1, 0, 1])
        np.testing.assert_almost_equal(circle.evaluate_at(.75), [0, -1, 1])
        np.testing.assert_almost_equal(circle.evaluate_at(1), [1, 0, 1])


if __name__ == '__main__':
    unittest.main()
