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


import enneper.curves as ec
import numpy as np


class NURBSCurve(unittest.TestCase):

    def setUp(self):
        self._args = [
            np.array([[0, 0, 1], [4, 4, 4], [3, 2, 1], [4, 1, 1], [5, -1, 1]]),
            np.array([0, 0, 0, 1, 2, 3, 3, 3]),
            2
        ]
        self._curve = ec.NURBSCurve(*self._args)

    def test_init(self):
        #---- testing ctrl_points assignment ----------------------------------
        actual = self._curve.ctrl_points
        desired = self._args[0]
        np.testing.assert_equal(actual, desired)
        #---- testing knots assignment ----------------------------------------
        actual = self._curve.knots
        desired = self._args[1]
        np.testing.assert_equal(actual, desired)
        #---- testing degree assignment ---------------------------------------
        actual = self._curve.degree
        desired = self._args[2]
        np.testing.assert_equal(actual, desired)
        #---- testing n + p + 1 != m ------------------------------------------
        self._args[2] = 3
        self.assertRaises(ValueError, ec.NURBSCurve, *self._args)

    def test_initializer_not_available(self):
        self.assertRaises(TypeError, ec.NURBSCurve)

    def test_evaluate_at(self):
        """The NURBS Book 2nd edition: example 4.1 (page 122 ff.)"""
        actual = self._curve.evaluate_at(1)
        desired = np.array([1.4, 1.2])
        np.testing.assert_equal(actual, desired)
        actual = self._curve.evaluate_at(1, True)
        desired = np.array([3.5, 3, 2.5])
        np.testing.assert_equal(actual, desired)
