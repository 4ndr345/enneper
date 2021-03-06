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

import enneper.foundation as fdn


class TestFoundation(unittest.TestCase):

    def setUp(self):
        self.degree = 2
        self.knots = np.array([0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5])
        self.weights = np.array([1, .7, 1, .7, 1, .7, 1, .7, 1])
        self.pnts = np.array([[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0],
                              [-1, -1], [0, -1], [1, -1], [1, 0]])

    def test_find_span(self):
        """The NURBS Book 2nd edition: example 2.4 (page 71)"""
        actual = fdn.nurbs.find_span(2.5, self.degree, self.knots)
        np.testing.assert_equal(actual, 4)

    def test_get_basis_functions(self):
        """The NURBS Book 2nd edition: example 2.4 (page 71)"""
        actual = fdn.nurbs.get_basis_functions(4, 2.5, self.degree, self.knots)
        np.testing.assert_equal(actual, np.array([0.125, 0.75, 0.125]))

    def test_get_homogeneous_points_of_a_point(self):
        pnts = np.array([1, 1])
        weights = np.array([2])
        actual = fdn.nurbs.get_homogeneous_points(pnts, weights)
        desired = np.array([2, 2, 2])
        np.testing.assert_equal(actual, desired)

    def test_get_homogeneous_points_of_a_curve(self):
        pnts = np.array([[1, 0], [1, 1]])
        weights = np.array([1, 2])
        actual = fdn.nurbs.get_homogeneous_points(pnts, weights)
        desired = np.array([[1, 0, 1], [2, 2, 2]])
        np.testing.assert_equal(actual, desired)

    def test_get_homogeneous_points_of_a_surface(self):
        pnts = np.array([[[1, 0], [1, 1]], [[0, 1], [-1, 1]]])
        weights = np.array([[1, 2], [1, 2]])
        actual = fdn.nurbs.get_homogeneous_points(pnts, weights)
        desired = np.array([[[1, 0, 1], [2, 2, 2]], [[0, 1, 1], [-2, 2, 2]]])
        np.testing.assert_equal(actual, desired)

    def test_get_cartesian_points_of_a_point(self):
        ctrl_pnts = np.array([2, 2, 2])
        actual = fdn.nurbs.get_cartesian_points(ctrl_pnts)
        desired = np.array([1, 1])
        np.testing.assert_equal(actual, desired)

    def test_get_cartesian_points_of_a_curve(self):
        ctrl_pnts = np.array([[1, 0, 1], [2, 2, 2]])
        actual = fdn.nurbs.get_cartesian_points(ctrl_pnts)
        desired = np.array([[1, 0], [1, 1]])
        np.testing.assert_equal(actual, desired)

    def test_get_cartesian_points_of_a_surface(self):
        ctrl_pnts = np.array([[[1, 0, 1], [2, 2, 2]], [[0, 1, 1], [-2, 2, 2]]])
        actual = fdn.nurbs.get_cartesian_points(ctrl_pnts)
        desired = np.array([[[1, 0], [1, 1]], [[0, 1], [-1, 1]]])
        np.testing.assert_equal(actual, desired)
