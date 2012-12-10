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
        self.kwargs = dict()
        ctrl_points = [[0, 0, 1], [4, 4, 4], [3, 2, 1], [4, 1, 1], [5, -1, 1]]
        self.kwargs['ctrl_points'] = ctrl_points
        self.kwargs['knots'] = [0, 0, 0, 1, 2, 3, 3, 3]
        self.kwargs['degree'] = 2

    def test_initializer_not_available(self):
        self.assertRaises(TypeError, ec.NURBSCurve, dummy=2)

    def test_empty_initializer(self):
        ec.NURBSCurve()

    def test_copy_initializer(self):
        copied_curve = ec.NURBSCurve(**self.kwargs)
        curve = ec.NURBSCurve(curve=copied_curve)
        #---- testing ctrl_points assignment ----------------------------------
        actual = curve.ctrl_points
        desired = copied_curve.ctrl_points
        np.testing.assert_equal(actual, desired)
        self.assertNotEqual(id(actual), id(desired))
        #---- testing knots assignment ----------------------------------------
        actual = curve.knots
        desired = copied_curve.knots
        np.testing.assert_equal(actual, desired)
        self.assertNotEqual(id(actual), id(desired))
        #---- testing degree assignment ---------------------------------------
        actual = curve.degree
        desired = copied_curve.degree
        np.testing.assert_equal(actual, desired)

    def test_designated_initializer(self):
        curve = ec.NURBSCurve(**self.kwargs)
        #---- testing ctrl_points assignment ----------------------------------
        actual = curve.ctrl_points
        desired = self.kwargs['ctrl_points']
        np.testing.assert_equal(actual, desired)
        #---- testing knots assignment ----------------------------------------
        actual = curve.knots
        desired = self.kwargs['knots']
        np.testing.assert_equal(actual, desired)
        #---- testing degree assignment ---------------------------------------
        actual = curve.degree
        desired = self.kwargs['degree']
        np.testing.assert_equal(actual, desired)
        #---- testing n + p + 1 != m ------------------------------------------
        desired = self.kwargs['degree'] = 3
        self.assertRaises(ValueError, ec.NURBSCurve, **self.kwargs)
    
    def test_resize(self):
        curve = ec.NURBSCurve(shape=(5, 2), degree=2)
        np.testing.assert_equal(curve.ctrl_points.shape, (5, 2))
        np.testing.assert_equal(curve.degree, 2)
        np.testing.assert_equal(curve.knots.shape, (8,))

    def test_evaluate_at(self):
        """The NURBS Book 2nd edition: example 4.1 (page 122 ff.)"""
        curve = ec.NURBSCurve(**self.kwargs)
        actual = curve.evaluate_at(1)
        desired = np.asarray([1.4, 1.2])
        np.testing.assert_equal(actual, desired)
        actual = curve.evaluate_at(1, homogeneous=True)
        desired = np.asarray([3.5, 3, 2.5])
        np.testing.assert_equal(actual, desired)
