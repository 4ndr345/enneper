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


import numpy as np

import enneper.foundation.nurbs as efn


__all__ = ['NURBSCurve']


class NURBSCurve(object):

    _init_args = {
        ('NURBSCurve',): '_copy',
        ('ndarray', 'ndarray', 'int'): '_init',
    }

    def __init__(self, *args):
        key = tuple(type(arg).__name__ for arg in args)
        if key not in NURBSCurve._init_args:
            raise TypeError
        self._ctrl_points = None
        self._degree = None
        self._knots = None
        getattr(self, NURBSCurve._init_args[key])(*args)

    def _init(self, ctrl_points, knots, degree):
        if ctrl_points.shape[0] + degree + 1 != knots.size:
            raise ValueError
        self._ctrl_points = ctrl_points
        self._degree = degree
        self._knots = knots

    def _copy(self, curve):
        self._ctrl_points = curve.ctrl_points.copy()
        self._degree = curve.degree
        self._knots = curve.knots.copy()

    @ property
    def ctrl_points(self):
        return self._ctrl_points

    @ property
    def degree(self):
        return self._degree

    @ property
    def knots(self):
        return self._knots

    def evaluate_at(self, t, homogenous=False):
        span = efn.find_span(t, self.degree, self.knots)
        n = efn.get_basis_functions(span, t, self.degree, self.knots)
        point = np.zeros(self.ctrl_points.shape[1])
        # todo: optimize loop: use np.sum and np.multiply
        for i in range(self.degree + 1):
            point += n[i] * self.ctrl_points[span - self.degree + i]
        if homogenous:
            return point
        point /= point[-1]
        return point[0:-1]
