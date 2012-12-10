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

    _init_kwargs = {
        (): '_empty_initializer',
        ('curve',): '_copy_initializer',
        ('degree', 'shape'): 'resize',
        ('ctrl_points', 'degree', 'knots'): '_designated_initializer'
    }

    def __init__(self, **kwargs):
        key = tuple(sorted(kwargs.keys()))
        if key not in NURBSCurve._init_kwargs:
            raise TypeError
        self._ctrl_points = None
        self._degree = None
        self._knots = None
        getattr(self, NURBSCurve._init_kwargs[key])(**kwargs)

    def _empty_initializer(self):
        pass

    def _copy_initializer(self, curve):
        self._ctrl_points = curve.ctrl_points.copy()
        self._degree = curve.degree
        self._knots = curve.knots.copy()

    def _designated_initializer(self, ctrl_points, knots, degree):
        ctrl_points = np.asarray(ctrl_points)
        knots = np.asarray(knots)
        if ctrl_points.shape[0] + degree + 1 != knots.size:
            raise ValueError
        self._ctrl_points = ctrl_points
        self._degree = degree
        self._knots = knots

    @ property
    def ctrl_points(self):
        return self._ctrl_points

    @ property
    def degree(self):
        return self._degree

    @ property
    def knots(self):
        return self._knots

    def elevate_degree(self, t):
        raise NotImplementedError

    def evaluate_at(self, u, **opt):
        homogeneous = opt.get('homogeneous', False)
        span = efn.find_span(u, self.degree, self.knots)
        n = efn.get_basis_functions(span, u, self.degree, self.knots)
        point = np.zeros(self.ctrl_points.shape[1])
        # todo: optimize loop: use np.sum and np.multiply
        for i in range(self.degree + 1):
            point += n[i] * self.ctrl_points[span - self.degree + i]
        if homogeneous:
            return point
        point /= point[-1]
        return point[0:-1]

    def merge_knots(self, knots):
        knots = np.asarray(knots)
        # Find the knots to insert
        done, i, ia, ib = False, 0, 0, 0
        knots_to_insert = np.zeros(knots.shape)
        while not done:
            if knots[ib] == self.knots[ia]:
                ib, ia = ib + 1, ia + 1
            else:
                knots_to_insert[i], i, ib = knots[ib], i + 1, ib + 1
            done = ia >= self.knots.size or ib >= knots.size
        # Refine the curve if necessary
        if knots_to_insert.size > 0:
            self.refine_knots(knots_to_insert[:i])

    def refine_knots(self, knots):
        raise NotImplementedError

    def resize(self, shape, degree):
        self._ctrl_points = np.zeros(shape)
        self._degree = degree
        self._knots = np.zeros(shape[0] + degree + 1)
