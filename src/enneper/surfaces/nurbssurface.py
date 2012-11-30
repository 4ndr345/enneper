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


__all__ = ['NURBSSurface']


class NURBSSurface(object):

    _init_args = {
        ('tuple', 'int', 'int'): 'resize',
        ('ndarray', 'ndarray', 'ndarray', 'int', 'int'): '_init',
    }

    def __init__(self, *args):
        key = tuple(type(arg).__name__ for arg in args)
        if key not in NURBSSurface._init_args:
            raise TypeError
        self._ctrl_points = None
        self._degree_u = None
        self._degree_v = None
        self._knots_u = None
        self._knots_v = None
        getattr(self, NURBSSurface._init_args[key])(*args)

    def _init(self, ctrl_points, knots_u, knots_v, degree_u, degree_v):
        if ctrl_points.shape[0] + degree_u + 1 != knots_u.size:
            raise ValueError
        if ctrl_points.shape[1] + degree_v + 1 != knots_v.size:
            raise ValueError
        self._ctrl_points = ctrl_points
        self._degree_u = degree_u
        self._degree_v = degree_v
        self._knots_u = knots_u
        self._knots_v = knots_v

    @ property
    def ctrl_points(self):
        return self._ctrl_points

    @ property
    def degree_u(self):
        return self._degree_u

    @ property
    def degree_v(self):
        return self._degree_v

    @ property
    def knots_u(self):
        return self._knots_u

    @ property
    def knots_v(self):
        return self._knots_v

    def evaluate_at(self, u, v, homogenous=False):
        span_u = efn.find_span(u, self.degree_u, self.knots_u)
        n_u = efn.get_basis_functions(span_u, u, self.degree_u, self.knots_u)
        span_v = efn.find_span(v, self.degree_v, self.knots_v)
        n_v = efn.get_basis_functions(span_v, v, self.degree_v, self.knots_v)
        index_u, index_v = span_u - self.degree_u, span_v - self.degree_v
        tmp = np.zeros((self.degree_v + 1, self.ctrl_points.shape[2]))
        for i in range(self.degree_v + 1):
            for j in range(self.degree_u + 1):
                tmp[i] += n_u[j] * self.ctrl_points[index_u + j, index_v + i]
        point = np.zeros(self.ctrl_points.shape[2])
        # todo: optimize loop: use np.sum and np.multiply
        for i in range(self.degree_v + 1):
            point += n_v[i] * tmp[i]
        if homogenous:
            return point
        point /= point[-1]
        return point[0:-1]

    def resize(self, shape, degree_u, degree_v):
        self._ctrl_points = np.zeros(shape)
        self._degree_u = degree_u
        self._degree_v = degree_v
        self._knots_u = np.zeros(shape[0] + degree_u + 1)
        self._knots_v = np.zeros(shape[1] + degree_v + 1)
