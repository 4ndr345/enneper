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
        m = self.knots.shape[0] - 1
        n = self.ctrl_points.shape[0] - 1
        p = self.degree
        r = knots.shape[0] - 1
        curve = NURBSCurve(curve=self)
        self.resize((r + n + 2, self.ctrl_points.shape[1]), p)
        a = efn.find_span(knots[0], curve.degree, curve.knots)
        b = efn.find_span(knots[r], curve.degree, curve.knots) + 1
        l, u = 0, a - p + 1
        self.ctrl_points[l:u] = curve.ctrl_points[l:u]
        l, u = b - 1, n + 1
        self.ctrl_points[(l + r + 1):(u + r + 1)] = curve.ctrl_points[l:u]
        l, u = 0, a + 1
        self.knots[l:u] = curve.knots[l:u]
        l, u = b + p, m + 1
        self.knots[(l + r + 1):(u + r + 1)] = curve.knots[l:u]
        i, k = b + p - 1, b + p + r
        for j in range(r, -1, -1):
            while knots[j] <= curve.knots[i] and i > a:
                self.ctrl_points[k - p - 1] = curve.ctrl_points[i - p - 1]
                self.knots[k] = curve.knots[i]
                k, i = k - 1, i - 1
            self.ctrl_points[k - p - 1] = self.ctrl_points[k - p]
            for l in range(1, p + 1):
                pre = k - p
                ind = pre + 1
                alpha = self.knots[k + l] - knots[j]
                if abs(alpha) < 1e-7:
                    self.ctrl_points[pre] = self.ctrl_points[ind]
                else:
                    alpha /= self.knots[k + l] - curve.knots[i - p + l]
                self.ctrl_points[pre] *= alpha
                self.ctrl_points[pre] += (1 - alpha) * self.ctrl_points[ind]
            self.knots[k] = knots[j]
            k -= 1

    def resize(self, shape, degree):
        self._ctrl_points = np.zeros(shape)
        self._degree = degree
        self._knots = np.zeros(shape[0] + degree + 1)
