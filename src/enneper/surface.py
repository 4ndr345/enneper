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


import numpy as np

import foundation as fdn
import cfoundation as cfdn


class Surface(object):

    def __init__(self, ctrl_pnts, knots_u, knots_v, deg_u, deg_v):
        """designated initializer"""
        if len(ctrl_pnts) + deg_u + 1 != len(knots_u):
            raise fdn.NURBSException(len(ctrl_pnts), len(knots_u), deg_u)
        if len(ctrl_pnts[0]) + deg_v + 1 != len(knots_v):
            raise fdn.NURBSException(len(ctrl_pnts[0]), len(knots_v), deg_v)
        self._ctrl_pnts = np.asarray(ctrl_pnts)
        self._knots_u = np.asarray(knots_u)
        self._knots_v = np.asarray(knots_v)
        self._deg_u = deg_u
        self._deg_v = deg_v

###############################################################################
# constructors
###############################################################################

    @classmethod
    def from_surface(cls, surface):
        ctrl_pnts = surface.ctrl_pnts.copy()
        knots_u = surface.knots_u.copy()
        knots_v = surface.knots_v.copy()
        deg_u = surface.deg_u
        deg_v = surface.deg_v
        return cls(ctrl_pnts, knots_u, knots_v, deg_u, deg_v)

###############################################################################
# properties
###############################################################################

    @property
    def ctrl_pnts(self):
        return self._ctrl_pnts

    @property
    def knots_u(self):
        return self._knots_u

    @property
    def knots_v(self):
        return self._knots_v

    @property
    def deg_u(self):
        return self._deg_u

    @property
    def deg_v(self):
        return self._deg_v

###############################################################################
# miscellaneous methods
###############################################################################

    def resize(self, n_u, n_v, dim, deg_u, deg_v):
        self._ctrl_pnts = np.zeros((n_u, n_v, dim))
        self._knots_u = np.zeros(n_u + deg_u + 1)
        self._knots_v = np.zeros(n_v + deg_v + 1)
        self._deg_u = deg_u
        self._deg_v = deg_v

    def evaluate_at(self, u, v):
        ctrl_pnts = self.ctrl_pnts
        knots_u, knots_v = self._knots_u, self._knots_v
        deg_u, deg_v = self._deg_u, self._deg_v
        span_u = cfdn.find_span(ctrl_pnts.shape[0], deg_u, u, knots_u)
        basis_funs_u = np.empty((deg_u + 1, 1, 1), dtype=np.double)
        cfdn.basis_funs(span_u, u, deg_u, knots_u, basis_funs_u[:, 0, 0])
        span_v = cfdn.find_span(ctrl_pnts.shape[1], deg_v, v, knots_v)
        basis_funs_v = np.empty((deg_v + 1, 1), dtype=np.double)
        cfdn.basis_funs(span_v, v, deg_v, knots_v, basis_funs_v[:, 0])
        lb_u, ub_u = span_u - deg_u, span_u + 1
        lb_v, ub_v = span_v - deg_v, span_v + 1
        tmp = np.sum(ctrl_pnts[lb_u:ub_u, lb_v:ub_v] * basis_funs_u, 0)
        return np.sum(tmp * basis_funs_v, 0)
