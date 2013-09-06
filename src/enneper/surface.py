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


class Surface(object):

    def __init__(self):
        """designated initializer"""
        self._ctrl_pnts = None
        self._knots_u = None
        self._knots_u = None
        self._deg_u = None
        self._deg_v = None

###############################################################################
# constructors
###############################################################################

    @classmethod
    def from_surface(cls, surface):
        n_u, n_v, dim = surface.ctrl_pnts.shape
        copied_surface = cls()
        copied_surface.resize(n_u, n_v, dim, surface.deg_u, surface.deg_v)
        copied_surface.ctrl_pnts[:] = surface.ctrl_pnts
        copied_surface.knots_u[:] = surface.knots_u
        copied_surface.knots_v[:] = surface.knots_v
        return copied_surface

    @classmethod
    def from_parameters(cls, ctrl_pnts, knots_u, knots_v, deg_u, deg_v):
        surface = cls()
        n_u, n_v, dim = np.asarray(ctrl_pnts).shape
        surface.resize(n_u, n_v, dim, deg_u, deg_v)
        surface.ctrl_pnts[:] = ctrl_pnts
        surface.knots_u[:] = knots_u
        surface.knots_v[:] = knots_v
        return surface

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
        deg_u, deg_v = self._deg_u, self._deg_v
        knots_u, knots_v = self._knots_u, self._knots_v
        span_u = fdn.find_span(u, deg_u, knots_u)
        basis_funcs_u = fdn.get_basis_funcs(span_u, u, deg_u, knots_u)[:, None]
        span_v = fdn.find_span(v, deg_v, knots_v)
        basis_funcs_v = fdn.get_basis_funcs(span_v, v, deg_v, knots_v)[:, None]
        lb_u, ub_u = span_u - self.deg_u, span_u + 1
        lb_v, ub_v = span_v - self.deg_v,  span_v + 1
        tmp = np.sum(self.ctrl_pnts[lb_u:ub_u, lb_v:ub_v] * basis_funcs_u, 0)
        return np.sum(tmp * basis_funcs_v, 0)
