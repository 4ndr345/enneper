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


__all__ = ['Surface']


class Surface(object):

    def __init__(self, ctrl_pnts, knots_u, knots_v, deg_u, deg_v):
        """designated initializer"""

        if len(ctrl_pnts) + deg_u + 1 != len(knots_u):
            raise fdn.NURBSException(len(ctrl_pnts), len(knots_u), deg_u)

        if len(ctrl_pnts[0]) + deg_v + 1 != len(knots_v):
            raise fdn.NURBSException(len(ctrl_pnts[0]), len(knots_v), deg_v)

        self._ctrl_pnts = np.asarray(ctrl_pnts, dtype=np.double)
        self._knots_u = np.asarray(knots_u, dtype=np.double)
        self._knots_v = np.asarray(knots_v, dtype=np.double)
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

    def evaluate_at(self, u, v):
        knots_u, knots_v = self._knots_u, self._knots_v
        deg_u, deg_v = self._deg_u, self._deg_v

        index_u = cfdn.get_index(u, deg_u, knots_u)
        basis_funs_u = np.empty((deg_u + 1, 1, 1), dtype=np.double)
        cfdn.calc_basis_funs(index_u, u, deg_u, knots_u, basis_funs_u[:, 0, 0])

        index_v = cfdn.get_index(v, deg_v, knots_v)
        basis_funs_v = np.empty((deg_v + 1, 1), dtype=np.double)
        cfdn.calc_basis_funs(index_v, v, deg_v, knots_v, basis_funs_v[:, 0])

        lb_u, ub_u = index_u - deg_u, index_u + 1
        lb_v, ub_v = index_v - deg_v, index_v + 1
        tmp = np.sum(self.ctrl_pnts[lb_u:ub_u, lb_v:ub_v] * basis_funs_u, 0)
        return np.sum(tmp * basis_funs_v, 0)
