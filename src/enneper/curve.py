#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ***************************************************************************
# *   Copyright (C) 2012 by Andreas Kührmann [andreas.kuehrmann@gmail.com]  *
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
import foundation_cython as cfdn



__all__ = ['Curve']


class Curve(object):

    def __init__(self, ctrl_pnts, knots, deg):
        """designated initializer"""
        if len(ctrl_pnts) + deg + 1 != len(knots):
            raise fdn.NURBSException(len(ctrl_pnts), len(knots), deg)
        self._ctrl_pnts = np.asarray(ctrl_pnts)
        self._knots = np.asarray(knots)
        self._deg = deg

###############################################################################
# constructors
###############################################################################

    @classmethod
    def from_curve(cls, curve):
        copied_curve = cls()
        n, dim = curve.ctrl_pnts.shape
        copied_curve.resize(n, dim, curve.deg)
        copied_curve.ctrl_pnts[:] = curve.ctrl_pnts
        copied_curve.knots[:] = curve.knots
        return copied_curve

###############################################################################
# properties
###############################################################################

    @property
    def ctrl_pnts(self):
        return self._ctrl_pnts

    @property
    def knots(self):
        return self._knots

    @property
    def deg(self):
        return self._deg

###############################################################################
# miscellaneous methods
###############################################################################

    def evaluate_at(self, u):
        deg, knots, ctrl_pnts = self.deg, self.knots, self.ctrl_pnts
        span = cfdn.find_span(ctrl_pnts.shape[0], deg, u, knots)
        basis_funs = np.empty((deg + 1, 1), dtype=np.double)
        cfdn.basis_funs(span, u, deg, knots, basis_funs[:, 0])
        lb, ub = span - deg, span + 1
        return np.sum(ctrl_pnts[lb:ub] * basis_funs, 0)

    def resize(self, n, dim, deg):
        self._ctrl_pnts = np.zeros((n, dim))
        self._knots = np.zeros(n + deg + 1)
        self._deg = deg
