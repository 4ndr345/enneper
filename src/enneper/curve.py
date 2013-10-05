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


__all__ = ['Curve']


class Curve(object):

    def __init__(self, ctrl_pnts, knots, deg):
        """designated initializer"""

        if len(ctrl_pnts) + deg + 1 != len(knots):
            raise fdn.NURBSException(len(ctrl_pnts), len(knots), deg)

        self._ctrl_pnts = np.asarray(ctrl_pnts, dtype=np.double)
        self._knots = np.asarray(knots, dtype=np.double)
        self._deg = deg

###############################################################################
# constructors
###############################################################################

    @classmethod
    def from_curve(cls, curve):
        ctrl_pnts = curve.ctrl_pnts.copy()
        knots = curve.knots.copy()
        deg = curve.deg
        return cls(ctrl_pnts, knots, deg)

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
        knots = self.knots
        deg = self.deg

        index = cfdn.get_index(u, deg, knots)
        basis_funs = np.empty((deg + 1, 1), dtype=np.double)
        cfdn.calc_basis_funs(index, u, deg, knots, basis_funs[:, 0])

        lb, ub = index - deg, index + 1
        return np.sum(self.ctrl_pnts[lb:ub] * basis_funs, 0)
