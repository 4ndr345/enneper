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


# standard packages
import json

# 3rd party packages
import numpy as np

# project packages
import foundation as fdn
import cfoundation as cfdn


__all__ = ['Curve']


class Curve(object):

    def __init__(self, ctrl_pnts, knots):
        self.ctrl_pnts = np.asarray(ctrl_pnts, dtype=np.double)
        self.knots = np.asarray(knots, dtype=np.double)

###############################################################################
# constructors
###############################################################################

    @classmethod
    def from_curve(cls, curve):
        return cls(curve.ctrl_pnts.copy(), curve.knots.copy())

    @classmethod
    def from_json(cls, flo):
        data = json.load(flo)
        return cls(data['ctrl_pnts'], data['knots'])

###############################################################################
# properties
###############################################################################

    @property
    def deg(self):
        return len(self.knots) - len(self.ctrl_pnts) - 1

###############################################################################
# miscellaneous methods
###############################################################################

    def evaluate_at(self, u):

        # save local to avoid looking up twice or more
        knots = self.knots
        deg = self.deg

        # low level nurbs function written in cython
        index = cfdn.get_index(u, deg, knots)
        basis_funs = np.empty((deg + 1, 1), dtype=np.double)
        cfdn.calc_basis_funs(index, u, deg, knots, basis_funs[:, 0])

        # calc homogeneous point
        lb, ub = index - deg, index + 1
        return np.sum(self.ctrl_pnts[lb:ub] * basis_funs, 0)

    def transform(self, matrix):
        self.ctrl_pnts = np.dot(matrix, self.ctrl_pnts)

    def export(self, flo, indent=None):

        # json can't handle ndarrays
        data = dict(ctrl_pnts=self.ctrl_pnts.tolist(), knots=self.knots.tolist())
        json.dump(data, flo, indent=indent)
