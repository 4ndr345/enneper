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

    def __init__(self, ctrl_pnts, knots, deg):

        # raise a NURBSError when len(ctrl_pnts) + deg + 1 != len(knots)
        fdn.check_nurbs_condition(len(ctrl_pnts), len(knots), deg)
        
        # assign attributes
        self._ctrl_pnts = np.asarray(ctrl_pnts, dtype=np.double)
        self._knots = np.asarray(knots, dtype=np.double)
        self._deg = deg

###############################################################################
# constructors
###############################################################################

    @classmethod
    def from_curve(cls, curve):

        # nothig special copy and call designated initializer
        return cls(curve.ctrl_pnts.copy(), curve.knots.copy(), curve.deg)

    @classmethod
    def from_json(cls, flo):

        # nothing special load data and call designated initializer
        data = json.load(flo)
        return cls(data['ctrl_pnts'], data['knots'], data['deg'])

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

    def export(self, flo, indent=None):

        # json can't handle ndarrays
        ctrl_pnts = self.ctrl_pnts.tolist()
        knots = self.knots.tolist()

        # export curve
        data = dict(ctrl_pnts=ctrl_pnts, knots=knots, deg=self.deg)
        json.dump(data, flo, indent=indent)
