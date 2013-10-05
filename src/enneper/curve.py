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


from lxml import etree
import json
import numpy as np

import foundation as fdn
import cfoundation as cfdn


__all__ = ['Curve']


class Curve(object):

    def __init__(self, ctrl_pnts, knots, deg):
        """designated initializer"""

        fdn.check_nurbs_condition(len(ctrl_pnts), len(knots), deg)

        self._ctrl_pnts = np.asarray(ctrl_pnts, dtype=np.double)
        self._knots = np.asarray(knots, dtype=np.double)
        self._deg = deg

###############################################################################
# constructors
###############################################################################

    @classmethod
    def from_curve(cls, curve):
        return cls(curve.ctrl_pnts.copy(), curve.knots.copy(), curve.deg)

    @classmethod
    def from_json(cls, filename):

        if filename[-5:] != '.json':
            filename = ''.join((filename, '.json'))

        with open(filename, 'r') as open_file:
            data = json.load(open_file)

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

        knots = self.knots
        deg = self.deg

        index = cfdn.get_index(u, deg, knots)
        basis_funs = np.empty((deg + 1, 1), dtype=np.double)
        cfdn.calc_basis_funs(index, u, deg, knots, basis_funs[:, 0])

        lb, ub = index - deg, index + 1
        return np.sum(self.ctrl_pnts[lb:ub] * basis_funs, 0)

###############################################################################
# i/o methods
###############################################################################

    def export_json(self, filename, verbose=False):

        if filename[-5:] != '.json':
            filename = ''.join((filename, '.json'))

        ctrl_pnts = self.ctrl_pnts.tolist()
        knots = self.knots.tolist()

        data = {'ctrl_pnts': ctrl_pnts, 'knots': knots, 'deg':self.deg}

        with open(filename, 'w') as open_file:
            indent = 4 if verbose else None
            json.dump(data, open_file, indent=indent)

    def import_json(self, filename):

        if filename[-5:] != '.json':
            filename = ''.join((filename, '.json'))

        with open(filename, 'r') as open_file:
            data = json.load(open_file)

        ctrl_pnts = data['ctrl_pnts']
        knots = data['knots']
        deg = data['deg']

        fdn.check_nurbs_condition(len(ctrl_pnts), len(knots), deg)

        self._ctrl_pnts = np.asarray(ctrl_pnts, dtype=np.double)
        self._knots = np.asarray(knots, dtype=np.double)
        self._deg = deg
