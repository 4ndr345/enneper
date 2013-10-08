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


__all__ = ['Surface']


class Surface(object):

    def __init__(self, ctrl_pnts, knots_u, knots_v):
        self.ctrl_pnts = np.asarray(ctrl_pnts, dtype=np.double)
        self.knots_u = np.asarray(knots_u, dtype=np.double)
        self.knots_v = np.asarray(knots_v, dtype=np.double)

###############################################################################
# constructors
###############################################################################

    @classmethod
    def from_surface(cls, surface):

        # copy data
        ctrl_pnts = surface.ctrl_pnts.copy()
        knots_u = surface.knots_u.copy()
        knots_v = surface.knots_v.copy()

        # call designated initializer
        return cls(ctrl_pnts, knots_u, knots_v)

    @classmethod
    def from_json(cls, file_like_obj):
        return cls(**json.load(file_like_obj))

###############################################################################
# properties
###############################################################################

    @property
    def deg_u(self):
        return self.knots_u.size - self.ctrl_pnts.shape[0] - 1

    @property
    def deg_v(self):
        return self.knots_v.size - self.ctrl_pnts.shape[1] - 1

###############################################################################
# miscellaneous methods
###############################################################################

    def evaluate_at(self, u, v):

        # save local to avoid looking up twice or more
        knots_u, knots_v = self.knots_u, self.knots_v
        deg_u, deg_v = self.deg_u, self.deg_v

        # low level nurbs function written in cython
        index_u = cfdn.get_index(u, deg_u, knots_u)
        basis_funs_u = np.empty((deg_u + 1, 1, 1), dtype=np.double)
        cfdn.calc_basis_funs(index_u, u, deg_u, knots_u, basis_funs_u[:, 0, 0])

        # low level nurbs function written in cython
        index_v = cfdn.get_index(v, deg_v, knots_v)
        basis_funs_v = np.empty((deg_v + 1, 1), dtype=np.double)
        cfdn.calc_basis_funs(index_v, v, deg_v, knots_v, basis_funs_v[:, 0])

        # calc homogeneous point
        lb_u, ub_u = index_u - deg_u, index_u + 1
        lb_v, ub_v = index_v - deg_v, index_v + 1
        tmp = np.sum(self.ctrl_pnts[lb_u:ub_u, lb_v:ub_v] * basis_funs_u, 0)
        return np.sum(tmp * basis_funs_v, 0)

    def export(self, file_like_obj, indent=None):

        # json can't handle ndarrays
        obj_to_serialize = dict()
        obj_to_serialize['ctrl_pnts'] = self.ctrl_pnts.tolist()
        obj_to_serialize['knots_u'] =  self.knots_u.tolist()
        obj_to_serialize['knots_v'] =  self.knots_v.tolist()

        # export
        json.dump(obj_to_serialize, file_like_obj, indent=indent)
