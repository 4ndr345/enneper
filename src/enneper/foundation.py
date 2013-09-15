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


def find_span(u, deg, knots):
    """Determine the knot span index"""
    lb, ub = deg, knots.size - deg
    # Special case
    if u == knots[ub - 1]:
        return ub - 2
    # search
    return np.where(knots[lb:ub] <= u)[0][-1] + deg


def get_basis_funcs(i, u, deg, knots):
    """ Computes the nonvanishing basis functions"""
    basis_funcs = np.ones(deg + 1)
    left = np.ones(deg + 1)
    right = np.ones(deg + 1)
    for j in range(1, deg + 1):
        left[j] = u - knots[i + 1 - j]
        right[j] = knots[i + j] - u
        cache = 0
        for k in range(j):
            tmp = basis_funcs[k] / (right[k + 1] + left[j - k])
            basis_funcs[k] = cache + right[k + 1] * tmp
            cache = left[j - k] * tmp
        basis_funcs[j] = cache
    return basis_funcs


def get_cartesian_points(h_pnts):
    shape = list(h_pnts.shape)
    h_coords = h_pnts.reshape(-1, shape[-1]).transpose()
    coords = np.divide(h_coords[:-1], h_coords[-1])
    shape[-1] -= 1
    return coords.transpose().reshape(*shape)


def get_h_pnt(pnts, weights):
    return np.hstack((pnts * weights[:, None], weights))
 


def intersect_lines(pos_v, dir_v):
    p13 = pos_v[0] - pos_v[1]
    d1343 = np.dot(p13, dir_v[1])
    d4321 = np.dot(dir_v[1], dir_v[0])
    d1321 = np.dot(p13, dir_v[0])
    d4343 = np.dot(dir_v[1], dir_v[1])
    d2121 = np.dot(dir_v[0], dir_v[0])
    denom = d2121 * d4343 - d4321 * d4321
    numer = d1343 * d4321 - d1321 * d4343
    mua = numer / denom
    return pos_v[0] + mua * dir_v[0], 1


def p_norm(v, p):
    return np.power(np.sum(np.power(np.abs(v), p)), 1. / p)


def project_point_to_line(pnt, pos_v, dir_v):
    return pos_v + np.dot(pnt - pos_v, dir_v) / np.dot(dir_v, dir_v) * dir_v
