#!/usr/bin/env python
# -*- coding: utf-8 -*-

#***************************************************************************
#*   Copyright (C) 2012 by Andreas Kührmann [andreas.kuehrmann@gmail.com]  *
#*                                                                         *
#*   This program is free software; you can redistribute it and/or modify  *
#*   it under the terms of the GNU General Public License as published by  *
#*   the Free Software Foundation; either version 3 of the License, or     *
#*   (at your option) any later version.                                   *
#*                                                                         *
#*   This program is distributed in the hope that it will be useful,       *
#*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
#*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
#*   GNU General Public License for more details.                          *
#*                                                                         *
#*   You should have received a copy of the GNU General Public License     *
#*   along with this program; if not, write to the                         *
#*   Free Software Foundation, Inc.,                                       *
#*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
#***************************************************************************


import numpy as np


def find_span(u, degree, knots):
    """Determine the knot span index"""
    low, high = degree, knots.size - degree
    # Special case
    if u == knots[high - 1]:
        return high - 2
    # search
    return np.where(knots[low:high] <= u)[0][-1] + degree


def get_basis_functions(i, u, degree, knots):
    """ Computes the nonvanishing basis functions"""
    basis_func = np.ones(degree + 1)
    left = np.ones(degree + 1)
    right = np.ones(degree + 1)
    for j in range(1, degree + 1):
        left[j] = u - knots[i + 1 - j]
        right[j] = knots[i + j] - u
        cache = 0
        for k in range(j):
            tmp = basis_func[k] / (right[k + 1] + left[j - k])
            basis_func[k] = cache + right[k + 1] * tmp
            cache = left[j - k] * tmp
        basis_func[j] = cache
    return basis_func


def get_cartesian_points(h_pnts):
    shape = list(h_pnts.shape)
    h_coords = h_pnts.reshape(-1, shape[-1]).transpose()
    coords = np.divide(h_coords[:-1], h_coords[-1])
    shape[-1] -= 1
    return coords.transpose().reshape(*shape)


def get_homogeneous_point(pnts, weights):
    shape = list(pnts.shape)
    coords = pnts.reshape(-1, shape[-1]).transpose()
    weights = weights.reshape(-1)
    h_coords = np.vstack((np.multiply(coords, weights), weights))
    shape[-1] += 1
    return h_coords.transpose().reshape(*shape)
