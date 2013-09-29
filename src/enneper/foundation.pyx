#!/usr/bin/env cython
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

cimport cython
from libc.stdlib cimport malloc, free

def find_span(u, deg, knots):
    """Determine the knot span index"""
    lb, ub = deg, knots.size - deg
    # Special case
    if u == knots[ub - 1]:
        return ub - 2
    # search
    return np.where(knots[lb:ub] <= u)[0][-1] + deg



@cython.boundscheck(False)
def get_basis_funcs(int i, double u, int p, double[::1] U):
    
    cdef:
        int j, r, n
        float saved, temp
        double[::1] N
        double *left
        double *right
    
    n = p + 1
    
    N = np.empty(n, dtype=np.float64)
    # TODO: maybe use c array
    left = <double *>malloc(n*cython.sizeof(double))
    right = <double *>malloc(n*cython.sizeof(double))
    
    N[0] = 1
    
    for j in range(1, n):
        left[j] = u - U[i + 1 -j]
        right[j] = U[i + j] - u
        saved = 0
        for r in range(j):
            temp = N[r] / (right[r + 1] + left[j - r])
            N[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        N[j] = saved
    
    free(left)
    free(right)
    
    return np.asarray(N)
        