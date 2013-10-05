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


###############################################################################
# Cython cimports
###############################################################################

cimport cython
from libc.stdlib cimport malloc, free
from numpy cimport ndarray


###############################################################################
# python function
###############################################################################

@cython.boundscheck(False)
@cython.cdivision(True)
def get_index(
    double u,
    int deg,
    ndarray[double] knots,
    ):

    """Determine the knot span index."""

    cdef:
        int low, high, mid, n

    # number of control points
    n = knots.size - deg - 1

    # handle special case
    if u == knots[n]:
        return n - 1

    # do binary search
    low, high = deg, n
    mid = (low + high) / 2
    while u < knots[mid] or u >= knots[mid + 1]:
        if u < knots[mid]:
            high = mid
        else:
            low = mid
        mid = (low + high) / 2
    return mid


@cython.boundscheck(False)
@cython.cdivision(True)
def calc_basis_funs(
    int index,
    double u,
    int deg,
    ndarray[double] konts,
    ndarray[double] basis_funs,
    ):

    """Computes the nonvanishing basis functions."""

    cdef:
        int i, j
        double saved, temp
        double *left = <double *> malloc((deg + 1) * cython.sizeof(double))
        double *right = <double *> malloc((deg + 1) * cython.sizeof(double))

    # calculate basis functions
    basis_funs[0] = 1
    for i in xrange(1, deg + 1):
        left[i] = u - konts[index + 1 - i]
        right[i] = konts[index + i] - u
        saved = 0
        for j in xrange(i):
            temp = basis_funs[j] / (right[j + 1] + left[i - j])
            basis_funs[j] = saved + right[j + 1] * temp
            saved = left[i - j] * temp
        basis_funs[i] = saved

    # release memory
    free(left)
    free(right)
