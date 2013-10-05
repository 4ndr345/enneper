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
def find_span(int n, int p, double u, ndarray[double] U):
    """
    Determine the knot span index
    INPUT: n, p, u, U
    Return: the knot span index
    """
    ###########################################################################
    # define types 
    ###########################################################################
    cdef:
        int low, high, mid
    ###########################################################################
    # handle special case
    ###########################################################################
    if u == U[n]:
        return n - 1
    ###########################################################################
    # do binary search
    ###########################################################################
    low, high = p, n
    mid = (low + high) / 2
    while u < U[mid] or u >= U[mid + 1]:
        if u < U[mid]:
            high = mid
        else:
            low = mid
        mid = (low + high) / 2
    return mid


@cython.boundscheck(False)
@cython.cdivision(True)
def basis_funs(int i, double u, int p, ndarray[double] U, ndarray[double] N):
    """
    Computes the nonvanishing basis functions.
    INPUT: i, u, p, U
    OUTPUT: N
    """
    ###########################################################################
    # define types and allocate memory for c arrays
    ###########################################################################
    cdef:
        int j, r, n = p + 1
        double saved, temp
        double *left = <double *> malloc(n * cython.sizeof(double))
        double *right = <double *> malloc(n * cython.sizeof(double))
    ###########################################################################
    # calculate basis functions
    ###########################################################################
    N[0] = 1
    for j in xrange(1, n):
        left[j] = u - U[i + 1 - j]
        right[j] = U[i + j] - u
        saved = 0
        for r in xrange(j):
            temp = N[r] / (right[r + 1] + left[j - r])
            N[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        N[j] = saved
    ###########################################################################
    # release memory
    ###########################################################################
    free(left)
    free(right)
