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
