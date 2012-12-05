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
    a = np.array([[1, -1], [1, -1]])
    a[0] *= [np.dot(dir_v[0], dir_v[i]) for i in range(2)]
    a[1] *= [np.dot(dir_v[1], dir_v[i]) for i in range(2)]
    b = [np.dot(dir_v[i], (pos_v[1] - pos_v[0])) for i in range(2)]
    t = np.linalg.solve(np.array(a), np.array(b))
    return pos_v + np.array([dir_v[i] * t[i] for i in range(2)])


def p_norm(v, p):
    return np.power(np.sum(np.abs(v), p), 1. / p)


def project_point_to_line(pnt, pos_v, dir_v):
    return pnt + np.dot(pos_v - pnt, dir_v) / np.dot(dir_v, dir_v) * dir_v
