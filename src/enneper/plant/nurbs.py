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

import enneper.curves as ec
import enneper.surfaces as es


def get_ruled_surface(curves):
    curves = [ec.NURBSCurve(curve) for curve in curves]
    degree_u, degree_v = max(curves[0].degree, curves[1].degree), 1
    map(lambda curve: curve.elevate_degree(degree_u - curve.degree), curves)
    map(lambda i: curves[i % 2].merge_knots(curves[(i + 1) % 2].knots), [0, 1])
    knots_u, knots_v = curves[0].knots, np.array([0., 0., 1., 1.])
    ctrl_points = np.append(*[[curve.ctrl_points] for curve in curves])
    return es.NURBSSurface(ctrl_points, degree_u, degree_v, knots_u, knots_v)
