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
import enneper.foundation.geometry as efg


def get_ruled_surface(curves):
    curves = [ec.NURBSCurve(curve) for curve in curves]
    degree_u, degree_v = max(curves[0].degree, curves[1].degree), 1
    map(lambda curve: curve.elevate_degree(degree_u - curve.degree), curves)
    map(lambda i: curves[i % 2].merge_knots(curves[(i + 1) % 2].knots), [0, 1])
    knots_u, knots_v = curves[0].knots, np.array([0., 0., 1., 1.])
    ctrl_points = np.append(*[[curve.ctrl_points] for curve in curves])
    return es.NURBSSurface(ctrl_points, degree_u, degree_v, knots_u, knots_v)


def get_surface_of_revolution(position_vector, direction_vector, curve, theta):
    if 2 * np.pi < theta < 0:
        raise ValueError
    narcs = 4
    narcs = 3 if theta < 1.5 * np.pi else narcs
    narcs = 2 if theta < np.pi else narcs
    narcs = 1 if theta < 0.5 * np.pi else narcs
    dtheta = theta / narcs
    surface = es.NURBSSurface()
    if narcs == 2:
        surface.knots_u[3] = surface.knots_u[4] = 0
    if narcs == 3:
        surface.knots_u[3] = surface.knots_u[4] = 1 / 3.
        surface.knots_u[5] = surface.knots_u[6] = 2 / 3.
    if narcs == 4:
        surface.knots_u[3] = surface.knots_u[4] = .25
        surface.knots_u[5] = surface.knots_u[6] = .5
        surface.knots_u[7] = surface.knots_u[8] = .75
    surface.knots_u[:3] = np.zeros(3)
    surface.knots_u[(1 + 2 * narcs):(4 + 2 * narcs)] = np.ones(3)
    surface.knots_v[:] = curve.konts
    wm = np.cos(dtheta / 2)
    cosines = np.cos(dtheta * np.arange(1, narcs + 1))
    sines = np.sin(dtheta * np.arange(1, narcs + 1))
    for j, ctrl_point in enumerate(curve.ctrl_points):
        o = efg.project_to_line(position_vector, direction_vector, ctrl_point)
        x = ctrl_point - o
        y = np.cross(direction_vector, x)
        r = np.sqrt(np.dot(x, x))
        p0, t0 = surface.ctrl_points[0, j], _ = curve.ctrl_points[j], y
        for i in range(1, narcs):
            p2 = o + r * cosines[i] * x + r * sines[i] * y
            t2 = -sines[i] * x + cosines[i] * y
            pij = efg.intersect_line(p0, t0, p2, t2)
            wij = ctrl_point[j, -1] * wm
            surface.ctrl_points[i * 2 + 2, j] = p2, ctrl_point[j, -1]
            surface.ctrl_points[i * 2 + 1, j] = pij, wij
            p0, t0 = p2, t2 if i < narcs else p0, t0
