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
import enneper.foundation as fdn


npa = lambda array: np.array(array)


def get_ruled_surface(curves):
    curves = [ec.NURBSCurve(curve) for curve in curves]
    degree_u, degree_v = max(curves[0].degree, curves[1].degree), 1
    map(lambda curve: curve.elevate_degree(degree_u - curve.degree), curves)
    map(lambda i: curves[i % 2].merge_knots(curves[(i + 1) % 2].knots), [0, 1])
    knots_u, knots_v = curves[0].knots, npa([0., 0., 1., 1.])
    ctrl_points = np.append(*[[curve.ctrl_points] for curve in curves])
    return es.NURBSSurface(ctrl_points, degree_u, degree_v, knots_u, knots_v)


def get_surface_of_revolution(curve, pos_v, dir_v, phi):
    narcs = 3 if phi < 1.5 * np.pi else 4
    narcs = 2 if phi < np.pi else narcs
    narcs = 1 if phi < .5 * np.pi else narcs
    n = 2 * narcs + 1
    m = curve.ctrl_points.shape[0]
    dphi = phi / float(narcs)
    cosines = np.cos(dphi * np.arange(narcs + 1))
    sines = np.sin(dphi * np.arange(narcs + 1))
    wm = np.cos(.5 * dphi)
    pj = fdn.nurbs.get_cartesian_points(curve.ctrl_points)
    wj = curve.ctrl_points[:, -1]
    surface = es.NURBSSurface((n, m, 4), 2, curve.degree)
    for i in range(1, narcs):
        surface.knots_u[i * 2 + 1] = i / float(narcs)
        surface.knots_u[i * 2 + 2] = i / float(narcs)
    surface.knots_u[n:(n + 3)] = np.ones(3)
    surface.knots_v[:] = curve.knots
    for j in range(m):
        o = fdn.math.project_point_to_line(pj[j], pos_v, dir_v)
        x = pj[j] - o
        y = np.cross(dir_v, x)
        r = fdn.math.p_norm(x, 2)
        index = 0
        if r < 1e-7:
            h_point = fdn.nurbs.get_homogeneous_points(o, wj[j])
            for i in range(narcs * 2 + 1):
                surface.ctrl_points[i, j] = h_point
            continue
        h_point = fdn.nurbs.get_homogeneous_points(pj[j], wj[j])
        surface.ctrl_points[0, j] = h_point
        p0 = pj[j]
        t0 = y
        for i in range(1, narcs + 1):
            p2 = o + r * cosines[i] * x + r * sines[i] * y
            t2 = -sines[i] * x + cosines[i] * y
            h_point = fdn.nurbs.get_homogeneous_points(p2, wj[j])
            surface.ctrl_points[index + 2, j] = h_point
            pnt, _ = fdn.math.intersect_lines(npa([p0, p2]), npa([t0, t2]))
            h_point = fdn.nurbs.get_homogeneous_points(pnt, wm * wj[j])
            surface.ctrl_points[index + 1, j] = h_point
            index += 2
            if i < narcs:
                p0, t0 = p2, t2
    return surface
