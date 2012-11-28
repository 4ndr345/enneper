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


def get_ruled_surface(cams):
    cams = [ec.NURBSCurve(cam) for cam in cams]
    degree_u, degree_v = max(cams[0].degree, cams[1].degree), 1
    map(lambda cam: cam.elevate_degree(degree_u - cam.degree), cams)
    map(lambda i: cams[i % 2].merge_knots(cams[(i + 1) % 2].knots), range(2))
    knots_u, knots_v = cams[0].knots, np.array([0., 0., 1., 1.])
    ctrl_points = np.append(*[[cam.ctrl_points] for cam in cams])
    return es.NURBSSurface(ctrl_points, degree_u, degree_v, knots_u, knots_v)
