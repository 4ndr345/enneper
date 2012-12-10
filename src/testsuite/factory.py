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


import unittest

import numpy as np

import enneper.curves as ec
import enneper.factory.nurbs as fcty


class TestNURBSFactory(unittest.TestCase):

    def setUp(self):
        wm = 0.707106781185
        kwargs = dict()
        kwargs['ctrl_points'] = np.zeros((5, 4))
        kwargs['ctrl_points'][0] = [0., 0., 1., 1.]
        kwargs['ctrl_points'][1] = [-wm, 0., wm, wm]
        kwargs['ctrl_points'][2] = [-1., 0., 0., 1.]
        kwargs['ctrl_points'][3] = [-wm, 0., -wm, wm]
        kwargs['ctrl_points'][4] = [0., 0., -1., 1.]
        kwargs['knots'] = [0, 0, 0, .5, .5, 1, 1, 1]
        kwargs['degree'] = 2
        self.curve = ec.NURBSCurve(**kwargs)

    def test_surface_of_revolution(self):
        pos_v = [0, 0, 0]
        dir_v = [0, 0, 1]
        phi = 2 * np.pi
        sphere = fcty.get_surface_of_revolution(self.curve, pos_v, dir_v, phi)
        opt = {'vertex_count_u': 50, 'vertex_count_v': 50}
        sphere.export_as_vrml('spehre', **opt)
