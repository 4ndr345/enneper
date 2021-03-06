#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ***************************************************************************
# *   Copyright (C) 2013 by Andreas Kührmann [andreas.kuehrmann@gmail.com]  *
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


# standard packages
from __future__ import division
import math

# 3rd party packages
import numpy as np

# project packages
from curve import Curve


def CircularArc(length):

    if abs(length) <= np.pi * .5:
        narcs = 1
        knots = [0, 0, 0, 1, 1, 1]
    elif abs(length) <= np.pi:
        narcs = 2
        knots = [0, 0, 0, .5, .5, 1, 1, 1]
    elif abs(length) <= 1.5 * np.pi:
        narcs = 3
        knots = [0, 0, 0, .33333333, .33333333, .66666666, .66666666, 1, 1, 1]
    else:
        narcs = 4
        knots = [0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1]

    step_size  = length / narcs
    angels = [i * step_size for i in range(narcs + 1)]

    x = np.cos(angels)
    y = np.sin(angels)
    w = math.cos(step_size * 0.5)

    ctrl_pnts = np.ones((2 * narcs + 1, 3))
    ctrl_pnts[::2, 0] = x
    ctrl_pnts[::2, 1] = y

    tangents = np.ones((len(angels), 3))
    tangents[:, 0] = x - y
    tangents[:, 1] = y + x

    bounding_box = np.cross(ctrl_pnts[::2], tangents)
    ctrl_pnts[1::2] = np.cross(bounding_box[:-1], bounding_box[1:])

    ctrl_pnts[1::2] = ctrl_pnts[1::2] / ctrl_pnts[1, -1] * w

    return Curve(ctrl_pnts, knots)


def Line(start_pnt, end_pnt):

    ctrl_pnts = np.hstack((np.vstack((start_pnt, end_pnt)), np.ones((2, 1))))
    knots = [0, 0, 1, 1]

    return Curve(ctrl_pnts, knots)
