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
import math

# 3rd party packages
import numpy as np

# project packages
from curve import Curve


def Arc(length):
    raise NotImplementedError


def Line(start_pnt, end_pnt):

    ctrl_pnts = np.hstack((np.vstack((start_pnt, end_pnt)), np.ones((2, 1))))
    knots = [0, 0, 1, 1]

    return Curve(ctrl_pnts, knots)


def Circle():

    r05 = math.sqrt(.5)
    
    x = [1, r05, 0, -r05, -1, -r05, 0, r05, 1]
    y = [0, r05, 1, r05, 0, -r05, -1, -r05, 0]
    w = [1, r05, 1, r05, 1, r05, 1, r05, 1]

    ctrl_pnts = zip(x, y, w)
    knots = [0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1]

    return Curve(ctrl_pnts, knots)
