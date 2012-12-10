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

import enneper.foundation as fdn


__all__ = ['NURBSSurface']


class NURBSSurface(object):

    _init_args = {
        (): '_empty_initializer',
        ('surface',): '_copy_initializer',
        ('shape', 'degree_u', 'degree_v'): 'resize',
        ('ctrl_points', 'degree_u', 'degree_v', 'knots_u', 'knots_v'): '_init',
    }

    def __init__(self, *args):
        key = tuple(type(arg).__name__ for arg in args)
        if key not in NURBSSurface._init_args:
            raise TypeError
        self._ctrl_points = None
        self._degree_u = None
        self._degree_v = None
        self._knots_u = None
        self._knots_v = None
        getattr(self, NURBSSurface._init_args[key])(*args)

    def _empty_initializer(self):
        pass

    def _copy_initializer(self, surface):
        self._ctrl_points = surface.ctrl_points.copy()
        self._degree_u = surface.degree
        self._degree_v = surface.degree
        self._knots_u = surface.knots_u.copy()
        self._knots_v = surface.knots_v.copy()

    def _init(self, ctrl_points, knots_u, knots_v, degree_u, degree_v):
        ctrl_points = np.asarray(ctrl_points)
        knots_u = np.asarray(knots_u)
        knots_v = np.asarray(knots_v)
        if ctrl_points.shape[0] + degree_u + 1 != knots_u.size:
            raise ValueError
        if ctrl_points.shape[1] + degree_v + 1 != knots_v.size:
            raise ValueError
        self._ctrl_points = ctrl_points
        self._degree_u = degree_u
        self._degree_v = degree_v
        self._knots_u = knots_u
        self._knots_v = knots_v

    @ property
    def ctrl_points(self):
        return self._ctrl_points

    @ property
    def degree_u(self):
        return self._degree_u

    @ property
    def degree_v(self):
        return self._degree_v

    @ property
    def knots_u(self):
        return self._knots_u

    @ property
    def knots_v(self):
        return self._knots_v

    def export_as_vrml(self, fname, **opt):
        color = np.asarray(opt.get('color', [255, 255, 255]))
        n = opt.get('vertex_count_u', 20)
        m = opt.get('vertex_count_v', 20)
        u_s = opt.get('u_start', 0)
        u_e = opt.get('u_end', 1)
        v_s = opt.get('v_start', 0)
        v_e = opt.get('v_end', 1)
        fname = fname if fname[:-4] == '.wrl' else ''.join((fname, '.wrl'))
        n = n if n > 1 else 2
        m = m if m > 1 else 2
        p_max = self.evaluate_at(u_s, v_s)
        p_min = self.evaluate_at(u_s, v_s)
        if self.ctrl_points.shape[2] == 4:
            write_point = lambda p: "\t\t\t\t{0} {1} {2},\n".format(*p)
        elif self.ctrl_points.shape[2] == 3:
            write_point = lambda p: "\t\t\t\t{0} {1} {0},\n".format(*p)
        else:
            raise Exception
        with open(fname, 'w') as fout:
            fout.write('#VRML V2.0 utf8\n')
            fout.write('\nGroup {\n')
            fout.write('\n  children [\n')
            fout.write('\tDEF T Transform {\n')
            fout.write('\t  children [\n')
            fout.write('\t\tShape {\n')
            fout.write('\t\t appearance Appearance {\n')
            fout.write('\t\t\tmaterial Material{ diffuseColor ')
            color = np.divide(color, np.asarray([255, 255, 255]))
            fout.write('{0} {1} {2} '.format(*color))
            fout.write('}\n')
            fout.write('\t\t }\n')
            fout.write('\t\t geometry IndexedFaceSet {\n')
            fout.write('\t\t\tsolid FALSE\n')
            fout.write('\t\t\tcoord Coordinate {\n')
            fout.write('\t\t\t point [\n')
            for u in np.linspace(u_s, u_e, n):
                for v in np.linspace(v_s, v_e, m):
                    p = self.evaluate_at(u, v)
                    fout.write(write_point(p))
                    indices, = np.where(p < p_min)
                    p_min[indices] = p[indices]
                    indices, = np.where(p > p_max)
                    p_max[indices] = p[indices]
            fout.write('\t\t\t ]\n')
            fout.write('\t\t\t}\n')
            fout.write('\t\t\t coordIndex [\n')
            for i in range(n - 1):
                for j in range(m - 1):
                    fout.write('\t\t\t\t{0}, '.format(i * m + j))
                    fout.write('{0}, '.format(i * m + j + 1))
                    fout.write('{0}, -1,\n'.format((i + 1) * m + j))
                    fout.write('\t\t\t\t{0}, '.format(i * m + j + 1))
                    fout.write('{0}, '.format((i + 1) * m + j + 1))
                    fout.write('{0}, -1,\n'.format((i + 1) * m + j))
            fout.write('\t\t\t ]\n')
            fout.write('\t\t\t}\n')
            fout.write('\t\t}\n')
            fout.write('\t ]\n')
            fout.write('\t}\n')
            fout.write('  ]\n')
            fout.write('}\n')
            p_mid = (p_max + p_min) / 2.
            x_axis, y_axis = (p_max - p_min)[:2]
            axis = y_axis if np.all(x_axis < y_axis) else x_axis
            fout.write('Viewpoint {\n\t position ')
            fout.write('{0} {1} '.format(p_mid[0], p_mid[1]))
            fout.write('{0}'.format(p_max[2] + axis * 2))
            fout.write('\n\t description \"top\"\n}\n')
            fout.write('NavigationInfo { type \"EXAMINE\" }\n')

    def evaluate_at(self, u, v, homogenous=False):
        degree_u, degree_v, = self.degree_u, self.degree_v
        knots_u, knots_v = self.knots_u, self.knots_v
        span_u = fdn.nurbs.find_span(u, degree_u, knots_u)
        n_u = fdn.nurbs.get_basis_functions(span_u, u, degree_u, knots_u)
        span_v = fdn.nurbs.find_span(v, degree_v, knots_v)
        n_v = fdn.nurbs.get_basis_functions(span_v, v, degree_v, knots_v)
        index_u, index_v = span_u - degree_u, span_v - degree_v
        tmp = np.zeros((degree_v + 1, self.ctrl_points.shape[2]))
        for i in range(degree_v + 1):
            for j in range(degree_u + 1):
                tmp[i] += n_u[j] * self.ctrl_points[index_u + j, index_v + i]
        point = np.zeros(self.ctrl_points.shape[2])
        # todo: optimize loop: use np.sum and np.multiply
        for i in range(degree_v + 1):
            point += n_v[i] * tmp[i]
        if homogenous:
            return point
        point /= point[-1]
        return point[0:-1]

    def resize(self, shape, degree_u, degree_v):
        self._ctrl_points = np.zeros(shape)
        self._degree_u = degree_u
        self._degree_v = degree_v
        self._knots_u = np.zeros(shape[0] + degree_u + 1)
        self._knots_v = np.zeros(shape[1] + degree_v + 1)
