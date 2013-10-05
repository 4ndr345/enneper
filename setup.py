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

from setuptools import setup

#from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

extra_args = []
# Comment/Uncomment the following line to disable/enable OpenMP for GCC-ish
# compilers.
extra_args = ["-fopenmp"]

exts = [Extension("foundation_cython", 
                  ["src/enneper/foundation.pyx"],
                  extra_compile_args=extra_args,
                  extra_link_args=extra_args,
                  include_dirs=[numpy.get_include()])
        ]

setup(name='enneper',
      version='0.1dev',
      description='A curve and surface library',
      author='Andreas Kührmann',
      author_email='andreas.kuehrmann@gmail.com',
      packages=['enneper'],
      package_dir={'': 'src/enneper'},
      install_requires=['numpy >= 1.5.1'],
      cmdclass = {'build_ext': build_ext},
      ext_modules = exts
      )