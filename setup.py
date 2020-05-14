from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

sources = ['calc_delay.pyx', 'calc_delay_error.c']

setup(
  name='delay time error calculation',
  ext_modules=[Extension('_dt_error', sources, include_dirs=[numpy.get_include()])],
  cmdclass={'build_ext': build_ext},
)
