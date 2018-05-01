import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_module = Extension(name="gait2dpi",
                       sources=["gait2dpi.pyx",
                                "gait2dpi_al.c"],
                       include_dirs=[numpy.get_include()])

setup(name="gait2dpi",
      cmdclass = {'build_ext': build_ext},
      ext_modules = [ext_module])
