from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
ext_modules = [Extension("choquet_base", ["choquet_base.pyx"])]
setup(
    name = 'Choquet base functions',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
    )

