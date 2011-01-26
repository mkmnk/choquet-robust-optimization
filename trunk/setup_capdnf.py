from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
ext_modules = [Extension("cap_dnf", ["cap_dnf.pyx"])]
setup(
    name = 'Capacity to dnf functions',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
    )

