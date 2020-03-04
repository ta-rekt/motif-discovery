from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize('onionDecompose_cython.pyx', language_level = "3", annotate=True))
