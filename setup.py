from setuptools import setup
from Cython.Build import cythonize

setup(
        # ext_modules = cythonize(["test3.py", "VectorDistributions/BinaryTrellis.py"])
        ext_modules = cythonize(["VectorDistributions/*.pyx"], annotate=True)
        )
