from setuptools import setup

setup(name='latticegpm',
      version='0.1',
      description='Sequence Space in the Protein Lattice Model Landscape',
      author='Zach Sailer',
      author_email='zachsailer@gmail.com',
      packages=['latticegpm'],
      install_requires=[
          'latticeproteins',
          'numpy',
      ],
      zip_safe=False)