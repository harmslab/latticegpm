try: 
    from setuptools import setup
except:
    from distutils.core import setup

setup(name='latticegpm',
      version='0.1',
      description='Sequence Space in the Protein Lattice Model Landscape',
      author='Zach Sailer',
      author_email='zachsailer@gmail.com',
      packages=['latticegpm'],
      install_requires=[
          'latticeproteins',
          'numpy',
          'gpmap',
      ],
      zip_safe=False)
