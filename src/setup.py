"""
Setup script for PyCRFlux
$Id: setup.py,v 1.4 2009/02/14 04:06:44 oxon Exp $
"""

from distutils.core import setup

setup(name="PyCRFlux",
      version="1.0.3",
      description="Package for cosmic-ray flux calculation",
      author="Akira Okumura",
      author_email="oxon@ceres.phys.s.u-tokyo.ac.jp",
      packages=["cr_flux"],
      package_data={"cr_flux": ["data/*/*.dat", "data/*/*.gz"]}
      )
