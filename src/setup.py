"""
"""
# $Id: setup.py,v 1.1 2009/02/12 18:10:39 oxon Exp $

from distutils.core import setup

setup(name="PyCRFlux",
      version="1.0.0",
      description="Package for cosmic-ray flux calculation",
      author="Akira Okumura",
      author_email="oxon@ceres.phys.s.u-tokyo.ac.jp",
      packages=["cr_flux"],
      package_data={"cr_flux": ["data/*/*.dat", "data/*/*.gz"]}
      )
