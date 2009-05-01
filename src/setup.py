"""
Setup script for PyCRFlux
$Id: setup.py,v 1.9 2009/05/01 15:21:46 oxon Exp $
"""

from numpy.distutils.core import Extension, setup

ext = Extension(name="cr_flux.pp_meson", sources=["ext/pp_meson_mod.f",])

setup(name="PyCRFlux",
      version="1.2.2",
      description="Package for cosmic-ray flux calculation",
      author="Akira Okumura",
      author_email="oxon@ceres.phys.s.u-tokyo.ac.jp",
      packages=["cr_flux"],
      package_data={"cr_flux": ["data/*/*.dat", "data/*/*.gz", "data/mori/*/gamma*"]},
      ext_modules = [ext,]
      )
