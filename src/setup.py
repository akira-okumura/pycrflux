"""
Setup script for PyCRFlux
$Id: setup.py,v 1.11 2009/05/19 17:05:38 oxon Exp $
"""

from numpy.distutils.core import Extension, setup

ext = Extension(name="cr_flux.pp_meson", sources=["ext/pp_meson_mod.f",])

setup(name="PyCRFlux",
      version="1.3.1",
      description="Package for cosmic-ray flux calculation",
      author="Akira Okumura",
      author_email="oxon@ceres.phys.s.u-tokyo.ac.jp",
      packages=["cr_flux"],
      package_data={"cr_flux": ["data/*/*.dat", "data/*/*.gz", "data/mori/*/gamma*"]},
      ext_modules = [ext,]
      )
