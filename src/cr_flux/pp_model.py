"""
Define p-p interaction models.
$Id: pp_model.py,v 1.3 2009/02/23 15:52:18 oxon Exp $
"""

import math
import numpy
from cparamlib import ParamModel

import matter
import pp_meson
import spectrum

class PPModel(object):
    """
    Base class for p-p collision model.
    """
    def __init__(self, mul):
        """
        mul: Multiplicity
        """
        self.mul = mul

    def gamma(self, gas, spec_cr, Fg):
        """
        gas: GasDensity or Cloud
        spec_cr: Spectrum of cosmic rays. You can use tuple or list to give
        more than one spectrum.
        Fg: Gamma-ray spectrum flux of which will NOT be overwritten        
        """
        if not isinstance(gas, matter.Gas):
            raise TypeError, "gas must be a subclass of Gas"

        if not (isinstance(spec_cr, list) or isinstance(spec_cr, tuple)):
            spec_cr = [spec_cr,]
        
        gas.cr_ab.expand_flux_list(spec_cr)
        
        zero = numpy.zeros(Fg.E.size)
        F    = numpy.zeros(Fg.E.size)
        
        for spec in spec_cr:
            par_cr = spec.par
            for par_gas in gas.ism_ab.abtable.keys():
                sigma = self.sigma_gamma(par_cr, par_gas, Fg.E, spec.E)
                for iT in range(spec.E.size):
                    dTp = spec.dEh[iT] + spec.dEl[iT]
                    dNp = 4*math.pi*spec.F[iT]*dTp # [/s/sr/cm^2/(MeV/n)] => [/s/cm^2]
                    dF = sigma[:, iT]*dNp*1e-27 # [/s/MeV]
                    ism_ab = gas.ism_ab.ism_abundance(par_gas)
                    F += dF*ism_ab
        
        if isinstance(gas, matter.GasDensity):
            if isinstance(Fg, spectrum.AbsoluteSpectrum):
                factor = gas.density
            else:
                raise TypeError, "Cannot calculate combination of GasDensity and %s" % Fg.__class__
        elif isinstance(gas, matter.Cloud):
            mass_par_p = 0. # total mass of gas par one proton [MeV/c^2]
            for par_gas in gas.ism_ab.abtable.keys():
                mass_par_p += par_gas.mass*gas.ism_ab.ism_abundance(par_gas)
            # calculate how many protons are in the cloud
            # 1 solar mass = 1.98892e30 [kg]
            # c = 2.99792458e8 [m/s]
            # e = 1.60e-19
            kg2MeV = 1./(1e6*1.60e-19/2.99792458e8**2)
            factor = gas.mass*1.98892e30*kg2MeV/mass_par_p
            if isinstance(Fg, spectrum.AbsoluteSpectrum):
                pass
            elif isinstance(Fg, spectrum.PointSourceSpectrum):
                # 1 [pc] = 3.086e18 [cm]
                pc2cm = 3.0857e18
                factor /= 4*math.pi*(gas.dist*pc2cm)**2 # change to /cm^2
            else:
                raise TypeError, "Cannot calculate combination of GasDensity and %s" % Fg.__class__

        return Fg.__class__(Fg.par, Fg.E, Fg.dEl, Fg.dEh, F*factor, zero, zero)
        
class Kamae2006(PPModel):
    """
    T. Kamae, et al., The Astrophysical Journal 647 (2006) 692-708
    """
    def __init__(self, mul):
        PPModel.__init__(self, mul)
        self._sigma = None
        self._Eg    = numpy.zeros(0)
        self._Tp    = numpy.zeros(0)
        
    def sigma_gamma(self, par1, par2, Eg, Tp):
        """
        Calculate inclusive cross section of gamma rays in unit of [mb/MeV]

        par1: particle type 1
        par2: particle type 2
        Eg: An array of energies of gamma rays in unit of [MeV]
        Tp: An array of kinetic energies of protons in unit of [MeV]
        """
        if (self._sigma != None) and\
        (self._Eg.size == Eg.size) and numpy.all(self._Eg == Eg) and\
        (self._Tp.size == Tp.size) and numpy.all(self._Tp == Tp):
            # Return pre-calculated table to reduce CPU load
            return self._sigma*self.mul.multiplicity(par1, par2)
        
        self._sigma = numpy.zeros([Eg.size, Tp.size])
        self._Eg = numpy.copy(Eg)
        self._Tp = numpy.copy(Tp)

        for iT in range(Tp.size):
            model = ParamModel.ParamModel(Tp[iT]/1e3) # [MeV] to [GeV]
            for iE in range(Eg.size):
                # Differntial inclusive cross section
                # sigma(Tp)dN/dEg = sigma(Tp)dN/dln(Eg)/Eg [mb/GeV]
                self._sigma[iE, iT] = model.sigma_incl_tot(Eg[iE]/1e3, Tp[iT]/1e3)/(Eg[iE]/1e3)
        
        # [mb/GeV] => [mb/MeV]
        self._sigma /= 1e3

        return self._sigma*self.mul.multiplicity(par1, par2)

class Dermer1986(PPModel):
    """
    Based on pp_meson.f in GALPROP
    [1] Badhwar, Stephens, & Golden 1977, Phys. Rev. D 15, 820
    [2] Dermer 1986, A&A 157, 223;  ApJ 307, 47
    [3] Mori 1997, ApJ 478, 225
    [4] Moskalenko & Strong 1998, ApJ 493, 694
    [5] Stecker 1970, Astrophys. Spa. Sci. 6, 377
    [6] Stephens & Badhwar 1981, Astrophys. Spa. Sci. 76, 213
    """
    def __init__(self, mul):
        PPModel.__init__(self, mul)
        self._sigma = {}
        self.key1 = 0
        
    def sigma_gamma(self, par1, par2, Eg, Tp):
        """
        Calculate inclusive cross section of gamma rays in unit of [mb/MeV]

        par1: particle type 1
        par2: particle type 2
        Eg: An array of energies of gamma rays in unit of [MeV]
        Tp: An array of kinetic energies of protons in unit of [MeV]
        """
  
        sigma = numpy.zeros([Eg.size, Tp.size])

        for iT in range(Tp.size):
            for iE in range(Eg.size):
                # Differntial inclusive cross section
                # sigma(Tp)dN/dEg = sigma(Tp)dN/dln(Eg)/Eg [mb/GeV]
                Esec = Eg[iE]/1e3
                T_tot = Tp[iT]*par1.A/1e3 # [GeV]
                Pp1 = (T_tot**2 + 2*par1.mass/1e3*T_tot)**0.5 # [GeV/c/nucleUS]
                NA1 = par1.A
                NA2 = par2.A
                sigma[iE, iT] = pp_meson.pp_meson(Esec, Pp1, NA1, NA2, self.key1)

        # [b/GeV] => [mb/MeV]
        # do nothing
        
        return sigma
            
class Mori1997(PPModel):
    """
    M. Mori, The Astrophysical Journal 478 (1997) 225-232
    """
    
class Huang2007(PPModel):
    """
    """
