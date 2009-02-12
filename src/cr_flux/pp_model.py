"""
Define p-p interaction models.
$Id: pp_model.py,v 1.1 2009/02/12 18:15:15 oxon Exp $
"""

import math
import numpy
from cparamlib import ParamModel

import matter
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
        Pure virtual method
        """
        raise NotImplementedError
    
class Kamae2006(PPModel):
    """
    T. Kamae, et al., The Astrophysical Journal 647 (2006) 692-708
    """
    def sigma_gamma(self, Eg, Tp):
        """
        Calculate inclusive cross section of gamma rays in unit of [mb/MeV]

        Eg: An array of energies of gamma rays in unit of [MeV]
        Tp: An array of kinetic energies of protons in unit of [MeV]
        """
        sigma = numpy.zeros([Eg.size, Tp.size])

        for iT in range(Tp.size):
            model = ParamModel.ParamModel(Tp[iT]/1e3) # [MeV] to [GeV]
            for iE in range(Eg.size):
                # Differntial inclusive cross section
                # sigma(Tp)dN/dEg = sigma(Tp)dN/dln(Eg)/Eg [mb/GeV]
                sigma[iE, iT] = model.sigma_incl_tot(Eg[iE]/1e3, Tp[iT]/1e3)\
                                                     /(Eg[iE]/1e3)
        # [mb/GeV] => [mb/MeV]
        sigma /= 1e3
        
        return sigma

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
            spec_cr = (spec_cr,)
        
        gas.cr_ab.expand_flux_list(spec_cr)
        
        zero = numpy.zeros(Fg.E.size)
        F    = numpy.zeros(Fg.E.size)
        
        for spec in spec_cr:
            par_cr = spec.par
            sigma = self.sigma_gamma(Fg.E, spec.E)
            for iT in range(spec.E.size):
                dTp = spec.dEh[iT] + spec.dEl[iT]
                dNp = 4*math.pi*spec.F[iT]*dTp # [/s/sr/cm^2/(MeV/n)] => [/s/cm^2]
                dF = sigma[:, iT]*dNp*1e-27 # [/s/MeV]

                for par_gas in gas.ism_ab.abtable.keys():
                    ism_ab = gas.ism_ab.ism_abundance(par_gas)
                    F += dF*self.mul.multiplicity(par_cr, par_gas)*ism_ab
        
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
    
class Mori1997(PPModel):
    """
    M. Mori, The Astrophysical Journal 478 (1997) 225-232
    """
    
class Huang2007(PPModel):
    """
    """
