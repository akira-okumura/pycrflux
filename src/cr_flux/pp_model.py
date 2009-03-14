"""
Define p-p interaction models.
$Id: pp_model.py,v 1.6 2009/03/14 07:50:10 oxon Exp $
"""

import math
import numpy
from cparamlib import ParamModel
import ROOT
import pkg_resources

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
        self._sigma = None
        self._Eg    = numpy.zeros(0)
        self._Tp    = numpy.zeros(0)
        self.key1 = 0
        
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
            for iE in range(Eg.size):
                # Differntial inclusive cross section
                # sigma(Tp)dN/dEg = sigma(Tp)dN/dln(Eg)/Eg [mb/GeV]
                Esec = Eg[iE]/1e3
                Tp_ = Tp[iT]/1e3 # [GeV]
                Pp1 = (Tp_**2 + 2*matter.proton.mass/1e3*Tp_)**0.5 # [GeV/c/n]
                self._sigma[iE, iT] = pp_meson.pp_meson(Esec, Pp1, self.key1)

        # [b/GeV] => [mb/MeV]
        # do nothing
        
        return self._sigma*self.mul.multiplicity(par1, par2)
            
class Mori1997(PPModel):
    """
    M. Mori, The Astrophysical Journal 478 (1997) 225-232
    Porting from gamspec1.f developed by M. Mori (private communication)
    NOTE: Any bug in porting process is responsible to me not him!
    """
    def __init__(self, mul):
        PPModel.__init__(self, mul)
        self._sigma = None
        self._Eg    = numpy.zeros(0)
        self._Tp    = numpy.zeros(0)

        self.NBIN = 92
        self.NFILES = 76
        # Some arrays includes dummy entry because indices in the original
        # FORTRAN code has double naming conventions
        self.fact = [0., 1., 0.67, 0.33, 0., 0.] # first one is dummy
        self.dbin = 10**0.05 - 10**-0.05
        self.ext = ['dummy', '5m1','6m1','7m1','8m1','9m1', # first one is dummy
                    '0e0','1e0','2e0','3e0','4e0','5e0','6e0','7e0','8e0','9e0',
                    '0e1','1e1','2e1','3e1','4e1','5e1','6e1','7e1','8e1','9e1',
                    '0e2','1e2','2e2','3e2','4e2','5e2','6e2','7e2','8e2','9e2',
                    '0e3','1e3','2e3','3e3','4e3','5e3','6e3','7e3','8e3','9e3',
                    '0e4','1e4','2e4','3e4','4e4','5e4','6e4','7e4','8e4','9e4',
                    '0e5','1e5','2e5','3e5','4e5','5e5','6e5','7e5','8e5','9e5',
                    '0e6','1e6','2e6','3e6','4e6','5e6','6e6','7e6','8e6','9e6',
                    '0e7']
        self.spec = ROOT.TH2D("", "", self.NBIN + 1, -0.5, self.NBIN + 0.5, self.NFILES + 1, -0.5, self.NFILES + 0.5)
        self.spec.SetDirectory(0)

        # NUCRIN for <10GeV
        for i in range(1, 17):
            fname = pkg_resources.resource_filename("cr_flux", "data/mori/nucrin/gamma.hist." + self.ext[i])
            try:
                f = open(fname)
            except IOError:
                continue
            f.readline() # skip header
            lines = f.readlines()
            for j, line in enumerate(lines): 
                spec = float(line.split()[1])
                self.spec.Fill(j, i, spec)
        
        # for 3-8GeV, take average with Scaling
        fname = pkg_resources.resource_filename("cr_flux", "data/mori/gamma/gamma-sb.d.low")
        f = open(fname)
        for k in range(1, 6):
            f.readline() # skip header
            for j in range(1, self.NBIN):
                e1, s1 = [float(x) for x in f.readline().split()]
                s1 = s1*self.dbin*e1*2
                cont = self.spec.GetBinContent(j + 1, 10 + k + 1)
                #TODO: Something wrong ? should be checked by Mori-san
                cont = cont*self.fact[k] + s1*(1. - self.fact[k])
                self.spec.SetBinContent(j + 1, 10 + k + 1, cont)

        for i in range(17, self.NFILES + 1):
            fname = pkg_resources.resource_filename("cr_flux", "data/mori/lund/gamma.hist." + self.ext[i])
            f = open(fname)
            f.readline() # skip header
            for j, line in enumerate(f.readlines()):
                spec = float(line.split()[1])
                self.spec.Fill(j, i, spec)

    def Pi0Crs(self, tp):
        """
        Calculate pi0 (multiplicity)*(cross section) [mb]
        Port from Pi0Crs in gamspec1.f
        tp: Kinetic energy of proton [GeV]
        """
        PTHR = 0.78
        ampi0 = 0.1349764 # [GeV/c^2]
        amp = 0.938273 # [GeV/c^2]

        e = tp + amp
        p = (e**2 - amp**2)**0.5
        s = 2*amp*(e + amp)
        if p < PTHR:
            return 0.
        elif p < 0.96:
            eta = ((s - ampi0**2 - (2*amp)**2)**2 - 4*ampi0**2*(2*amp)**2)**0.5\
                  /(2*ampi0*s**0.5)
            return 0.032*eta**2 + 0.040*eta**6 + 0.047*eta**8
        elif p < 1.27:
            return 32.6*(p - 0.8)**3.21
        elif p < 8.0:
            return 5.40*(p - 0.8)**0.81
        elif p < 1000.:
            return 32.0*math.log(p) + 48.5/p**0.5 - 59.5
        else:
            return 163.08*(s/1876.)**0.21
    
    def sigma_gamma(self, par1, par2, Eg, Tp):
        if (self._sigma != None) and\
        (self._Eg.size == Eg.size) and numpy.all(self._Eg == Eg) and\
        (self._Tp.size == Tp.size) and numpy.all(self._Tp == Tp):
            # Return pre-calculated table to reduce CPU load
            return self._sigma*self.mul.multiplicity(par1, par2)
  
        self._sigma = numpy.zeros([Eg.size, Tp.size])
        self._Eg = numpy.copy(Eg)
        self._Tp = numpy.copy(Tp)

        # Ignore error messages in this for loops
        err = ROOT.gErrorIgnoreLevel
        ROOT.gErrorIgnoreLevel = ROOT.kSysError
        for iT in range(Tp.size):
            j = (math.log10(Tp[iT]/1e3) + 3.)/0.1 - 24.
            for iE in range(Eg.size):
                i = (math.log10(Eg[iE]/1e3) + 3.)/0.1 + 1.
                # Original value does not nclude [/MeV]. Add it in the last term
                self._sigma[iE, iT] = self.Pi0Crs(Tp[iT]/1e3)*self.spec.Interpolate(i, j)/(Eg[iE]/1e3*self.dbin)
        ROOT.gErrorIgnoreLevel = err
        
        self._sigma *= 1e-3 # [mb/GeV] to [mb/MeV]

        return self._sigma*self.mul.multiplicity(par1, par2)
    
class Huang2007(PPModel):
    """
    """
    #TODO: Implement this