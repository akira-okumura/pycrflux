"""
"""
#$Id: spectrum.py,v 1.20 2010/02/27 04:19:05 oxon Exp $

import copy
import decimal
import math
import numpy
import pkg_resources
import pyfits
import ROOT

import matter

class Spectrum(object):
    """
    Base class for energy spectra of particles. 
    """
    def __init__(self, par, E, dEl, dEh, F=None, dFl=None, dFh=None):
        """
        par (Particle): Particle type
        E (numpy.ndarray): Kinetic energy [MeV] or [MeV/n]
        dEl/dEh (numpy.ndarray): Asymmetric errors of E
        F (numpy.ndarray): Flux dN/dE [/cm^2/s/MeV] or [/cm^2/s/sr/MeV]
        dFl/dFh (numpy.ndarray): Asymmetric errors of E
        """
        if not isinstance(par, matter.Particle):
            raise TypeError, "par must be Particle class"
        if not isinstance(E, numpy.ndarray) or\
           not isinstance(dEl, numpy.ndarray) or\
           not isinstance(dEh, numpy.ndarray) or\
           not (isinstance(F, numpy.ndarray) or F == None) or\
           not (isinstance(dFl, numpy.ndarray) or dFl == None) or\
           not (isinstance(dFh, numpy.ndarray) or dFh == None): 
            raise TypeError, "numpy.ndarray must be given"
        """
        if not (E.size == dEl.size == dEh.size) or\
           not ((F   != None) and F.size == E.size) or\
           not ((dFl != None) and dFl.size == E.size) or\
           not ((dFh != None) and dFh.size == E.size):
            raise TypeError, "The lengths of arrays are different"
        """
        if not (E.size == dEl.size == dEh.size) or\
        ((F != None) and F.size != E.size) or\
        ((dFl != None) and dFl.size != E.size) or\
        ((dFh != None) and dFh.size != E.size):
            raise TypeError, "The lengths of arrays are different" 
        
        self.par = copy.copy(par)
        self.E   = copy.copy(E)
        self.dEl = copy.copy(dEl)
        self.dEh = copy.copy(dEh)
        self.F = numpy.zeros(self.E.size) if F == None else copy.copy(F) 
        self.dFl = numpy.zeros(self.E.size) if dFl == None else copy.copy(dFl) 
        self.dFh = numpy.zeros(self.E.size) if dFh == None else copy.copy(dFh) 
                
        # weight index for SED plot (i.e. idx == 2.0 means E^2 dN/dE)
        self.__idx = decimal.Decimal("2.0")
        self.show_xerr = False
        self.init_graph()
        
    def __add__(self, other):
        if self.par != other.par:
            raise TypeError, "Different particle types are given"
        elif not numpy.all(self.E == other.E) or\
             not numpy.all(self.dEl == other.dEl) or\
             not numpy.all(self.dEh == other.dEh):
            raise TypeError, "Different energy binnings are given"
        elif self.__class__ != other.__class__:
            raise TypeError, "Different spectrum classes are given"
        else:
            return self.__class__(self.par, self.E, self.dEl, self.dEh,
                                  self.F + other.F,
                                  (self.dFl**2 + other.dFl**2)**0.5,
                                  (self.dFh**2 + other.dFh**2)**0.5)
    
    def __mul__(self, factor):
        return self.__class__(self.par, self.E, self.dEl, self.dEh,
                              self.F*factor, self.dFl*factor, self.dFh*factor)

    def get_mean_energy(self):
        """
        Calculate weighted mean energy
        """
        n = 0.
        e = 0.
        for i in range(self.E.size):
            dn = self.F[i]*(self.dEl[i] + self.dEh[i]) # dN/dE*dE
            n += dn
            e += self.E[i]*dn # dN*E
        
        return e/n

    def init_graph(self):
        """
        """
        self.graph = ROOT.TGraphAsymmErrors(0)
        self.graph.SetTitle("")
        for i, (E, dEl, dEh, F, dFl, dFh) in \
        enumerate(zip(self.E, self.dEl, self.dEh, self.F, self.dFl, self.dFh)):
            weight = E**float(self.__idx)
            self.graph.SetPoint(i, E, F*weight)
            if self.show_xerr:
                self.graph.SetPointError(i, dEl, dEh, dFl*weight, dFh*weight)
            else:
                self.graph.SetPointError(i, 0, 0, dFl*weight, dFh*weight)
        
        xax = self.graph.GetXaxis()
        yax = self.graph.GetYaxis()
        xax.CenterTitle()
        yax.CenterTitle()
        
        # Change the unit of energy
        if isinstance(self.par, matter.Photon):
            eunit = "MeV" # per particle
            xax.SetTitle("Energy (MeV)")
        elif isinstance(self.par, matter.ChargedLepton):
            eunit = "MeV" # per nucleon
            xax.SetTitle("Kinetic Energy (MeV)")
        else:
            eunit = "(MeV/n)" # per nucleon
            xax.SetTitle("Kinetic Energy (MeV/n)")
        
        # Change the unit of acceptance
        if isinstance(self, AbsoluteSpectrum):
            aunit = "s^{-1}"
        elif isinstance(self, DiffuseSpectrum):
            aunit = "cm^{-2}s^{-1}sr^{-1}"
        elif isinstance(self, PointSourceSpectrum):
            aunit = "cm^{-2}s^{-1}"
        else:
            aunit = "s^{-1}"
          
        if self.__idx.normalize() == decimal.Decimal("0"):
            yax.SetTitle("dN/dE (%s%s^{-1})" % (aunit, eunit))
        elif self.__idx.normalize() == decimal.Decimal("1"):
            yax.SetTitle("E dN/dE (%s)" % aunit)
        elif self.__idx.normalize() == decimal.Decimal("2"):
            yax.SetTitle("E^{%s} dN/dE (%s%s)" % \
                         (str(self.__idx.normalize()), aunit, eunit))
        else:
            yax.SetTitle("E^{%s} dN/dE (%s%s^{%s})" % \
                         (str(self.__idx.normalize()), aunit,
                          eunit, str((self.__idx.normalize() - 1))))
            
    def scale(self, factor):
        """
        Scale the flux and its errors. 
        """
        self.F   *= factor
        self.dFl *= factor
        self.dFh *= factor
        self.init_graph()
    
    def set_index(self, idx):
        """
        """
        self.__idx = decimal.Decimal(str(idx))
        self.init_graph()

    def fill_power_law(self, normalization, scale, index):
        """
        Create power law spectrum
        dN/dE = N x (E/E0)^(index)
    
        normalization: N, unit depends on the class type
        scale: Energy scale E0 [MeV]
        index: Spectrum index (usually negative)
        """
        self.F   = normalization*(self.E/scale)**index
        self.dFh *= 0.
        self.dFl *= 0.
        
        self.init_graph()
        
    def solar_modulate(self, phi):
        """
        Calculate the effect from solar modulation
        
        Original code is modulate.cc in GALPROP
        L. J. Gleeson and W. I. Axford,
        The Astrophysical Journal 154 (1968) 1011-1026
        """
        if isinstance(self.par, matter.Nucleus):
            charge = self.par.Z
            m = self.par.mass/self.par.A
            T = self.E + abs(charge)*phi/self.par.A
        elif isinstance(self.par, matter.ChargedLepton):
            charge = abs(self.par.charge)
            m = self.par.mass
            T = self.E + abs(charge)*phi
        else:
            return

        spec = self.__class__(self.par, self.E, self.dEl, self.dEh, self.F, self.dFl, self.dFh)

        np = len(spec.F)
        cr_density = spec.F*spec.E**2
        density = numpy.zeros(np)

        for i in range(np):
            for j in range(np):
                # NOTE! j will be used outside the loop
                if T[i] < spec.E[j]:
                    break
            if j == np - 1:
                density[i] = cr_density[i]
                break
            y = cr_density[j - 1] + (T[i] - spec.E[j - 1])*(cr_density[j] - cr_density[j - 1])/(spec.E[j] - spec.E[j - 1])
            density[i] = y*spec.E[i]*(spec.E[i] + 2*m)/T[i]/(T[i] + 2*m)*(spec.E[i]/T[i])**2
            
        spec.F = density/spec.E**2
        
        spec.init_graph()
        
        return spec

class AbsoluteSpectrum(Spectrum):
    """
    """

class DiffuseSpectrum(Spectrum):
    """
    """
        
class PointSourceSpectrum(Spectrum):
    """
    """

class FluxModelArchive(object):
    """
    Archive of cosmic-ray flux of simulation and experimental data.
    """ 
    def galprop(self, par, R, z, galdef="54_59Xvarh7S"):
        """
        par: particle type
        R: Distance from the galactic center [kpc]
        z: Distance from the galactic plane [kpc]

        Original data is available on SLAC server
        /nfs/farm/g/glast/u33/diffuse/galprop/FITS/nuclei_full_54_59Xvarh7S.gz 
        """

        fname = pkg_resources.resource_filename("cr_flux", "data/galprop/nuclei_full_%s.gz" % galdef)

        if par == matter.proton: # H
            pid = 3
        elif par == matter.deutron: # D
            pid = 4
        elif par == matter.Nucleus(3, 2): # 3He
            pid = 5
        elif par == matter.alpha: # He
            pid = 6
        elif par == matter.positron:
            pid = 0
        elif par == matter.electron:
            pass
        else:
            raise TypeError, "Invalid particle type" 
        
        hdulist = pyfits.open(fname)
        hdu = hdulist[0]
        head = hdu.header
        if par == matter.electron:
            data = hdu.data[1] + hdu.data[2]
        else:
            data = hdu.data[pid]
        
        E0 = 10**head["CRVAL3"] # [MeV]
        Estep = 10**head["CDELT3"]
        E   = E0*Estep**numpy.arange(head["NAXIS3"])
        dEh = E*Estep**0.5 - E
        dEl = E - E/Estep**0.5
        F   = numpy.zeros(head["NAXIS3"])

        dR = head["CDELT1"] # [kpc]
        dz = head["CDELT2"] # [kpc]
        R0 = head["CRVAL1"]
        z0 = head["CRVAL2"]
        nR = head["NAXIS1"]
        nz = head["NAXIS2"]
        
        h2 = ROOT.TH2D("", "", nR, R0 - 0.5*dR, R0 + (nR - 0.5)*dR,\
                       nz, z0 - 0.5*dz, z0 + (nz - 0.5)*dz)
        h2.SetDirectory(0)

        if R < R0 or R > R0 + dR*(nR - 1):
            raise ValueError, "R is out of range"
        if z < z0 or z > z0 + dz*(nz - 1):
            raise ValueError, "z is out of range"
    
        nearR = int(round((R - R0)/dR)) # nearest index of R
        nearz = int(round((z - z0)/dz)) # nearest index of z
            
        for i in range(E.size):
            for iR in range(nearR - 1, nearR + 2):
                if iR < 0 or nR <= iR:
                    continue
                for iz in range(nearz - 1, nearz + 2):
                    if iz < 0 or nz <= iz:
                        continue
                    h2.SetBinContent(iR + 1, iz + 1, data[i][iz][iR])
                        
            F[i] = h2.Interpolate(R, z)/E[i]**2

        spec = DiffuseSpectrum(par, E, dEl, dEh, F)
        
        return spec

    def Mori1997(self, par, E, dEl, dEh):
        """
        M. Mori, The Astrophysical Journal 478 (1997) 225-232
        """
        if not isinstance(par, matter.Nucleus):
            raise TypeError, "Invalid particle type"
        if not (E.size == dEl.size == dEh.size):
            raise TypeError, "Size of E/dEl/dEh are different"
        
        F = numpy.zeros(E.size)

        m = matter.proton.mass
        if par != matter.proton:
            ab = matter.gaisser_table.cr_abundance(par)
        
        for i in range(E.size):
            # see eq (3)
            Tp = E[i]
            Ep = Tp + m
            if Ep <= 1e5: # <= 100[GeV] NOTE! Total energy
                p = (Ep**2 - m**2)**0.5 # [MeV/c]
                F[i] = 1.67*(p/1e3)**-2.7*(1 + (p/1e3/2.5)**-2.)**-0.5
            else:
                F[i] = 6.65e-6*(Ep/100e3)**-2.75
            
            if par != matter.proton:
                # see eq (5)
                if Tp < 1e5: # <= 100 [GeV/n] NOTE! Kinetic energy per nucleon
                    F[i] *= ab
                else:
                    F[i] *= ab*(Tp/1e5)**0.12
                
            F[i] *= 1e-3 # [/cm^2/s/sr/GeV] -> [/cm^2/s/sr/MeV]

        spec = DiffuseSpectrum(par, E, dEl, dEh, F)
                    
        return spec
    
    def Honda2004(self, par, E, dEl, dEh):
        """
        M. Honda et al. Physical Review D 70 (2004) 043008
        """
        if not isinstance(par, matter.Nucleus):
            raise TypeError, "Invalid particle type"
        if not (E.size == dEl.size == dEh.size):
            raise TypeError, "Size of E/dEl/dEh are different"
        
        F = numpy.zeros(E.size)

        if par == matter.proton:
            a, K, b, c = 2.74, 14900, 2.15, 0.21
        elif par == matter.alpha:
            a, K, b, c = 2.64, 600, 1.25, 0.14
        else:
            a, K, b, c = 0, 0, 0, 0
        
        for i in range(E.size):
            F[i] = K*(E[i]/1e3 + b*math.exp(-c*(E[i]/1e3)**0.5))**(-a)

        # [/m^2/sr/s/GeV] -> [/cm^2/sr/s/(MeV/n)] (x 1e-7)
        F *= 1e-7
        
        spec = DiffuseSpectrum(par, E, dEl, dEh, F)
                    
        return spec
        
    def ATIC2(self, par):
        """
        A. D. Panov, et al.,
        Bulletin of the Russian Academy of Sciences: Physics 71 (2007) 494-497 
        and private communication with Dr. Panov
        """
        if par == matter.proton:
            fname = pkg_resources.resource_filename("cr_flux", "data/atic/panov2007_proton.dat")
        elif par == matter.alpha:
            fname = pkg_resources.resource_filename("cr_flux", "data/atic/panov2007_alpha.dat")
        elif par == matter.electron:
            fname = pkg_resources.resource_filename("cr_flux", "data/atic/chang2008_electron.dat")
        else:
            raise TypeError, "Invalid particle type"

        f = open(fname)
        lines = f.readlines()[2:] # skip the header
        E   = numpy.zeros(len(lines))
        dEh = numpy.zeros(len(lines))
        dEl = numpy.zeros(len(lines))
        F   = numpy.zeros(len(lines))
        dFh = numpy.zeros(len(lines))
        dFl = numpy.zeros(len(lines))
               
        if par == matter.proton:
            for i, line in enumerate(lines):
                val = [float(x) for x in line.split()]
                # [GeV] -> [MeV/n] (x 1e3)
                # [/m^2/sr/s/GeV] -> [/cm^2/sr/s/(MeV/n)] (x 1e-7)
                E[i]   = val[2]*1e3
                dEl[i] = E[i] - val[0]*1e3
                dEh[i] = val[1]*1e3 - E[i]
                F[i]   = val[3]*1e-7
                dFl[i] = dFh[i] = val[4]*1e-7
        elif par == matter.alpha:
            for i, line in enumerate(lines):
                val = [float(x) for x in line.split()]
                # [GeV] -> [MeV/n] (x 1e3 / A)
                # [/m^2/s/sr/GeV] -> [/cm^2/sr/s/(MeV/n)] (x A x 1e-7)
                E[i] = val[2]*1e3/par.A # [GeV] -> [MeV/n]
                dEl[i] = E[i] - val[0]*1e3/par.A
                dEh[i] = val[1]*1e3/par.A - E[i]
                F[i]   = val[3]*par.A*1e-7
                dFl[i] = dFh[i] = val[4]*par.A*1e-7
        elif par == matter.electron:
            for i, line in enumerate(lines):
                val = [float(x) for x in line.split()]
                # [GeV] -> [MeV] (x 1e3)
                # [GeV^3/m^2/s/sr/GeV] -> [/cm^2/sr/s/MeV] (x 1e-7)
                E[i] = val[2]
                dEl[i] = E[i] - val[0]
                dEh[i] = val[1] - E[i]
                F[i]   = val[4]
                dFh[i] = val[5] - F[i]
                dFl[i] = F[i] - val[6]
                F[i]   *= 1e-7/E[i]**3
                dFl[i] *= 1e-7/E[i]**3
                dFh[i] *= 1e-7/E[i]**3
                E[i]   *= 1e3
                dEl[i] *= 1e3
                dEh[i] *= 1e3
            
        spec = DiffuseSpectrum(par, E, dEl, dEh, F, dFl, dFh)
                    
        return spec
    
    def AMS01(self, par):
        """
        Create spectra of proton and alpha from AMS result.
        J. Alcaraz, et al., Physics Letters B 490 (2000) 27-35
        J. Alcaraz, et al., Physics Letters B 494 (2000) 193-202
        J. Alcaraz, et al., Physics Letters B 484 (2000) 10-22 
        M. Aguilar, et al., Physics Reports 366 (2002) 331-405
        """
        if par == matter.proton:
            fname = pkg_resources.resource_filename("cr_flux", "data/ams/alcaraz2000_proton.dat")
        elif par == matter.alpha:
            fname = pkg_resources.resource_filename("cr_flux", "data/ams/alcaraz2000_alpha.dat")
        elif par == matter.deutron:
            fname = pkg_resources.resource_filename("cr_flux", "data/ams/aguilar2002_deutron.dat")
        elif par == matter.electron:
            fname = pkg_resources.resource_filename("cr_flux", "data/ams/alcaraz2000_electron.dat")
        elif par == matter.positron:
            fname = pkg_resources.resource_filename("cr_flux", "data/ams/alcaraz2000_positron.dat")
        else:
            raise TypeError, "Invalid particle type"

        f = open(fname)
        lines = f.readlines()[1:] # skip the header
        E   = numpy.zeros(len(lines))
        dEh = numpy.zeros(len(lines))
        dEl = numpy.zeros(len(lines))
        F   = numpy.zeros(len(lines))
        dFh = numpy.zeros(len(lines))
        dFl = numpy.zeros(len(lines))
        
        if par == matter.proton:
            for i, line in enumerate(lines):
                val = [float(x) for x in line.split()]
                # [GeV] -> [MeV/n] (x 1e3)
                # [/m^2/sr/s/(MeV/n)] -> [/cm^2/sr/s/(MeV/n)] (x 1e-4)
                val[0] *= 1e3 # [GeV] -> [MeV]
                val[1] *= 1e3 # [GeV] -> [MeV]
                E[i]   = (val[1] * val[0])**0.5
                dEl[i] = E[i] - val[0]
                dEh[i] = val[1] - E[i]
                F[i]   = val[2]*val[7]*1e-4
                dFl[i] = dFh[i] = (val[3]**2 + val[4]**2 + val[5]**2 + val[6]**2)**0.5*val[7]*1e-4
        elif par == matter.alpha:
            for i, line in enumerate(lines):
                val = [float(x) for x in line.split()]
                # [GV] -> [MeV/n] (x Z / A x 1e3)
                # [/m^2/s/sr/GV] -> [/cm^2/sr/s/(MeV/n)] (x A / Z x 1e-7)
                A, Z, m = float(par.A), float(par.Z), par.mass
                R1, R2 = val[0], val[1]
                dNdR = val[2] # dN/dR
                dN = dNdR*(R2 - R1)
                ddN = val[3]*(R2 - R1)
                # [MeV/n]
                E1 = ((m**2 + (R1*Z*1e3)**2)**0.5 - m)/A 
                E2 = ((m**2 + (R2*Z*1e3)**2)**0.5 - m)/A 
                dE = E2 - E1
                dNdE = dN/dE
                ddNdE = ddN/dE
                E[i]   = (E1 * E2)**0.5
                dEl[i] = E[i] - E1
                dEh[i] = E2 - E[i]
                F[i]   = dNdE*1e-4 # [m^2] -> [cm^2]
                dFl[i] = dFh[i] = ddNdE*1e-4
        elif par == matter.deutron:
            for i, line in enumerate(lines):
                val = [float(x) for x in line.split()]
                # [GeV/n] -> [MeV/n] (x 1e3)
                # [/m^2/sr/s/(MeV/n)] -> [/cm^2/sr/s/(MeV/n)] (x 1e-4)
                val[0] *= 1e3 # [GeV/n] -> [MeV/n]
                val[1] *= 1e3 # [GeV/n] -> [MeV/n]
                E[i]   = (val[1] * val[0])**0.5
                dEl[i] = E[i] - val[0]
                dEh[i] = val[1] - E[i]
                F[i]   = val[2]*val[5]*1e-4 # [/cm^2/sr/s/(MeV/n)]
                dFl[i] = dFh[i] = (val[3]**2 + val[4]**2)**0.5*val[5]*1e-4
        elif par == matter.electron or par == matter.positron:
            for i, line in enumerate(lines):
                val = [float(x) for x in line.split()]
                # [GeV/n] -> [MeV/n] (x 1e3)
                # [/m^2/sr/s/(GeV/n)] -> [/cm^2/sr/s/(MeV/n)] (x 1e-7)
                val[0] *= 1e3 # [GeV/n] -> [MeV/n]
                val[1] *= 1e3 # [GeV/n] -> [MeV/n]
                E[i]   = (val[1] * val[0])**0.5
                dEl[i] = E[i] - val[0]
                dEh[i] = val[1] - E[i]
                F[i]   = val[2]*1e-7 # [/cm^2/sr/s/(MeV/n)]
                dFl[i] = dFh[i] = val[3]*1e-7
        
        spec = DiffuseSpectrum(par, E, dEl, dEh, F, dFl, dFh)
                    
        return spec
        
    def AMS02(self, par):
        """
        L. Accardo, et al., Physical Review Letters 113 (2014) 121101
        M. Aguilar, et al., Physical Review Letters 113 (2014) 121102
        """
        if par == matter.electron:
            fname = pkg_resources.resource_filename("cr_flux", "data/ams/aguilar2014_electron.dat")
        elif par == matter.positron:
            fname = pkg_resources.resource_filename("cr_flux", "data/ams/aguilar2014_positron.dat")
        else:
            raise TypeError, "Invalid particle type"

        f = open(fname)
        lines = f.readlines()[1:] # skip the header
        E   = numpy.zeros(len(lines))
        dEh = numpy.zeros(len(lines))
        dEl = numpy.zeros(len(lines))
        F   = numpy.zeros(len(lines))
        dFh = numpy.zeros(len(lines))
        dFl = numpy.zeros(len(lines))

        if par == matter.electron or par == matter.positron:
            for i, line in enumerate(lines):
                val = [float(x) for x in line.split()]
                # [GeV/n] -> [MeV/n] (x 1e3)
                # [/m^2/sr/s/(GeV/n)] -> [/cm^2/sr/s/(MeV/n)] (x 1e-7)
                val[0] *= 1e3 # [GeV/n] -> [MeV/n]
                val[1] *= 1e3 # [GeV/n] -> [MeV/n]
                E[i]   = (val[1] * val[0])**0.5
                dEl[i] = E[i] - val[0]
                dEh[i] = val[1] - E[i]
                F[i]   = val[2]*1e-7 # [/cm^2/sr/s/(MeV/n)]
                dFl[i] = dFh[i] = ((val[3]**2 + val[4]**2)**0.5)*1e-7
        
        spec = DiffuseSpectrum(par, E, dEl, dEh, F, dFl, dFh)
                    
        return spec

    def BESS_TeV(self, par):
        """
        Create spectra of proton and alpha from BESS-TeV result.
        S. Haino, et al., Physics Letters B 594 (2004) 35-46        
        """
        if par == matter.proton:
            fname = pkg_resources.resource_filename("cr_flux", "data/bess/haino2004_proton.dat")
        elif par == matter.alpha:
            fname = pkg_resources.resource_filename("cr_flux", "data/bess/haino2004_alpha.dat")
        else:
            raise TypeError, "Invalid particle type"
        
        f = open(fname)
        lines = f.readlines()[1:] # skip the header
        E   = numpy.zeros(len(lines))
        dEh = numpy.zeros(len(lines))
        dEl = numpy.zeros(len(lines))
        F   = numpy.zeros(len(lines))
        dFh = numpy.zeros(len(lines))
        dFl = numpy.zeros(len(lines))
        
        for i, line in enumerate(lines):
            val = [float(x) for x in line.split()]
            # [GeV/n] -> [MeV/n] (x 1e3)
            # [/m^2/sr/s/(GeV/n)] -> [/cm^2/sr/s/(MeV/n)] (x 1e-7)
            E[i], dEl[i], dEh[i], F[i], dFl[i] =\
            val[2]*1e3, (val[2] - val[0])*1e3, (val[1] - val[2])*1e3,\
            val[3]*1e-7, (val[4]**2 + val[5]**2)**0.5*1e-7
            dFh[i] = dFl[i]
            
        spec = DiffuseSpectrum(par, E, dEl, dEh, F, dFl, dFh)
                    
        return spec

    def BESS_LIS(self, par, E, dEl, dEh):
        """
        Create spectra of proton and alpha from the long-term result of BESS.
        These spectra are "de-modulated" from observed ones.
        Y. Shikaze, et al., Astroparticle Physics 28 (2007) 154�167
        """
        if not (par == matter.proton or par == matter.alpha):
            raise TypeError, "Invalid particle type"
        
        def lis(spec, A, P1, P2):
            m = spec.par.mass
            for i, e in enumerate(spec.E):
                e *= spec.par.A # per particle
                g = (e + m)/m # gamma factor
                b = (1. - 1./g**2)**0.5
                p = ((e + m)**2 - m**2)**0.5
                spec.F[i] = (A/1e7)*b**P1*(p/spec.par.A*spec.par.Z/1e3)**-P2
                spec.dFl[i] = spec.dFh[i] = 0.

        F = numpy.zeros(E.size)
        spec = DiffuseSpectrum(par, E, dEl, dEh, F)
        
        if spec.par == matter.proton:
            lis(spec, 1.94e4, 0.70, 2.76)
        elif spec.par == matter.alpha:
            lis(spec, 7.10e3, 0.50, 2.78)

        spec.init_graph()
                    
        return spec

    def CAPRICE94(self, par):
        """
        Create spectra of electron and positron from CAPRICE94 result.
        M. Boezio, et al., The Astrophysical Journal 532 (2000) 653-669
        """
        if par == matter.electron:
            fname = pkg_resources.resource_filename("cr_flux", "data/caprice/boezio2000_electron.dat")
        elif par == matter.positron:
            fname = pkg_resources.resource_filename("cr_flux", "data/caprice/boezio2000_positron.dat")
        else:
            raise TypeError, "Invalid particle type"
        
        f = open(fname)
        lines = f.readlines()[2:] # skip the header
        E   = numpy.zeros(len(lines))
        dEh = numpy.zeros(len(lines))
        dEl = numpy.zeros(len(lines))
        F   = numpy.zeros(len(lines))
        dFh = numpy.zeros(len(lines))
        dFl = numpy.zeros(len(lines))
        
        for i, line in enumerate(lines):
            val = [float(x) for x in line.split()]
            # [GeV] -> [MeV] (x 1e3)
            # [/m^2/sr/s/GeV] -> [/cm^2/sr/s/MeV] (x 1e-7)
            E[i], dEl[i], dEh[i], F[i], dFl[i], dFh[i] =\
            val[2]*1e3, (val[2] - val[0])*1e3, (val[1] - val[2])*1e3,\
            val[3]*1e-7*10**val[6], val[4]*1e-7*10**val[6], val[5]*1e-7*10**val[6]
            
        spec = DiffuseSpectrum(par, E, dEl, dEh, F, dFl, dFh)
                    
        return spec

    def Fermi(self, par, year=2009):
        """
        Create spectra of electron (+positron) from Fermi/LAT result.
        A. A. Abdo, et al., Physical Review Letters 102 (2009) 181101
        """
        if par == matter.electron:
            if year == 2009:
                fname = pkg_resources.resource_filename("cr_flux", "data/fermi/abdo2009_electron.dat")
            elif year == 2012:
                fname = pkg_resources.resource_filename("cr_flux", "data/fermi/ackermann2012_electron.dat")
            else:
                raise TypeError, "Invalid year"
        elif par == matter.positron:
            fname = pkg_resources.resource_filename("cr_flux", "data/fermi/ackermann2012_positron.dat")
        else:
            raise TypeError, "Invalid particle type"
        
        f = open(fname)
        lines = f.readlines()[2:] # skip the header
        E   = numpy.zeros(len(lines))
        dEh = numpy.zeros(len(lines))
        dEl = numpy.zeros(len(lines))
        F   = numpy.zeros(len(lines))
        dFh = numpy.zeros(len(lines))
        dFl = numpy.zeros(len(lines))
        
        for i, line in enumerate(lines):
            val = [float(x) for x in line.split()]
            # [GeV] -> [MeV] (x 1e3)
            E[i] = (val[0] + val[1])*0.5*1e3
            dEl[i], dEh[i] = E[i] - val[0]*1e3, val[1]*1e3 - E[i]
            if year == 2009:
                # [GeV^3/m^2/sr/s/GeV] -> [/cm^2/sr/s/MeV]
                for j in range(2, 6):
                    val[j] = val[j]/(E[i]/1e3)**3*1e-7
            elif year == 2012:
                # [/m^2/sr/s/GeV] -> [/cm^2/sr/s/MeV]
                for j in range(2, 6):
                    val[j] = val[j]*1e-7

            F[i] = val[2]
            dFh[i] = (val[3]**2 + val[4]**2)**0.5
            dFl[i] = (val[3]**2 + val[5]**2)**0.5

        spec = DiffuseSpectrum(par, E, dEl, dEh, F, dFl, dFh)
                    
        return spec

    def HEAT(self, par):
        """
        Create spectra from HEAT result.
        S. W. Barwick, et al., The Astrophysical Journal 498 (1998) 779-789
        """
        if par == matter.electron:
            fname = pkg_resources.resource_filename("cr_flux", "data/heat/barwick1998_electron.dat")
        elif par == matter.positron:
            fname = pkg_resources.resource_filename("cr_flux", "data/heat/barwick1998_positron.dat")
        else:
            raise TypeError, "Invalid particle type"
        
        f = open(fname)
        lines = f.readlines()[1:] # skip the header
        E   = numpy.zeros(len(lines))
        dEh = numpy.zeros(len(lines))
        dEl = numpy.zeros(len(lines))
        F   = numpy.zeros(len(lines))
        dFh = numpy.zeros(len(lines))
        dFl = numpy.zeros(len(lines))
        
        for i, line in enumerate(lines):
            val = [float(x) for x in line.split()]
            # [GeV/n] -> [MeV/n] (x 1e3)
            # [/m^2/sr/s/(GeV/n)] -> [/cm^2/sr/s/(MeV/n)] (x 1e-7)
            E[i], dEl[i], dEh[i], F[i], dFh[i], dFl[i] =\
            val[2]*1e3, (val[2] - val[0])*1e3, (val[1] - val[2])*1e3,\
            val[3]*1e-7, val[4]*1e-7, val[5]*1e-7
            
        spec = DiffuseSpectrum(par, E, dEl, dEh, F, dFl, dFh)
                    
        return spec

    def HESS(self, par, year):
        """
        Create spectra from HESS result.
        F. A. Aharonian, et al., Physical Review Letters 101 (2008) 261104
        F. A. Aharonian, et al., A&A 508 (2009) 561-564 
        """
        if par == matter.electron:
            if year == 2008:
                fname1 = pkg_resources.resource_filename("cr_flux", "data/hess/aharonian2008_electron.dat")
                fname2 = pkg_resources.resource_filename("cr_flux", "data/hess/aharonian2008_electron_sys.dat")
            elif year == 2009:
                fname1 = pkg_resources.resource_filename("cr_flux", "data/hess/aharonian2009_electron.dat")
                fname2 = pkg_resources.resource_filename("cr_flux", "data/hess/aharonian2009_electron_sys.dat")
            else:
                raise TypeError, "Invalid year"
        else:
            raise TypeError, "Invalid particle type"
        
        gra1 = ROOT.TGraph(0) # upper systematic error 
        gra2 = ROOT.TGraph(0) # lower systematic error 
        
        f = open(fname2)
        lines = f.readlines()[2:] # skip the header

        for i, line in enumerate(lines):
            # [GeV^2/m^2/sr/s] -> [/m^2/sr/s/GeV]
            val = [float(x) for x in line.split()]
            gra1.SetPoint(i, val[0], val[1]/val[0]**3)
            gra2.SetPoint(i, val[0], val[2]/val[0]**3)
        
        f = open(fname1)
        lines = f.readlines()[2:] # skip the header
        E   = numpy.zeros(len(lines))
        dEh = numpy.zeros(len(lines))
        dEl = numpy.zeros(len(lines))
        F   = numpy.zeros(len(lines))
        dFh = numpy.zeros(len(lines))
        dFl = numpy.zeros(len(lines))
        
        for i, line in enumerate(lines):
            val = [float(x) for x in line.split()]
            # [GeV] -> [MeV] (x 1e3)
            # [/m^2/sr/s/GeV] -> [/cm^2/sr/s/MeV] (x 1e-7)
            E[i] = val[0]*1e3
            F[i] = val[1]/val[0]**3*1e-7
            dFl[i] = 1e-7*((val[2]/val[0]**3)**2 + (val[1]/val[0]**3 - gra2.Eval(val[0]))**2)**0.5
            dFh[i] = 1e-7*((val[3]/val[0]**3)**2 + (gra1.Eval(val[0]) - val[1]/val[0]**3)**2)**0.5
            
        spec = DiffuseSpectrum(par, E, dEl, dEh, F, dFl, dFh)
                    
        return spec

    def IMAX(self, par):
        """
        Create spectra of proton and alpha from IMAX result.
        W. Menn, et al., The Astrophysical Journal 533 (2000) 281-297
        """
        if par == matter.proton:
            fname = pkg_resources.resource_filename("cr_flux", "data/imax/menn2000_proton.dat")
        elif par == matter.alpha:
            fname = pkg_resources.resource_filename("cr_flux", "data/imax/menn2000_alpha.dat")
        else:
            raise TypeError, "Invalid particle type"
        
        f = open(fname)
        lines = f.readlines()[2:] # skip the header
        E   = numpy.zeros(len(lines))
        dEh = numpy.zeros(len(lines))
        dEl = numpy.zeros(len(lines))
        F   = numpy.zeros(len(lines))
        dFh = numpy.zeros(len(lines))
        dFl = numpy.zeros(len(lines))
        
        for i, line in enumerate(lines):
            val = [float(x) for x in line.split()]
            # [GeV/n] -> [MeV/n] (x 1e3)
            # [/m^2/sr/s/(GeV/n)] -> [/cm^2/sr/s/(MeV/n)] (x 1e-7)
            E[i], dEl[i], dEh[i], F[i], dFl[i] =\
            val[2]*1e3, (val[2] - val[0])*1e3, (val[1] - val[2])*1e3,\
            val[3]*1e-7*10**val[6], (val[4]**2 + val[5]**2)**0.5*1e-7*10**val[6]
            dFh[i] = dFl[i]
            
        spec = DiffuseSpectrum(par, E, dEl, dEh, F, dFl, dFh)
                    
        return spec

    def MASS(self, par):
        """
        Create spectrum of proton from MASS result.
        R. Bellotti, et al. Phsyical Review D 60 (1999) 052002
        """
        if par == matter.proton:
            fname = pkg_resources.resource_filename("cr_flux", "data/mass/bellotti1999_proton.dat")
        else:
            raise TypeError, "Invalid particle type"
        
        f = open(fname)
        lines = f.readlines()[2:] # skip the header
        E   = numpy.zeros(len(lines))
        dEh = numpy.zeros(len(lines))
        dEl = numpy.zeros(len(lines))
        F   = numpy.zeros(len(lines))
        dFh = numpy.zeros(len(lines))
        dFl = numpy.zeros(len(lines))
        
        for i, line in enumerate(lines):
            val = [float(x) for x in line.split()]
            # [GeV/n] -> [MeV/n] (x 1e3)
            # [/m^2/sr/s/GeV] -> [/cm^2/sr/s/MeV] (x 1e-7)
            E1 = 1e3*val[0]
            E2 = 1e3*val[1]
            E[i] = 1e3*val[2]
            dEl[i] = E[i] - E1
            dEh[i] = E2 - E[i]
            F[i] = val[3]*1e-7
            dFh[i] = dFl[i] = val[4]*1e-7
            
        spec = DiffuseSpectrum(par, E, dEl, dEh, F, dFl, dFh)
                    
        return spec

    def PAMELA(self, par, year):
        """
        Create spectra of electron (+positron) from Fermi/LAT result.
        A. A. Abdo, et al., Physical Review Letters 102 (2009) 181101
        """
        if par == matter.electron:
            if year == 2011:
                fname = pkg_resources.resource_filename("cr_flux", "data/pamela/adriani2011_electron.dat")
            else:
                raise TypeError, "Invalid year"
        elif par == matter.positron:
            if year == 2013:
                fname = pkg_resources.resource_filename("cr_flux", "data/pamela/adriani2013_positron.dat")
            else:
                raise TypeError, "Invalid year"
        else:
            raise TypeError, "Invalid particle type"
        
        f = open(fname)
        lines = f.readlines()[2:] # skip the header
        E   = numpy.zeros(len(lines))
        dEh = numpy.zeros(len(lines))
        dEl = numpy.zeros(len(lines))
        F   = numpy.zeros(len(lines))
        dFh = numpy.zeros(len(lines))
        dFl = numpy.zeros(len(lines))
        
        for i, line in enumerate(lines):
            val = [float(x) for x in line.split()]
            # [GeV] -> [MeV] (x 1e3)
            E[i] = (val[0] + val[1])*0.5*1e3
            dEl[i], dEh[i] = E[i] - val[0]*1e3, val[1]*1e3 - E[i]
            # [/m^2/sr/s/GeV] -> [/cm^2/sr/s/MeV]
            if year == 2013:
                for j in range(2, 5):
                    val[j] = val[j]*1e-7
            elif year == 2011:
                for j in range(2, 6):
                    val[j] = val[j]*1e-7

            F[i] = val[2]
            if year == 2011:
                dFh[i] = (val[3]**2 + val[4]**2)**0.5
                dFl[i] = (val[4]**2 + val[5]**2)**0.5
            elif year == 2013:
                dFh[i] = (val[3]**2 + val[4]**2)**0.5
                dFl[i] = dFh[i]

        spec = DiffuseSpectrum(par, E, dEl, dEh, F, dFl, dFh)
                    
        return spec
    
    def RUNJOB(self, par):
        """
        Create spectra of proton and alpha from RUNJOB result.
        A.V. Apanasenko et al. Astroparticle Physics 16 (2001) 13-46
        """
        if par == matter.proton:
            fname = pkg_resources.resource_filename("cr_flux", "data/runjob/apanasenko2001_proton.dat")
        elif par == matter.alpha:
            fname = pkg_resources.resource_filename("cr_flux", "data/runjob/apanasenko2001_alpha.dat")
        else:
            raise TypeError, "Invalid particle type"
        
        f = open(fname)
        lines = f.readlines()[1:] # skip the header
        E   = numpy.zeros(len(lines))
        dEh = numpy.zeros(len(lines))
        dEl = numpy.zeros(len(lines))
        F   = numpy.zeros(len(lines))
        dFh = numpy.zeros(len(lines))
        dFl = numpy.zeros(len(lines))
        
        for i, line in enumerate(lines):
            val = [float(x) for x in line.split()]
            # [GeV/n] -> [MeV/n] (x 1e3)
            # [/m^2/sr/s/(GeV/n)] -> [/cm^2/sr/s/(MeV/n)] (x 1e-7)
            E1 = 1e3*val[0]*10**val[2]
            E2 = 1e3*val[1]*10**val[2]
            E[i] = (E1*E2)**0.5
            dEl[i] = E[i] - E1
            dEh[i] = E2 - E[i]
            F[i] = val[3]*1e-7*10**val[6]
            dFh[i] = val[4]*1e-7*10**val[6]
            dFl[i] = -val[5]*1e-7*10**val[6]
            
        spec = DiffuseSpectrum(par, E, dEl, dEh, F, dFl, dFh)
                    
        return spec

def create_log_binned_energies(emin, emax, nbins):
    """
    Create binned energies (E, dEl, dEh)
    emin: Minimum energy [MeV]
    emax: Maximum energy [MeV]
    nbins: Number of bins between emin and emax
    """
    step = (float(emax)/float(emin))**(1./nbins)
    E   = emin*step**0.5*step**numpy.arange(nbins)
    dEl = E - E/step**0.5
    dEh = E*step**0.5 - E

    return E, dEl, dEh