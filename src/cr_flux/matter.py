"""
Definition of particles and interstellar medium
$Id: matter.py,v 1.2 2009/02/15 14:37:42 oxon Exp $
"""

import pkg_resources

import spectrum

class Particle(object):
    """
    Base class for particles
    """
    def __init__(self):
        raise NotImplementedError

class ChargedLepton(Particle):
    """
    Class for charged leptons
    """
    def __init__(self, charge, generation):
        """
        charge: +1/-1
        generation: 0, 1, 2 (e/mu/tau)
        
        e.g. electron = (-1, 0)
             anti-tauon = (+1, 2)             
        """
        if not (charge == -1 or charge == 1):
            raise TypeError, "Charge must be +1, 0 or -1"
        if not (generation == 0 or generation == 1 or generation == 2):
            raise TypeError, "Generation must be 0, 1 or 2"
        
        self.charge = charge
        self.generation = generation
        
    def __cmp__(self, other):
        return cmp([self.charge, self.generation], [other.charge, other.generation])

    def __hash__(self):
        if self.charge > 0:
            return 100*(self.charge*1000 + self.generation) + 1
        else:
            return 100*(-self.charge*1000 + self.generation) + 2

class Nucleus(Particle):
    """
    Class for nuclei. This class is also used to define proton, neutron and
    anti nuclei.
    """
    def __init__(self, A, Z):
        """
        A: Mass number of the nucleus (A > 0)
        Z: Charge of the nucleus
        
        e.g. proton      Nucleus(1, 1)
             neutron     Nucleus(1, 0)
             anti-proton Nucleus(1, -1)
             alpha       Nucleus(4, 2)
        """
        if A < 1:
            raise TypeError, "Mass number must be greater than zero"
        if Z == 0 and A != 1:
            raise TypeError, "Only neutron can be used as non-charged"
        if A < abs(Z):
            raise TypeError, "A must be greater than or equals to Z"
        
        self.A = A
        self.Z = Z
        # Rough estimation of the mass
        # Mass defect (~1%) is not taken into account 
        self.mass = 938.27203*abs(Z) + 939.56536*(A - abs(Z)) # [MeV]

    def __cmp__(self, other):
        return cmp([self.Z, self.A], [other.Z, other.A])
    
    def __hash__(self):
        if self.Z > 0:
            return 100*(self.A*1000 + self.Z) + 3
        else:
            return 100*(self.A*1000 - self.Z) + 4

class Photon(Particle):
    """
    Class for photon
    """
    def __init__(self):
        self.mass = 0.
    
    def __hash__(self):
        return 5

class CRAbundance(object):
    """
    Class for cosmic-ray abundance. The abundance of proton is always normalized
    to 1.
    """
    def __init__(self):
        """
        The abundance of proton is set to 1.
        """
        self.abtable = {}
        self.abtable[proton] = 1.0

    def cr_abundance(self, par):
        """
        Returns the abundance of par.
        """
        try:
            ab = self.abtable[par]
        except KeyError:
            ab = 0.

        return ab

    def add_cr(self, par, ab):
        """
        Add a new particle type and its abundance.
        """
        self.abtable[par] = ab
        
    def expand_flux_list(self, f_list):
        """
        Expand flux list assuming CR abundance. cr_list[-1] will be the shape
        of other cosmic-ray flux.
        f_list: list of cosmic-ray flux
        """
        par_list = []
        for f in f_list:
            par_list.append(f.par)

        base_spec = f_list[-1]
        
        for par in self.abtable.keys():
            try:
                par_list.index(par)
            except ValueError:
                norm = self.cr_abundance(par)/self.cr_abundance(base_spec.par)
                spec = spectrum.DiffuseSpectrum(par, base_spec.E, base_spec.dEl,
                                                base_spec.dEh, base_spec.F*norm,
                                                base_spec.dFl*norm,
                                                base_spec.dFh*norm)
                f_list.append(spec)

class ISMAbundance(object):
    """
    Class for interstellar medium abundance. The abundance of hydrogen is always
    normalized to 1.
    """
    def __init__(self):
        """
        The abundance of hydrogen is set to 1.
        """
        self.abtable = {}
        self.abtable[proton] = 1.0
    
    def ism_abundance(self, par):
        """
        Returns the abundance of par.
        """
        try:
            ab = self.abtable[par]
        except KeyError:
            ab = 0.
            
        return ab
    
    def add_ism(self, par, ab):
        """
        Add a new particle type and its abundance.
        """
        self.abtable[par] = ab

class Multiplicity(object):
    """
    Class for multiplicity of nucl-nucl interactions. Cross section of p-p
    collision is normalized to 1. This idea is based on
    T. K. Gaisser and R. K. Schaefer
    The Astrophysical Journal (1992) 394, 174-183
    """
    def __init__(self):
        """
        Set the multiplicity of p-p collision to 1. If self.force == True,
        then multiplicities of not-registered collision is the product of
        mass numbers of particles.
        """
        self.multable = {}
        self.multable[(proton, proton)] = 1.
        self.force = True
        
    def multiplicity(self, par1, par2):
        """
        par1: Type of particle in collision
        par2: Type of particle in collision
        (par1, par2) and (par2, par1) are the same meaning.
        """
        try:
            m = self.multable[(par1, par2)]
        except KeyError:
            try:
                # search an invsersed pair. physically same meaning.
                m = self.multable[(par2, par1)]
            except:
                if self.force:
                    m = par1.A * par2.A # ignore shadowing
                else:
                    m = 0.
            
        return m
    
    def add_mul(self, par1, par2, mul):
        """
        Add a new particle pair and its multiplicity.
        """
        self.multable[(par1, par2)] = mul
        
class GaisserTable(CRAbundance, Multiplicity):
    """
    See the table 1 of
    T. K. Gaisser and R. K. Schaefer
    The Astrophysical Journal (1992) 394, 174-183
    """
    __instance = None
    def __new__(cls, *args, **kwargs):
        """
        Singleton
        """
        if cls != type(cls.__instance):
            cls.__instance = object.__new__(cls, *args, **kwargs)
            
        return cls.__instance

    def __init__(self):
        CRAbundance.__init__(self)
        Multiplicity.__init__(self)
        self.read_table()

    def read_table(self):
        """
        Read the data of table 1 from ASCII.
        """
        fname = pkg_resources.resource_filename("cr_flux", "data/abundance/gaisser1992.dat")

        f = open(fname)
        lines = f.readlines()[1:] # skip the header
        
        for line in lines:
            val = line.split()
            A = int(val[1])
            Z = int(val[2])
            # 0:Name 1:A 2:Z 3:ni/nH|CR 4:mip 5:eH 6:mia 7:eHe
            self.add_cr(Nucleus(A, Z), float(val[4]))
            self.add_mul(Nucleus(1, 1), Nucleus(A, Z), float(val[5]))
            self.add_mul(Nucleus(4, 2), Nucleus(A, Z), float(val[7]))

class Gas(object):
    """
    Base class for gas
    """
    def __init__(self, cr_ab, ism_ab):
        """
        cr_ab: Cosmic-ray abundance in the gas
        ism_ab: ISM abundance in the gas
        """
        self.cr_ab  = cr_ab
        self.ism_ab = ism_ab

class GasDensity(Gas):
    """
    Class for gas density. It is used when calculating gamma-ray emissivity
    par unit density.
    """
    def __init__(self, cr_ab, ism_ab, density):
        """
        density: n(H) [/cm^3]
        """
        Gas.__init__(self, cr_ab, ism_ab)
        self.density = density
    
class Cloud(Gas):
    """
    Class for molecular clouds.
    """
    def __init__(self, cr_ab, ism_ab, mass, dist):
        """
        mass: Total mass  of the cloud in unit of solar mass
        dist: Distance between the Sun and the Cloud [pc]
        """
        Gas.__init__(self, cr_ab, ism_ab)
        self.mass = mass
        self.dist = dist

proton = Nucleus(1, 1)
alpha = Nucleus(4, 2)
gaisser_table = GaisserTable()