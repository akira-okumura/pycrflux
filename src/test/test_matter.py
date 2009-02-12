"""
Unit test script for matter.py
$Id: test_matter.py,v 1.1 2009/02/12 18:12:35 oxon Exp $
"""

import numpy
import unittest

import cr_flux

class TestParticle(unittest.TestCase):
    """
    Unit test for class Particle
    """
    def setUp(self):
        self.p1 = cr_flux.Particle()
        self.p2 = cr_flux.Particle()
        
    def testCmp(self):
        self.assertEqual(self.p1, self.p2)
        
    def testHash(self):
        dic = {}
        dic[self.p1] = 1
        self.assertEqual(dic[self.p2], 1)

class TestChargedLepton(unittest.TestCase):        
    """
    Unit test for class ChargedLepton
    """
    def setUp(self):
        self.e    = cr_flux.ChargedLepton(-1, 0)
        self.ebar = cr_flux.ChargedLepton(+1, 0)
        self.mass    = cr_flux.ChargedLepton(-1, 1)
        self.mbar = cr_flux.ChargedLepton(+1, 1)
        self.t    = cr_flux.ChargedLepton(-1, 2)
        self.tbar = cr_flux.ChargedLepton(+1, 2)
        
    def testCmp(self):
        self.assertEqual(self.e,    cr_flux.ChargedLepton(-1, 0))
        self.assertEqual(self.ebar, cr_flux.ChargedLepton(+1, 0))
        self.assertEqual(self.mass,    cr_flux.ChargedLepton(-1, 1))
        self.assertEqual(self.mbar, cr_flux.ChargedLepton(+1, 1))
        self.assertEqual(self.t,    cr_flux.ChargedLepton(-1, 2))
        self.assertEqual(self.tbar, cr_flux.ChargedLepton(+1, 2))
        
    def testHash(self):
        dic = {}
        dic[self.e] = 1
        self.assertEqual(dic[cr_flux.ChargedLepton(-1, 0)], 1)

class TestNucleus(unittest.TestCase):        
    """
    Unit test for class Nucleus
    """
    def setUp(self):
        self.p    = cr_flux.Nucleus(1, 1)
        self.pbar = cr_flux.Nucleus(1, -1)
        self.a    = cr_flux.Nucleus(4, 2)
        self.n    = cr_flux.Nucleus(1, 0)
    
    def testInit(self):
        """
        Test for __init__
        """
        self.assertRaises(TypeError, cr_flux.Nucleus, 0, 0)
        self.assertRaises(TypeError, cr_flux.Nucleus, 4, 5)
        self.assertRaises(TypeError, cr_flux.Nucleus, 1, -2)
        self.assertRaises(TypeError, cr_flux.Nucleus, 2, 0)
        
        self.assertEqual(self.p.A, 1)
        self.assertEqual(self.p.Z, 1)
        self.assertEqual(self.p.mass, 938.27203)
        self.assertEqual(self.pbar.Z, -1)
        self.assertEqual(self.pbar.A, 1)
        self.assertEqual(self.pbar.mass, 938.27203)
        self.assertEqual(self.n.A, 1)
        self.assertEqual(self.n.Z, 0)
        self.assertEqual(self.n.mass, 939.56536)
        self.assertEqual(self.a.A, 4)
        self.assertEqual(self.a.Z, 2)
        self.assertEqual(self.a.mass, (938.27203 + 939.56536)*2)

    def testComp(self):
        """
        Test for __cmp__
        """
        self.assertEqual(self.p, self.p)
        self.assertEqual(self.a, self.a)
        self.assertNotEqual(self.p, self.a)
        
        self.assertEqual(self.p, cr_flux.Nucleus(1, 1))
        self.assertEqual(self.a, cr_flux.Nucleus(4, 2))
        self.assertEqual(self.n, cr_flux.Nucleus(1, 0))
        self.assertEqual(self.pbar, cr_flux.Nucleus(1, -1))
        
    def testHash(self):
        """
        Test for __hash__
        """
        d = {}
        d[cr_flux.Nucleus(1, 1)] = 1
        d[cr_flux.Nucleus(4, 2)] = 4
        self.assertEqual(d[cr_flux.Nucleus(1, 1)], 1)
        self.assertEqual(d[cr_flux.Nucleus(4, 2)], 4)

class TestPhoton(unittest.TestCase):        
    """
    Unit test for class Photon
    """
    def setUp(self):
        self.photon = cr_flux.Photon()

    def testComp(self):
        """
        Test for __cmp__
        """
        self.assertEqual(self.photon, cr_flux.Photon())
        
    def testHash(self):
        """
        Test for __hash__
        """
        d = {}
        d[self.photon] = 1
        self.assertEqual(d[cr_flux.Photon()], 1)
        
class TestCRAbundance(unittest.TestCase):
    """
    Unit test for class CRAbundance
    """
    def setUp(self):
        self.ab = cr_flux.CRAbundance()
        self.ab.add_cr(cr_flux.Nucleus(4, 2), 0.1)
        self.ab.add_cr(cr_flux.Nucleus(7, 3), 0.02)
    
    def testCRAbundance(self):
        """
        Test for cr_abundance
        """
        self.assertEqual(self.ab.cr_abundance(cr_flux.proton), 1)
        self.assertEqual(self.ab.cr_abundance(cr_flux.alpha), 0.1)
        self.assertEqual(self.ab.cr_abundance(cr_flux.Nucleus(6, 1)), 0)
        
    def testExpandFluxList(self):
        """
        Test for expand_flux_list
        """
        arc = cr_flux.FluxModelArchive()
        Fp_galp = arc.galprop(cr_flux.proton, 8.5, 0.)
        Fa_galp = arc.galprop(cr_flux.alpha, 8.5, 0.)
        cr_list = [Fp_galp, Fa_galp]
        self.ab.expand_flux_list(cr_list)
        self.assertEqual(True, numpy.all(Fp_galp.F == cr_list[0].F))
        self.assertEqual(True, numpy.all(Fa_galp.F == cr_list[1].F))
        self.assertEqual(True, numpy.all(Fa_galp.F*(0.02/0.1) == cr_list[2].F))
        
class TestISMAbundance(unittest.TestCase):
    """
    Unit test for class ISMAbundance
    """
    def setUp(self):
        self.ab = cr_flux.ISMAbundance()
        self.ab.add_ism(cr_flux.Nucleus(4, 2), 0.1)
    
    def testISMAbundance(self):
        """
        Test for ism_abundance
        """
        self.assertEqual(self.ab.ism_abundance(cr_flux.proton), 1)
        self.assertEqual(self.ab.ism_abundance(cr_flux.alpha), 0.1)
        self.assertEqual(self.ab.ism_abundance(cr_flux.Nucleus(6, 1)), 0)
        
class TestMultiplicity(unittest.TestCase):
    """
    Unit test for class Multiplicity
    """
    def setUp(self):
        self.mul = cr_flux.Multiplicity()
        self.mul.add_mul(cr_flux.proton, cr_flux.alpha, 3)
        
    def testMultiplicity(self):
        self.assertEqual(self.mul.multiplicity(cr_flux.proton, cr_flux.proton), 1)
        self.assertEqual(self.mul.multiplicity(cr_flux.proton, cr_flux.alpha), 3)
        self.assertEqual(self.mul.multiplicity(cr_flux.alpha, cr_flux.proton), 3)
        self.assertEqual(self.mul.multiplicity(cr_flux.alpha, cr_flux.alpha), 16)
        self.mul.force = False
        self.assertEqual(self.mul.multiplicity(cr_flux.alpha, cr_flux.alpha), 0.)

class TestGaisserTable(unittest.TestCase):
    """
    Unit test for class GaisserTable
    """
    def setUp(self):
        self.gaisser = cr_flux.gaisser_table
        self.H  = cr_flux.Nucleus(1, 1)
        self.He = cr_flux.Nucleus(4, 2)
        self.Fe = cr_flux.Nucleus(56, 26)
    
    def testCRAbundance(self):
        self.assertEqual(self.gaisser.cr_abundance(self.H), 1)
        self.assertEqual(self.gaisser.cr_abundance(self.He), 0.042)
        self.assertEqual(self.gaisser.cr_abundance(self.Fe), 1.1e-4)

    def testMultiplicity(self):
        self.assertEqual(self.gaisser.multiplicity(self.H, self.H), 1)
        self.assertEqual(self.gaisser.multiplicity(self.H, self.He), 3.57)
        self.assertEqual(self.gaisser.multiplicity(self.H, self.Fe), 40)
        self.assertEqual(self.gaisser.multiplicity(self.He, self.H), 3.57)
        self.assertEqual(self.gaisser.multiplicity(self.He, self.He), 12.6)
        self.assertEqual(self.gaisser.multiplicity(self.He, self.Fe), 135)
        self.assertEqual(self.gaisser.multiplicity(self.Fe, self.H), 40)
        self.assertEqual(self.gaisser.multiplicity(self.Fe, self.He), 135)
        self.assertEqual(self.gaisser.multiplicity(self.Fe, self.Fe), 56*56)

class TestGas(unittest.TestCase):
    """
    Unit test for class Gas
    """
    def setUp(self):
        self.gas = cr_flux.Gas(cr_flux.CRAbundance(), cr_flux.ISMAbundance())        

class TestGasDensity(unittest.TestCase):
    """
    Unit test for class GasDensity
    """
    def setUp(self):
        self.gas = cr_flux.GasDensity(cr_flux.CRAbundance(), cr_flux.ISMAbundance(), 1)        
        
    def testInit(self):
        self.assertEqual(self.gas.density, 1)

class TestCloud(unittest.TestCase):
    """
    Unit test for class Cloud
    """
    def setUp(self):
        self.cloud = cr_flux.Cloud(cr_flux.CRAbundance(), cr_flux.ISMAbundance(), 1e5, 500)        
        
    def testInit(self):
        self.assertEqual(self.cloud.mass, 1e5)
        self.assertEqual(self.cloud.dist, 500)
        
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestParticle))
    suite.addTest(unittest.makeSuite(TestChargedLepton))
    suite.addTest(unittest.makeSuite(TestNucleus))
    suite.addTest(unittest.makeSuite(TestPhoton))
    suite.addTest(unittest.makeSuite(TestCRAbundance))
    suite.addTest(unittest.makeSuite(TestISMAbundance))
    suite.addTest(unittest.makeSuite(TestGaisserTable))
    suite.addTest(unittest.makeSuite(TestGas))
    suite.addTest(unittest.makeSuite(TestGasDensity))
    suite.addTest(unittest.makeSuite(TestCloud))

    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())            
