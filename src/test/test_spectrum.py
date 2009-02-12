"""
Unit test script for spectrum.py
$Id: test_spectrum.py,v 1.1 2009/02/12 18:12:35 oxon Exp $ 
"""

import numpy
import unittest

import cr_flux

photon = cr_flux.Photon()
proton = cr_flux.Nucleus(1, 1)
alpha  = cr_flux.Nucleus(4, 2)

class TestSpectrum(unittest.TestCase):
    """
    Unit test for class Spectrum
    """
    def setUp(self):
        self.n = 5
        self.class_list = (cr_flux.Spectrum, cr_flux.AbsoluteSpectrum,
                           cr_flux.DiffuseSpectrum, cr_flux.PointSourceSpectrum)    
        self.E = 10**numpy.arange(self.n)
        self.e = 0.5*10**numpy.arange(self.n)
        self.F = 2*10**numpy.arange(self.n)
        self.Fg = []
        self.Fp = []
        self.Fa = []
        for f, p in (self.Fg, photon), (self.Fp, proton), (self.Fa, alpha):
            for cls in self.class_list:
                f.append(cls(p, self.E, self.e, self.e, self.F, self.e, self.e))
             
    def testInit(self):
        """
        Test function for __init__
        """
        self.assertRaises(TypeError, cr_flux.Spectrum, 1, self.E, self.e, self.e, self.F, self.e, self.e)
        self.assertRaises(TypeError, cr_flux.Spectrum, proton, self.E[1:], self.e, self.e, self.F, self.e, self.e)
        self.assertRaises(TypeError, cr_flux.Spectrum, proton, 1, self.e, self.e, self.F, self.e, self.e)

    def testAdd(self):
        """
        Test function for __add__
        """
        self.assertRaises(TypeError, self.Fg[0].__add__,  self.Fg[0], self.Fg[1])
        self.assertRaises(TypeError, self.Fg[0].__add__, self.Fg[0], self.Fp[0])
        self.assertRaises(TypeError, self.Fp[0].__add__, self.Fp[0], self.Fa[0])

        tmpspec0 = cr_flux.Spectrum(photon, self.E, self.e, self.e, self.F, self.e, self.e)
        tmpspec1 = cr_flux.Spectrum(photon, self.E*0.1, self.e, self.e, self.F, self.e, self.e)
        tmpspec2 = cr_flux.Spectrum(photon, self.E, self.e*0.1, self.e, self.F, self.e, self.e)
        tmpspec3 = cr_flux.Spectrum(photon, self.E, self.e, self.e*0.1, self.F, self.e, self.e)

        self.assertRaises(TypeError, tmpspec0.__add__, tmpspec0, tmpspec1)
        self.assertRaises(TypeError, tmpspec0.__add__, tmpspec0, tmpspec2)
        self.assertRaises(TypeError, tmpspec0.__add__, tmpspec0, tmpspec3)

        res = self.Fg[0] + self.Fg[0]
        
        self.assertEqual(True, numpy.all(res.F==self.Fg[0].F*2))
        self.assertEqual(True, numpy.all(res.dFl==(self.Fg[0].dFl**2 + self.Fg[0].dFl**2)**0.5))
        self.assertEqual(True, numpy.all(res.dFh==(self.Fg[0].dFh**2 + self.Fg[0].dFh**2)**0.5))
        
    def testMul(self):
        """
        Test function for __mul__
        """
        spec = cr_flux.Spectrum(photon, self.E, self.e, self.e, self.F, self.e, self.e)
        spec2 = spec*2.
        
        self.assertEqual(spec.F[0]*2, spec2.F[0])
        self.assertEqual(spec.dFl[0]*2, spec2.dFl[0])
        self.assertEqual(spec.dFh[0]*2, spec2.dFh[0])
        
    def testTitle(self):
        """
        Test for init_graph and set_index
        """
        
        # Check X titles
        for i in range(len(self.class_list)):
            self.assertEqual(self.Fg[i].graph.GetXaxis().GetTitle(), "Energy [MeV]")
            self.assertEqual(self.Fp[i].graph.GetXaxis().GetTitle(),\
                             "Kinetic Energy [MeV/n]")
            self.assertEqual(self.Fa[i].graph.GetXaxis().GetTitle(),\
                             "Kinetic Energy [MeV/n]")
        
        # Check Y titles
        self.assertEqual(self.Fg[0].graph.GetYaxis().GetTitle(),\
                         "E^{2} dN/dE [MeV^{2}/s/MeV]")
        self.assertEqual(self.Fg[1].graph.GetYaxis().GetTitle(),\
                         "E^{2} dN/dE [MeV^{2}/s/MeV]")
        self.assertEqual(self.Fg[2].graph.GetYaxis().GetTitle(),\
                         "E^{2} dN/dE [MeV^{2}/cm^{2}/s/sr/MeV]")
        self.assertEqual(self.Fg[3].graph.GetYaxis().GetTitle(),\
                         "E^{2} dN/dE [MeV^{2}/cm^{2}/s/MeV]")
        
        for i in range(len(self.class_list)):
            self.Fg[i].set_index(2.50)

        # Check Y again
        self.assertEqual(self.Fg[0].graph.GetYaxis().GetTitle(),\
                        "E^{2.5} dN/dE [MeV^{2.5}/s/MeV]")
        self.assertEqual(self.Fg[1].graph.GetYaxis().GetTitle(),\
                         "E^{2.5} dN/dE [MeV^{2.5}/s/MeV]")
        self.assertEqual(self.Fg[2].graph.GetYaxis().GetTitle(),\
                         "E^{2.5} dN/dE [MeV^{2.5}/cm^{2}/s/sr/MeV]")
        self.assertEqual(self.Fg[3].graph.GetYaxis().GetTitle(),\
                         "E^{2.5} dN/dE [MeV^{2.5}/cm^{2}/s/MeV]")
            
        for f in self.Fp, self.Fa:
            self.assertEqual(f[0].graph.GetYaxis().GetTitle(),\
                             "E^{2} dN/dE [(MeV/n)^{2}/s/(MeV/n)]")
            self.assertEqual(f[1].graph.GetYaxis().GetTitle(),\
                             "E^{2} dN/dE [(MeV/n)^{2}/s/(MeV/n)]")
            self.assertEqual(f[2].graph.GetYaxis().GetTitle(),\
                             "E^{2} dN/dE [(MeV/n)^{2}/cm^{2}/s/sr/(MeV/n)]")
            self.assertEqual(f[3].graph.GetYaxis().GetTitle(),\
                             "E^{2} dN/dE [(MeV/n)^{2}/cm^{2}/s/(MeV/n)]")
            
            for i in range(len(self.class_list)):
                f[i].set_index(2.50)
                
            self.assertEqual(f[0].graph.GetYaxis().GetTitle(),\
                             "E^{2.5} dN/dE [(MeV/n)^{2.5}/s/(MeV/n)]")
            self.assertEqual(f[1].graph.GetYaxis().GetTitle(),\
                             "E^{2.5} dN/dE [(MeV/n)^{2.5}/s/(MeV/n)]")
            self.assertEqual(f[2].graph.GetYaxis().GetTitle(),\
                             "E^{2.5} dN/dE [(MeV/n)^{2.5}/cm^{2}/s/sr/(MeV/n)]")
            self.assertEqual(f[3].graph.GetYaxis().GetTitle(),\
                             "E^{2.5} dN/dE [(MeV/n)^{2.5}/cm^{2}/s/(MeV/n)]")
            
class TestFluxModelArchive(unittest.TestCase):
    """
    Unit test for FluxModelArchive
    """
    def setUp(self):
        self.arc = cr_flux.FluxModelArchive()
        
    def testGalprop(self):
        """
        Teset for galprop
        """
        self.assertRaises(TypeError, self.arc.galprop, cr_flux.Photon(), 0, 0)
        f00 = self.arc.galprop(cr_flux.Nucleus(1, 1), 0, 0)
        f01 = self.arc.galprop(cr_flux.Nucleus(1, 1), 0, 0.1)
        f10 = self.arc.galprop(cr_flux.Nucleus(1, 1), 1, 0)
        f11 = self.arc.galprop(cr_flux.Nucleus(1, 1), 1, 0.1)
        
        # f0_ = <f0?>, f_0 = <f?0>, f1_ = <f1?>, f_1 = <f?1>
        f0_ = self.arc.galprop(cr_flux.Nucleus(1, 1), 0, 0.05)
        f_0 = self.arc.galprop(cr_flux.Nucleus(1, 1), 0.5, 0)
        f1_ = self.arc.galprop(cr_flux.Nucleus(1, 1), 1, 0.05)
        f_1 = self.arc.galprop(cr_flux.Nucleus(1, 1), 0.5, 0.1)
        
        self.assertEqual(f00.E[0], 100.) # 100 [MeV]
        self.assertAlmostEqual(f00.F[0], 2.586473/100./100., 10)
        self.assertAlmostEqual(f01.F[0], 2.551192/100./100., 10)
        self.assertAlmostEqual(f10.F[0], 8.815813/100./100., 10)
        self.assertAlmostEqual(f11.F[0], 8.279984/100./100., 10)
        
        self.assertAlmostEqual((f00.F[0] + f01.F[0])/2., f0_.F[0], 10)
        self.assertAlmostEqual((f00.F[0] + f10.F[0])/2., f_0.F[0], 10)
        self.assertAlmostEqual((f01.F[0] + f11.F[0])/2., f_1.F[0], 10)
        self.assertAlmostEqual((f10.F[0] + f11.F[0])/2., f1_.F[0], 10)
        
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSpectrum))
    suite.addTest(unittest.makeSuite(TestFluxModelArchive))

    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())            
