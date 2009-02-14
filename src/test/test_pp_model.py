"""
Unit test script for pp_model.py
$Id: test_pp_model.py,v 1.2 2009/02/14 07:12:58 oxon Exp $
"""

import unittest

import cr_flux

class TestKamae2006(unittest.TestCase):
    """
    Unit test for class Kamae2006
    """
    def setUp(self):
        self.cr_ab = cr_flux.CRAbundance()
        self.ism_ab= cr_flux.ISMAbundance()
        self.multi = cr_flux.Multiplicity()
        self.kamae = cr_flux.Kamae2006(self.multi)
        arc = cr_flux.FluxModelArchive()
        self.Fp = arc.galprop(cr_flux.proton, 8.5, 0.)
        
    def testGamma(self):
        # CR = (p), ISM = (H)
        gas = cr_flux.GasDensity(self.cr_ab, self.ism_ab, 1)
        Fg1 = cr_flux.AbsoluteSpectrum(cr_flux.Photon(), self.Fp.E/100., self.Fp.dEl/100., self.Fp.dEh/100.)
        Fg1 = self.kamae.gamma(gas, self.Fp, Fg1)
        
        # CR = (p, alpha), ISM = (H)
        self.cr_ab.add_cr(cr_flux.alpha, 0.1)
        Fg2 = self.kamae.gamma(gas, self.Fp, Fg1)
        
        # CR = (p, alpha), ISM = (H, He)
        self.ism_ab.add_ism(cr_flux.alpha, 0.01)
        Fg3 = self.kamae.gamma(gas, self.Fp, Fg1)

        # 0.05 times smaller flux for F_alpha
        self.Fa = cr_flux.DiffuseSpectrum(cr_flux.alpha, self.Fp.E, self.Fp.dEl, self.Fp.dEh, self.Fp.F*0.05)
        Fg4 = self.kamae.gamma(gas, [self.Fp, self.Fa], Fg1)
        
        # CR = (p, alpha), ISM = (H, He), sigma_pa = 3
        self.multi.add_mul(cr_flux.proton, cr_flux.alpha, 3)
        Fg5 = self.kamae.gamma(gas, self.Fp, Fg1)
        
        for i in range(Fg1.E.size):
            self.assertAlmostEqual(Fg2.F[i]/Fg1.F[i], (1. + 0.1*4)/1., 5)
            self.assertAlmostEqual(Fg3.F[i]/Fg1.F[i], (1. + 0.1*4 + 0.01*4 + 0.01*0.1*4*4)/1., 5)
            self.assertAlmostEqual(Fg4.F[i]/Fg1.F[i], (1. + 0.05*4 + 0.01*4 + 0.01*0.05*4*4)/1., 5)
            self.assertAlmostEqual(Fg5.F[i]/Fg1.F[i], (1. + 0.1*3 + 0.01*3 + 0.01*0.1*4*4)/1., 5)
        
        Fg3 = self.kamae.gamma(gas, self.Fp, Fg1)

def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestKamae2006))

    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())            
