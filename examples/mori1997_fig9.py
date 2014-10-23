"""
Reproducibility test for Mori 1997 model
$Id: mori1997_fig9.py,v 1.1 2009/03/03 06:47:39 oxon Exp $
"""

import cr_flux
import ROOT

arc = cr_flux.FluxModelArchive()

# proton to Fe in CR, H and He in ISM
F_gal = arc.galprop(cr_flux.proton, 8.5, 0.) # dummy
Fcr = [arc.Mori1997(cr_flux.proton, F_gal.E, F_gal.dEl, F_gal.dEh),
       arc.Mori1997(cr_flux.alpha, F_gal.E, F_gal.dEl, F_gal.dEh)]
mori = cr_flux.Mori1997(cr_flux.gaisser_table)
ism_ab = cr_flux.ISMAbundance()
ism_ab.add_ism(cr_flux.proton, 1.)
ism_ab.add_ism(cr_flux.alpha, 0.07/0.93)

gas = cr_flux.GasDensity(cr_flux.gaisser_table, ism_ab, 1) # 1 [p/cm^3]

Fg_base = cr_flux.AbsoluteSpectrum(cr_flux.photon, Fcr[0].E/100., Fcr[0].dEl/100., Fcr[0].dEh/100.)
Fg = mori.gamma(gas, Fcr, Fg_base)
Fg.set_index(0.)

can1 = ROOT.TCanvas("can1", "can1", 600, 600)

Fg.graph.SetTitle("Gamma-ray emissivity per unit density (n(H) = 1 cm^{-3})")
Fg.graph.GetXaxis().SetLimits(1e1, 1e5) # 10 [MeV] to 100 [GeV]
Fg.graph.GetHistogram().SetMinimum(1e-33)
Fg.graph.GetHistogram().SetMaximum(1e-27)
Fg.graph.Draw("al")

ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()