"""
Example of plotting gamma-ray flux from GMC
$Id: gmc_plot.py,v 1.2 2009/02/14 07:12:56 oxon Exp $
"""

import cr_flux
import ROOT

# Ignore nuclei heavier than proton
multi = cr_flux.Multiplicity()
kamae = cr_flux.Kamae2006(multi)

arc = cr_flux.FluxModelArchive()
# Assume GALPROP flux at the GMC position (R = 8.5 [kpc], z = 0 [kpc]) 
Fp_galp = arc.galprop(cr_flux.proton, 8.5, 0.)

# We calculate gamma-ray spectrum at the Earth
# Assume GMC is concentrated at a point
# Energy binning start from 100 times smaller value than GALPROP
Fg = cr_flux.PointSourceSpectrum(cr_flux.Photon(), Fp_galp.E/100., Fp_galp.dEl/100., Fp_galp.dEh/100.)

# Assume hevier particles do not exist
# Use default cosmic-ray and ISM abundances that consist of only proton and hydrogen
cr_ab = cr_flux.CRAbundance()
ism_ab = cr_flux.ISMAbundance()

# Orion A, mass = 100 k solar mass at 400 [pc] from the Earth
orion = cr_flux.Cloud(cr_ab, ism_ab, 1e5, 400)

Fg = kamae.gamma(orion, Fp_galp, Fg)

can = ROOT.TCanvas("can", "")

Fg.graph.GetXaxis().SetLimits(10, 1e5) # 10 [MeV] to 100 [GeV]
Fg.graph.GetHistogram().SetMaximum(1e-4)
Fg.graph.GetHistogram().SetMinimum(1e-8)
Fg.graph.Draw("al")

ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()
