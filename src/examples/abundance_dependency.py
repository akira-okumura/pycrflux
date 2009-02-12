"""
Example of plotting abundance dependency of gamma-ray emission
$Id: abundance_dependency.py,v 1.2 2009/02/12 18:21:32 oxon Exp $
"""

import cr_flux
import ROOT

# Definition of nuclei based on Gaisser and Schaefer (1992)
H = cr_flux.Nucleus(1, 1)
He= cr_flux.Nucleus(4, 2)

gaisser = cr_flux.gaisser_table

arc = cr_flux.FluxModelArchive()

# only proton
F1 = [arc.galprop(H, 8.5, 0.)]
cr_ab1 = cr_flux.CRAbundance()
ism_ab1= cr_flux.ISMAbundance()
multi1 = cr_flux.Multiplicity()
multi1.force = False # Do not calculate other than p-p
kamae1 = cr_flux.Kamae2006(multi1)
gas1 = cr_flux.GasDensity(cr_ab1, ism_ab1, 1) # 1 [p/cm^3]

# only proton and alpha (multiplicity = 4)
F2 = [arc.galprop(H, 8.5, 0.), arc.galprop(He, 8.5, 0.)]
cr_ab2 = cr_flux.CRAbundance()
ism_ab2= cr_flux.ISMAbundance()
ism_ab2.add_ism(He, 0.07/0.93) # Mori 1997
# assume multiplicity of p-alpha is 4
multi2 = cr_flux.Multiplicity()
kamae2 = cr_flux.Kamae2006(multi2)
gas2 = cr_flux.GasDensity(cr_ab2, ism_ab2, 1) # 1 [p/cm^3]

# only proton and alpha (multiplicity = 3.57)
F3 = [arc.galprop(H, 8.5, 0.), arc.galprop(He, 8.5, 0.)]
cr_ab3 = cr_flux.CRAbundance()
ism_ab3= cr_flux.ISMAbundance()
ism_ab3.add_ism(He, 0.07/0.93) # Mori 1997
# assume multiplicity of p-alpha is 3.57, alpha-alpha is 12.6
multi3 = cr_flux.Multiplicity()
multi3.add_mul(H, He, 3.57)
multi3.add_mul(He, He, 12.6)
kamae3 = cr_flux.Kamae2006(multi3)
gas3 = cr_flux.GasDensity(cr_ab3, ism_ab3, 1) # 1 [p/cm^3]

# proton to Fe in CR, H and He in ISM
F4 = [arc.galprop(H, 8.5, 0.), arc.galprop(He, 8.5, 0.)]
cr_ab4 = gaisser
ism_ab4= cr_flux.ISMAbundance()
ism_ab4.add_ism(He, 0.07/0.93) # Mori 1997
multi4 = gaisser
kamae4 = cr_flux.Kamae2006(multi4)
gas4 = cr_flux.GasDensity(cr_ab4, ism_ab4, 1) # 1 [p/cm^3]

Fg = cr_flux.AbsoluteSpectrum(cr_flux.Photon(), F1[0].E/100, F1[0].dEl/100, F1[0].dEh/100, F1[0].F*0, F1[0].dFl*0, F1[0].dFh*0)
Fg1 = kamae1.gamma(gas1, F1, Fg)
Fg2 = kamae2.gamma(gas2, F2, Fg)
Fg3 = kamae3.gamma(gas3, F3, Fg)
Fg4 = kamae4.gamma(gas4, F4, Fg)

can1 = ROOT.TCanvas("can1", "can1")

Fg1.graph.SetTitle("Gamma-ray emissivity per unit density (n(H) = 1 cm^{-3})")
Fg1.graph.GetXaxis().SetLimits(10, 1e5) # 10 [MeV] to 100 [GeV]
Fg1.graph.GetHistogram().SetMaximum(1e-22)
Fg1.graph.GetHistogram().SetMinimum(1e-26)
Fg1.graph.Draw("al")
Fg2.graph.SetLineStyle(2)
Fg2.graph.Draw("l same")
Fg3.graph.SetLineStyle(3)
Fg3.graph.Draw("l same")
Fg4.graph.SetLineStyle(4)
Fg4.graph.Draw("l same")

ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()
ROOT.gPad.Update()

leg = ROOT.TLegend(0.22, 0.15, 0.85, 0.35)
leg.SetFillColor(0)
en1 = leg.AddEntry(Fg1.graph, "CR = p, ISM = H", "l")
en2 = leg.AddEntry(Fg2.graph, "CR = (p, #alpha), ISM = (H, He), #sigma_{p#alpha}/#sigma_{pp}=4, #sigma_{#alpha#alpha}/#sigma_{pp}=16", "l")
en3 = leg.AddEntry(Fg3.graph, "CR = (p, #alpha), ISM = (H, He), #sigma_{p#alpha}/#sigma_{pp}=3.57, #sigma_{#alpha#alpha}/#sigma_{pp}=12.6", "l")
en4 = leg.AddEntry(Fg4.graph, "CR = (p, #alpha, ... Fe), ISM = (H, He), Gaisser's multiplicity", "l")
leg.Draw()

can2 = ROOT.TCanvas("can2", "can2")
ratio2 = ROOT.TGraph(Fg1.E.size, Fg1.E, Fg2.F/Fg1.F)
ratio3 = ROOT.TGraph(Fg1.E.size, Fg1.E, Fg3.F/Fg1.F)
ratio4 = ROOT.TGraph(Fg1.E.size, Fg1.E, Fg4.F/Fg1.F)
ratio2.SetTitle(";Energy [MeV];Ratio of Emissivity")
ratio2.GetXaxis().CenterTitle()
ratio2.GetYaxis().CenterTitle()
ratio2.GetXaxis().SetLimits(10, 1e5) # 10 [MeV] to 100 [GeV]
ratio2.GetHistogram().SetMinimum(1.4)
ratio2.GetHistogram().SetMaximum(1.6)
ratio2.SetLineStyle(2)
ratio2.Draw("al")
ratio3.SetLineStyle(3)
ratio4.SetLineStyle(4)
ratio3.Draw("l same")
ratio4.Draw("l same")
ROOT.gPad.SetLogx()

leg2 = ROOT.TLegend(0.22, 0.15, 0.85, 0.35)
leg2.SetFillColor(0)
leg2.AddEntry(ratio2, en2.GetLabel(), "l")
leg2.AddEntry(ratio3, en3.GetLabel(), "l")
leg2.AddEntry(ratio4, en4.GetLabel(), "l")
leg2.Draw()
