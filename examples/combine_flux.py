"""
Example of combining CR flux data
$Id: combine_flux.py,v 1.1 2009/02/25 04:58:27 oxon Exp $
"""

import cr_flux
import ROOT

proton = H  = cr_flux.Nucleus(1, 1)
alpha  = He = cr_flux.Nucleus(4, 2)

gaisser = cr_flux.gaisser_table

arc = cr_flux.FluxModelArchive()

F_galp = [arc.galprop(proton, 8.5, 0.), arc.galprop(alpha, 8.5, 0.)]
F_atic = [arc.ATIC2(proton), arc.ATIC2(alpha)]

ism_ab= cr_flux.ISMAbundance()
ism_ab.add_ism(He, 0.11)
cr_ab = cr_flux.CRAbundance()
kamae = cr_flux.Kamae2006(gaisser)
gas = cr_flux.GasDensity(cr_ab, ism_ab, 1.) # 1 [p/cm^3]
Fg = cr_flux.AbsoluteSpectrum(cr_flux.Photon(), F_galp[0].E/10, F_galp[0].dEl/10, F_galp[0].dEh/10)

# Calculate gamma-ray emissivity using GALPROP CR flux
Fg1 = kamae.gamma(gas, [F_galp[0], F_galp[1]], Fg)

# Remove data points that are in the energy range of ATIC
for i in range(2):
    size = F_galp[i].E.size
    Emin = F_atic[i].E[0] - F_atic[i].dEl[0]
    Emax = F_atic[i].E[-1] + F_atic[i].dEh[-1]
    for j in range(size):
        E1 = F_galp[i].E[j] - F_galp[i].dEl[j]
        E2 = F_galp[i].E[j] + F_galp[i].dEh[j]
        if Emin <= E1 and E2 <= Emax:
            F_galp[i].F[j] = 0

for f in F_galp[0], F_galp[1], F_atic[0], F_atic[1]:
    f.set_index(2.5)

can = ROOT.TCanvas("can", "can")
F_galp[0].graph.GetXaxis().SetLimits(1e2, 1e8)
F_galp[0].graph.Draw("ap")
F_galp[1].graph.Draw("p same")
F_atic[0].graph.Draw("p same")
F_atic[1].graph.Draw("p same")
ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()

# Calculate gamma-ray emissivity using GALPROP CR flux and ATIC data
Fg2 = kamae.gamma(gas, [F_galp[0], F_galp[1], F_atic[0], F_atic[1]], Fg)

can2 = ROOT.TCanvas("can2", "can2")

Fg1.graph.GetXaxis().SetLimits(1e2, 1e8)
Fg1.graph.GetHistogram().SetMinimum(1e-26)
Fg1.graph.GetHistogram().SetMaximum(3e-23)
Fg1.graph.Draw("al")
Fg2.graph.SetLineStyle(2)
Fg2.graph.Draw("l same")
ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()

leg = ROOT.TLegend(0.15, 0.15, 0.4, 0.3)
leg.SetFillColor(0)
leg.AddEntry(Fg1.graph, "GALPROP", "l")
leg.AddEntry(Fg2.graph, "GALPROP + ATIC", "l")
leg.Draw()
