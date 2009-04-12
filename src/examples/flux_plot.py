"""
Example of plotting various data of proton and alpha flux
$Id: flux_plot.py,v 1.3 2009/04/12 02:48:26 oxon Exp $
"""

import cr_flux
import ROOT

arc = cr_flux.FluxModelArchive()

# Create proton and alpha flux from various experimetal data
F_bess = [arc.BESS_TeV(cr_flux.proton), arc.BESS_TeV(cr_flux.alpha)]
F_ams  = [arc.AMS01(cr_flux.proton), arc.AMS01(cr_flux.alpha)]
F_atic = [arc.ATIC2(cr_flux.proton), arc.ATIC2(cr_flux.alpha)]
F_runj = [arc.RUNJOB(cr_flux.proton), arc.RUNJOB(cr_flux.alpha)]
# R = 8.5 [kpc], z = 0 [kpc]
F_galp = [arc.galprop(cr_flux.proton, 8.5, 0.), arc.galprop(cr_flux.alpha, 8.5, 0.)]
F_mori = [arc.Mori1997(cr_flux.proton, F_galp[0].E, F_galp[0].dEl, F_galp[0].dEh),
          arc.Mori1997(cr_flux.alpha, F_galp[0].E, F_galp[0].dEl, F_galp[0].dEh)]
F_honda= [arc.Honda2004(cr_flux.proton, F_galp[0].E, F_galp[0].dEl, F_galp[0].dEh),
          arc.Honda2004(cr_flux.alpha, F_galp[0].E, F_galp[0].dEl, F_galp[0].dEh)]

# Chanege the flux from dN/dE to E^2.5 x dN/dE
for f in F_bess, F_ams, F_atic, F_runj, F_galp, F_mori, F_honda:
    f[0].set_index(2.5)
    f[1].set_index(2.5)

leg = ROOT.TLegend(0.55, 0.67, 0.88, 0.88)
leg.SetFillColor(0)

for F, c, title in (F_bess, ROOT.kOrange + 7, "BESS-TeV"),\
                   (F_ams,  ROOT.kGreen + 3, "AMS01"),\
                   (F_atic, ROOT.kMagenta + 2, "ATIC-2"),\
                   (F_runj, ROOT.kBlue + 2, "RUNJOB"):
    F[0].graph.SetMarkerStyle(21)
    F[1].graph.SetMarkerStyle(21)
    F[0].graph.SetMarkerColor(c)
    F[1].graph.SetMarkerColor(c)

    leg.AddEntry(F[0].graph, title, "p")

for F, c, title in (F_galp, 2, "GALPROP (R = 8.5 [kpc])"),\
                   (F_mori, 4, "Mori 1997 'median'"),\
                   (F_honda, 6, "Honda et al. 2004"):   
    F[0].graph.SetLineColor(c)
    F[1].graph.SetLineColor(c)
    
    leg.AddEntry(F[0].graph, title, "l")

F_bess[0].graph.GetHistogram().SetMinimum(1e2)
F_bess[0].graph.GetHistogram().SetMaximum(1e5)
F_bess[0].graph.GetXaxis().SetLimits(0.8*1e3, 1e6*1e3)

can = ROOT.TCanvas("can", "", 600, 600)

F_bess[0].graph.Draw("ap")
F_bess[1].graph.Draw("p same")

for F in F_ams, F_atic, F_runj:
    F[0].graph.Draw("p same")
    F[1].graph.Draw("p same")
    
for F in F_galp, F_mori, F_honda:
    F[0].graph.Draw("l same")
    F[1].graph.Draw("l same")

ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()

leg.Draw()

tex = ROOT.TLatex(1e4, 2e3, "#alpha")
tex.Draw()
tex2 = ROOT.TLatex(1e4, 4e4, "proton")
tex2.Draw()
