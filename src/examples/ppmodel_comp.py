"""
Example of comparison between Kamae et al. 2006 and Dermer 1986
$Id: ppmodel_comp.py,v 1.1 2009/02/23 16:07:47 oxon Exp $
"""

import cr_flux
import ROOT

# Definition of nuclei based on Gaisser and Schaefer (1992)
H = cr_flux.Nucleus(1, 1)
He= cr_flux.Nucleus(4, 2)

arc = cr_flux.FluxModelArchive()

# proton to Fe in CR, H and He in ISM
Fcr = [arc.galprop(H, 8.5, 0.), arc.galprop(He, 8.5, 0.)]
ism_ab = cr_flux.ISMAbundance()
cr_ab = cr_flux.CRAbundance()
multi = cr_flux.Multiplicity()
multi.add_mul(H, He, 3.57)
multi.add_mul(He, He, 12.6)
kamae = cr_flux.Kamae2006(multi)
dermer= cr_flux.Dermer1986(multi)
gas = cr_flux.GasDensity(cr_ab, ism_ab, 1) # 1 [p/cm^3]

Fg_base = cr_flux.AbsoluteSpectrum(cr_flux.Photon(), Fcr[0].E/100., Fcr[0].dEl/100., Fcr[0].dEh/100.)
Fg1 = kamae.gamma(gas, Fcr, Fg_base)
Fg2 = dermer.gamma(gas, Fcr, Fg_base)

# p-p interaction
cr_ab.add_cr(H, 1.)
cr_ab.add_cr(He, 0.)
ism_ab.add_ism(He, 0.)
ism_ab.add_ism(H, 1.)
Fk_pp = kamae.gamma(gas, Fcr[0], Fg_base) # p-p
Fd_pp = dermer.gamma(gas, Fcr[0], Fg_base) # p-p
Fm_pp = Fk_pp

cr_ab.add_cr(H, 1.)
cr_ab.add_cr(He, 0.)
ism_ab.add_ism(He, 0.11)
ism_ab.add_ism(H, 0.)
Fk_pa = kamae.gamma(gas, Fcr[0], Fg_base) # p-He
Fd_pa = dermer.gamma(gas, Fcr[0], Fg_base) # p-He
Fm_pa = Fd_pa

cr_ab.add_cr(H, 0.)
cr_ab.add_cr(He, 1.)
ism_ab.add_ism(He, 0.)
ism_ab.add_ism(H, 1.)
Fk_ap = kamae.gamma(gas, Fcr[1], Fg_base) # alpha-H
Fd_ap = dermer.gamma(gas, Fcr[1], Fg_base) # alpha-H
Fm_ap = Fd_ap

cr_ab.add_cr(H, 0.)
cr_ab.add_cr(He, 1.)
ism_ab.add_ism(He, 0.11)
ism_ab.add_ism(H, 0.)
Fk_aa = kamae.gamma(gas, Fcr[1], Fg_base) # alpha-He
Fd_aa = dermer.gamma(gas, Fcr[1], Fg_base) # alpha-He
Fm_aa = Fd_aa

Fk = Fk_pp + Fk_pa + Fk_ap + Fk_aa
Fd = Fd_pp + Fd_pa + Fd_ap + Fd_aa
Fm = Fm_pp + Fm_pa + Fm_ap + Fm_aa

can1 = ROOT.TCanvas("can1", "can1", 600, 600)

Fk.graph.SetTitle("Gamma-ray emissivity per unit density (n(H) = 1 cm^{-3})")
Fk.graph.GetXaxis().SetLimits(50, 1e4) # 50 [MeV] to 10 [GeV]
Fk.graph.GetHistogram().SetMinimum(1e-26)
Fk.graph.GetHistogram().SetMaximum(3e-23)
Fk.graph.Draw("al")

Fd.graph.SetLineStyle(2)
Fd.graph.Draw("l same")
Fm.graph.SetLineStyle(3)
Fm.graph.Draw("l same")

for f in Fk_pp, Fk_pa, Fk_ap, Fk_aa:
    f.graph.SetLineColor(2)
    f.graph.Draw("l same")

for f in Fd_pp, Fd_pa, Fd_ap, Fd_aa:
    f.graph.SetLineColor(4)
    f.graph.Draw("l same")

ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()
ROOT.gPad.Update()

leg = ROOT.TLegend(0.22, 0.15, 0.85, 0.32)
leg.SetFillColor(0)
en1 = leg.AddEntry(Fk.graph, "CR = (p, #alpha), ISM = (H, He), Kamae et al. 2006", "l")
en2 = leg.AddEntry(Fd.graph, "CR = (p, #alpha), ISM = (H, He), Dermer 1986", "l")
en3 = leg.AddEntry(Fm.graph, "CR = (p, #alpha), ISM = (H, He), Mix", "l")
en4 = leg.AddEntry(Fk_pp.graph, "Kamae et al. 2006", "l")
en5 = leg.AddEntry(Fd_pp.graph, "Dermer 1986", "l")
leg.Draw()

tex1 = ROOT.TLatex(650, 9e-24, "p-H")
tex2 = ROOT.TLatex(650, 3e-24, "p-He")
tex3 = ROOT.TLatex(650, 1.3e-24, "#alpha-H")
tex4 = ROOT.TLatex(650, 5e-25, "#alpha-He")

for tex in tex1, tex2, tex3, tex4:
    tex.SetTextSize(0.035)
    tex.Draw()

can2 = ROOT.TCanvas("can2", "can2")
ratio1 = ROOT.TGraph(Fk.E.size, Fk.E, Fk.F/Fk_pp.F)
ratio2 = ROOT.TGraph(Fd.E.size, Fd.E, Fd.F/Fk_pp.F)
ratio3 = ROOT.TGraph(Fm.E.size, Fm.E, Fm.F/Fk_pp.F)
ratio2.SetTitle(";Energy [MeV];Ratio to Kamae 2006 p-p")
ratio2.GetXaxis().CenterTitle()
ratio2.GetYaxis().CenterTitle()
ratio2.GetXaxis().SetLimits(10, 1e5) # 10 [MeV] to 100 [GeV]
ratio2.GetHistogram().SetMinimum(0.)
ratio2.GetHistogram().SetMaximum(2.)
ratio2.SetLineStyle(2)
ratio2.Draw("al")
ratio3.SetLineStyle(3)
ratio3.Draw("l same")
ROOT.gPad.SetLogx()

leg2 = ROOT.TLegend(0.22, 0.15, 0.85, 0.35)
leg2.SetFillColor(0)
leg2.AddEntry(ratio1, en1.GetLabel(), "l")
leg2.AddEntry(ratio2, en2.GetLabel(), "l")
leg2.AddEntry(ratio3, en3.GetLabel(), "l")

ratio_kpa = ROOT.TGraph(Fk_pa.E.size, Fk_pa.E, Fk_pa.F/Fk_pp.F)
ratio_kap = ROOT.TGraph(Fk_ap.E.size, Fk_ap.E, Fk_ap.F/Fk_pp.F)
ratio_kaa = ROOT.TGraph(Fk_aa.E.size, Fk_aa.E, Fk_aa.F/Fk_pp.F)
ratio_dpa = ROOT.TGraph(Fd_pa.E.size, Fd_pa.E, Fd_pa.F/Fd_pp.F)
ratio_dap = ROOT.TGraph(Fd_ap.E.size, Fd_ap.E, Fd_ap.F/Fd_pp.F)
ratio_daa = ROOT.TGraph(Fd_aa.E.size, Fd_aa.E, Fd_aa.F/Fd_pp.F)

for ratio in ratio_kpa, ratio_kap, ratio_kaa:
    ratio.SetLineColor(2)
    ratio.Draw("l same")

for ratio in ratio_dpa, ratio_dap, ratio_daa:
    ratio.SetLineColor(4)
    ratio.Draw("l same")

leg2.AddEntry(ratio_kpa, en4.GetLabel(), "l")
leg2.AddEntry(ratio_dpa, en5.GetLabel(), "l")

leg2.Draw()