"""
Example of plotting local intestellar spectrum (LIS)
$Id: lis.py,v 1.1 2009/05/17 11:32:46 oxon Exp $
"""

import cr_flux
import ROOT

arc = cr_flux.FluxModelArchive()

phi = [491., 591., 658., 1300., 1109]

Fp_galp = [arc.galprop(cr_flux.proton, 8.5, 0.),]
Fa_galp = [arc.galprop(cr_flux.alpha, 8.5, 0.),]

Fp_bess = [arc.BESS_LIS(cr_flux.proton, Fp_galp[0].E, Fp_galp[0].dEl, Fp_galp[0].dEh),]
Fa_bess = [arc.BESS_LIS(cr_flux.alpha,  Fa_galp[0].E, Fa_galp[0].dEl, Fa_galp[0].dEh),]

for p in phi:
    Fp_galp.append(Fp_galp[0].solar_modulate(p))
    Fa_galp.append(Fa_galp[0].solar_modulate(p))
    Fp_bess.append(Fp_bess[0].solar_modulate(p))
    Fa_bess.append(Fa_bess[0].solar_modulate(p))

can = ROOT.TCanvas("can", "can")
can.Divide(2, 1, 1e-10, 1e-10)

can.cd(1)
ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()

for i, f in enumerate(Fp_galp):
    f.set_index(0)
    if i == 0:
        f.graph.GetXaxis().SetLimits(1e2, 1e5) # 100 [MeV] to 100 [GeV]
        f.graph.GetHistogram().SetMaximum(1e-4)
        f.graph.GetHistogram().SetMinimum(1e-8)
        f.graph.Draw("al")
    else:
        f.graph.Draw("l same")

for i, f in enumerate(Fp_bess):
    f.set_index(0)
    f.graph.SetLineColor(2)
    f.graph.Draw("l same")

can.cd(2)
ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()

for i, f in enumerate(Fa_galp):
    f.set_index(0)
    if i == 0:
        f.graph.GetXaxis().SetLimits(1e2, 1e5) # 100 [MeV] to 100 [GeV]
        f.graph.GetHistogram().SetMaximum(1e-4)
        f.graph.GetHistogram().SetMinimum(1e-8)
        f.graph.Draw("al")
    else:
        f.graph.Draw("l same")

for i, f in enumerate(Fa_bess):
    f.set_index(0)
    f.graph.SetLineColor(2)
    f.graph.Draw("l same")
