import ROOT
from UserCode.HiggsAnalysis.PlotUtils import saveHists, drawMultiple


fPU = ROOT.TFile.Open("Higgs-analysis-PU.root")
fnoPU = ROOT.TFile.Open("Higgs-analysis-noPU.root")

colors = [ROOT.kBlue, ROOT.kBlack, ROOT.kRed, ROOT.kGreen]

h1 = fPU.Get("analysis/en_diff_endcap_r9")
h2 = fPU.Get("analysis/en_diff_barrel_r9")
drawMultiple([h1,h2],["Endcap","Barrel"],"en_diff.png",colors = colors, norm = True, xtitle = "(E_{gen} - E_{reco})/E_{gen}", ytitle = "Percent Events")


h1 = fPU.Get("analysis/diphoton_massE+B_r9")
h2 = fnoPU.Get("analysis/diphoton_massE+B_r9")

print type(h1),type(h2)
drawMultiple([h2,h1],["No Pileup","Pileup"],"diphoton_E+B.png",colors = colors, norm = True, ytitle = "% Events/2 GeV", rebin = 2)



h1 = fPU.Get("analysis/diphoton_massB+B_r9")
h2 = fnoPU.Get("analysis/diphoton_massB+B_r9")
print type(h1),type(h2)
drawMultiple([h2,h1],["No Pileup","Pileup"],"diphoton_B+B.png",colors = colors, norm = True, ytitle = "% Events/GeV")


