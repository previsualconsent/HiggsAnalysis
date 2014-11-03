import ROOT

"""
Style options mostly from CMS's tdrStyle.C
"""
def customROOTstyle() :
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(False)
    ROOT.gStyle.SetPadTopMargin(0.06);
    ROOT.gStyle.SetPadBottomMargin(0.13);
    ROOT.gStyle.SetPadLeftMargin(0.12);
    ROOT.gStyle.SetPadRightMargin(0.02);
    ROOT.gStyle.SetLabelColor(1, "XYZ");
    ROOT.gStyle.SetLabelFont(42, "XYZ");
    ROOT.gStyle.SetLabelOffset(0.007, "XYZ");
    ROOT.gStyle.SetLabelSize(0.05, "XYZ");
    ROOT.gStyle.SetTitleSize(0.05, "XYZ");
    ROOT.gStyle.SetAxisColor(1, "XYZ");
    ROOT.gStyle.SetStripDecimals(True);
    ROOT.gStyle.SetTickLength(0.03, "XYZ");
    ROOT.gStyle.SetNdivisions(510, "XYZ");
    ROOT.gStyle.SetPadTickX(0);
    ROOT.gStyle.SetPadTickY(0);
    ROOT.gStyle.SetMarkerStyle(20);
    ROOT.gStyle.SetHistLineColor(1);
    ROOT.gStyle.SetHistLineStyle(0);
    ROOT.gStyle.SetHistLineWidth(1);
    ROOT.gStyle.SetFrameBorderMode(0);
    ROOT.gStyle.SetFrameBorderSize(1);
    ROOT.gStyle.SetFrameFillColor(0);
    ROOT.gStyle.SetFrameFillStyle(0);
    ROOT.gStyle.SetFrameLineColor(1);
    ROOT.gStyle.SetFrameLineStyle(1);
    ROOT.gStyle.SetFrameLineWidth(1);

def saveHists(file,prefix=""):
    customROOTstyle()
    ROOT.gROOT.SetBatch(True)
    histObjectNames = ["TH1F", "TH2F", "TH1", "TH2"]
    for key in file.GetListOfKeys():
        if key.IsFolder():
            dir = file.Get(key.GetName())
            saveHists(dir,prefix=key.GetName())
        if key.GetClassName() in histObjectNames:
            hist = file.Get(key.GetName())
            drawHist(hist,prefix + '-' + key.GetName() + ".png")

def drawHist(hist,name,width=500,height=500):
    customROOTstyle()
    c = ROOT.TCanvas("c","c",width,height)
    #hist.SetLineWidth(2)
    hist.Draw()
    c.SaveAs(name)

def drawMultiple(hists,labels,filename,colors=[], width = 500, height = 500):
    customROOTstyle()
    hist_max = 0
    if not colors:
        colors = [ROOT.kRed, ROOT.kBlue, ROOT.kBlack]
        colors = colors[:len(hists)]
    for h in hists:
        if h.GetMaximum() > hist_max:
            hist_max = h.GetMaximum()

    canv = ROOT.TCanvas("c","c",width,height)
    first = True

    x1 = ROOT.gStyle.GetPadLeftMargin();
    x2 = 1 - ROOT.gStyle.GetPadRightMargin();
    y2 = 1 - ROOT.gStyle.GetPadTopMargin();
    y1 = y2*.9

    leg = ROOT.TLegend(x1,y1,x2,y2)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetNColumns(len(hists))
    for h,l,c in zip(hists,labels,colors):
        h.SetMaximum(1.2 * hist_max)
        h.SetTitle(l)
        h.SetLineColor(c)
        #h.SetLineWidth(2)
        #h.SetOptStat(0)
        if first:
            h.Draw()
            first = False
        else:
            h.Draw("same")

        leg.AddEntry(h,l)

    leg.Draw()
    canv.SaveAs(filename)
