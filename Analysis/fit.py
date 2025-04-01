import ROOT
from setTDRStyle import setTDRStyle


def fitAndPlot(ws, obj):
	
	fitResult = ws.pdf("model").fitTo(ws.data("hist%sRoo"%obj), ROOT.RooFit.Save())


	c1 = ROOT.TCanvas("c1","c1",700,700)
	c1.cd()	
	plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
	style = setTDRStyle()
	ROOT.gStyle.SetOptStat(0)
	plotPad.UseCurrentStyle()
	plotPad.Draw()	
	plotPad.cd()

	ws.var("mass").setBins(40)
	frame = ws.var('mass').frame(ROOT.RooFit.Title('Invariant mass of dimuon pairs'))
	frame.GetXaxis().SetTitle('m_{#mu#mu} [GeV]')
	frame.GetYaxis().SetTitle("Events / 0.01 GeV")
	ROOT.RooAbsData.plotOn(ws.data('hist%sRoo'%obj), frame,ROOT.RooFit.Name("hist%sRoo"%obj))
	ws.pdf('model').plotOn(frame,ROOT.RooFit.Name("model"))
	frame.Draw()


	latex = ROOT.TLatex()
	latex.SetTextFont(42)
	latex.SetTextAlign(31)
	latex.SetTextSize(0.04)
	latex.SetNDC(True)
	latexCMS = ROOT.TLatex()
	latexCMS.SetTextFont(61)
	latexCMS.SetTextSize(0.08)
	latexCMS.SetNDC(True)
	latexCMSExtra = ROOT.TLatex()
	latexCMSExtra.SetTextFont(52)
	latexCMSExtra.SetTextSize(0.045)
	latexCMSExtra.SetNDC(True) 
		
	latex.DrawLatex(0.95, 0.96, "(13.6 TeV)")

	cmsExtra = "#splitline{Simulation}{Preliminary}"
	latexCMS.DrawLatex(0.19,0.85,"CMS")
	if "Simulation" in cmsExtra:
		yLabelPos = 0.78	
	else:
		yLabelPos = 0.78	

	latexCMSExtra.DrawLatex(0.19,yLabelPos,"%s"%(cmsExtra))				

	latexCMSExtra.DrawLatex(0.19,yLabelPos,"%s"%(cmsExtra))				

	latexEta = ROOT.TLatex()
	latexEta.SetTextFont(42)
	latexEta.SetTextAlign(31)
	latexEta.SetTextSize(0.035)
	latexEta.SetNDC(True)		
	if obj == "Trk":
		latexEta.DrawLatex(0.55,0.7,"Track fired 2Mu1Tk trigger")
	else:	
		latexEta.DrawLatex(0.55,0.7,"Track fired 3Mu trigger")

	c1.SaveAs("DsPhiMuNuFit_Data_%s.pdf"%obj)
	

tree = ROOT.TChain()
tree.Add("AnalysedTree_data_phimunu_2022_phimunu_merged_test5.root/FinalTree")

histTrk = ROOT.TH1F('trkHist','trkHist',40, 0.8,1.2)
histMu = ROOT.TH1F('muHist','muHist',40, 0.8,1.2)

#tree.Draw("phi_mass>>trkHist", "mu3Matched2Mu1Tk > 0", "goff", -1)
#tree.Draw("phi_mass>>muHist", "mu3Matched3Mu > 0", "goff", -1)

for ev in tree:
	if ev.mu3Matched2Mu1Tk > 0:
		histTrk.Fill(ev.phi_mass)
	if ev.mu3Matched3Mu > 0:
		histMu.Fill(ev.phi_mass)

print (histTrk.GetEntries())
print (histMu.GetEntries())

ws = ROOT.RooWorkspace("w")
mass = ROOT.RooRealVar('mass','mass',1.0, 0.8, 1.2 )
getattr(ws,'import')(mass,ROOT.RooCmdArg())			
histTrkRoo = ROOT.RooDataHist("histTrkRoo","histTrkRoo",ROOT.RooArgList(ws.var("mass")),histTrk)
getattr(ws,'import')(histTrkRoo, ROOT.RooCmdArg())
histMuRoo = ROOT.RooDataHist("histMuRoo","histMuRoo",ROOT.RooArgList(ws.var("mass")),histMu)
getattr(ws,'import')(histMuRoo, ROOT.RooCmdArg())

ws.factory("RooCBShape::cb(mass, mean[1.0,-0.8,1.2], sigma[0.02,0,10], alphaL[3,-25,25], nL[5,-25,25])")
ws.factory("RooExponential::exp(mass, tau[-1.0, -10, 0])")

fsig = ROOT.RooRealVar("fsig","signal fraction",0.5,0.,1.)

# model(x) = fsig*sig(x) + (1-fsig)*bkg(x)
model = ROOT.RooAddPdf("model","model",ROOT.RooArgList(ws.pdf("cb"),ws.pdf("exp")),ROOT.RooArgList(fsig))
getattr(ws,'import')(model, ROOT.RooCmdArg())

fitAndPlot(ws, "Trk")
fitAndPlot(ws, "Mu")
