import ROOT
import numpy as np
import sys
import array

fn_peds = {}

fin = open(sys.argv[1])
for line in fin:
  fn_peds[line.split()[0]] = float(line.split()[1].rstrip())

thresholds = []
widths = []
measured_thresh = []

draw = False

hs = []
for fn, thresh in sorted(fn_peds.iteritems(), key=lambda (k,v): (v,k)):
  print "%s: %s" % (fn, thresh)
  hs.append(ROOT.TH1F("h_%s" % fn,"h_%s" % fn,100,125,525))
  f = ROOT.TFile(fn)
  t = f.Get("makeSD/sddiag")
  for i in range(t.GetEntries()):
    t.GetEntry(i)
    if t.mcproc != 40:
      continue
    if t.straw < 88 or t.panel != 1:
      continue

    adc = [t.adc[j] for j in range(15)]
    peak = np.max(adc)
    ped = np.mean(adc[:3])
    hs[-1].Fill((peak-ped))
  hs[-1].SetDirectory(0)
  f.Close()

import scipy.stats
from scipy.optimize import curve_fit
if draw:
  import matplotlib.pyplot as plt

def fcn(x,mu,sig):
  return scipy.stats.norm(loc=mu,scale=sig).cdf(x)

for i in range(len(hs)):
  hs[i].SetLineColor(i+1)
  if i==0:
    pass
  else:
#    hs[0].Draw()
#    hs[i].Draw("same")
#    raw_input()
    hs[i].Divide(hs[0])
    x = np.zeros(hs[i].GetNbinsX())
    y = np.zeros(hs[i].GetNbinsX())
    for j in range(hs[i].GetNbinsX()):
      x[j] = hs[i].GetXaxis().GetBinCenter(j+1)
      y[j] = hs[i].GetBinContent(j+1)
  
    popt, pcov = curve_fit(fcn,x,y,p0=[200,10])
    yfit = scipy.stats.norm(loc=popt[0],scale=popt[1]).cdf(x)
    hfit = ROOT.TH1F("","",100,125,525)
    for j in range(hfit.GetNbinsX()):
      hfit.SetBinContent(j+1,yfit[j])
  
    print popt
    print pcov
    print popt[1]/popt[0]
  
    if i > 0:
      hs[i].Draw()
      hfit.SetLineColor(ROOT.kRed)
      hfit.Draw("same")
      raw_input()
    thresholds.append(popt[0])
    widths.append(popt[1])

    if draw:
      plt.plot(x,y)
      plt.plot(x,yfit)
      plt.show()

for fn, thresh in sorted(fn_peds.iteritems(), key=lambda (k,v): (v,k)):
  measured_thresh.append(thresh)

c1 = ROOT.TCanvas("c1","c1",800,600)
g = ROOT.TGraph(len(thresholds),np.array(measured_thresh[1:]),np.array(thresholds))
g.SetMarkerStyle(21)
#fr = mg.Fit("pol1","SQ")
g.SetTitle("Fe55 Threshold")
g.Draw("ALP")
#l.Draw()
g.GetXaxis().SetTitle("Set threshold (mV)")
g.GetYaxis().SetTitle("Output signal threshold (peak minus pedestal counts)")
#params = fr.GetParams()
#errors = fr.GetErrors()
#print "Signal threshold (counts pmp): (",params[0],"+-",errors[0],") + (",params[1],"+-",errors[1],") * thresh [mV]"
#p = ROOT.TPaveText(0.5,0.2,0.7,0.3,"NDC NB")
#p.AddText("signal cutoff (pmp) = (%.1f #pm %.1f) + (%.1f #pm %.1f) * thresh [mV]" % (params[0],errors[0],params[1],errors[1]))
#p.SetFillStyle(0)
#p.SetBorderSize(0)
#p.Draw()

c2 = ROOT.TCanvas("c2","c2",800,600)
g = ROOT.TGraph(len(thresholds),np.array(measured_thresh[1:]),np.array(widths))
g.SetMarkerStyle(21)
#fr = mg.Fit("pol1","SQ")
g.SetTitle("Fe55 Threshold width")
g.Draw("ALP")
g.GetXaxis().SetTitle("Set threshold (mV)")
g.GetYaxis().SetTitle("Signal cutoff width (ADC counts)")
#params2 = fr.GetParams()
#errors2 = fr.GetErrors()
#print "Cutoff width (counts): (",params2[0],"+-",errors2[0],") + (",params2[1],"+-",errors2[1],") * thresh [mV]"
#p = ROOT.TPaveText(0.5,0.2,0.7,0.3,"NDC NB")
#p.AddText("Cutoff width (counts) = (%.1f #pm %.1f) + (%.1f #pm %.1f) * thresh [mV]" % (params2[0],errors2[0],params2[1],errors2[1]))
#p.SetFillStyle(0)
#p.SetBorderSize(0)
#p.Draw()

fout = ROOT.TFile(sys.argv[2],"RECREATE")
t_measured_thresh = array.array('d',[0])
t_thresh = array.array('d',[0])
t_width = array.array('d',[0])
tout = ROOT.TTree("fe55threshold","fe55threshold")
tout.Branch("measured_thresh",t_measured_thresh,"measured_thresh/D")
tout.Branch("thresh",t_thresh,"thresh/D")
tout.Branch("width",t_width,"width/D")
for i in range(1,len(measured_thresh)):
  t_measured_thresh[0] = measured_thresh[i]
  t_thresh[0] = thresholds[i-1]
  t_width[0] = widths[i-1]
  tout.Fill()
tout.Write()

#threshoffset = array.array('d',[params[0]])
#threshslope = array.array('d',[params[1]])
#widthoffset = array.array('d',[params2[0]])
#widthslope = array.array('d',[params2[1]])
#threshoffseterr = array.array('d',[errors[0]])
#threshslopeerr = array.array('d',[errors[1]])
#widthoffseterr = array.array('d',[errors2[0]])
#widthslopeerr = array.array('d',[errors2[1]])
#
#tout = ROOT.TTree("fe55thresholdfit","fe55thresholdfit")
#tout.Branch("threshoffset",threshoffset,"threshoffset/D")
#tout.Branch("threshslope",threshslope,"threshslope/D")
#tout.Branch("widthoffset",widthoffset,"widthoffset/D")
#tout.Branch("widthslope",widthslope,"widthslope/D")
#tout.Branch("threshoffseterr",threshoffseterr,"threshoffseterr/D")
#tout.Branch("threshslopeerr",threshslopeerr,"threshslopeerr/D")
#tout.Branch("widthoffseterr",widthoffseterr,"widthoffseterr/D")
#tout.Branch("widthslopeerr",widthslopeerr,"widthslopeerr/D")
#tout.Fill()
#tout.Write()

raw_input()
