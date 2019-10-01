import ROOT
import numpy as np

ROOT.gSystem.Load("event/Dict.so")

filenames = ["run_2683_2684.root","run_2681.root","run_2682.root","run_2685.root","run_2686.root","run_2687.root","run_2688.root","run_2689.root","run_2690.root"]
#fn_peds = {}

def calc_threshold(fn,straws):
  f = ROOT.TFile("/data/HD2/rbonventre/prototypeData/lblPrototype/%s" % fn)
  t = f.Get("T")
  pedestals = {straw: [] for straw in straws}
  for i in range(t.GetEntries()):
    t.GetEntry(i)
    for j in range(t.events.straws.size()):
      for straw in straws:
        if t.events.straws[j].channel == straw:
          pedestals[straw].append(t.events.straws[j].pedestal)
  threshold = {straw: (4096/2.-np.median(pedestals[straw]))/4096.*1500./20.*1.1-0.4 for straw in straws}
  print fn,threshold
  return threshold

#print fn_peds
#straws = [3,5]
straws = [5]
#fn_peds = {fn: calc_threshold(fn,straws) for fn in filenames}
fn_peds = {'run_2682.root': {3: 12.376490077376367, 5: 12.416773280501367},
           'run_2683_2684.root': {3: 9.9796394914388653, 5: 10.0132080078125},
           'run_2685.root': {3: 15.001610508561136, 5: 14.988183593750001},
           'run_2686.root': {3: 17.405175781250001, 5: 17.358177891373636},
           'run_2687.root': {3: 17.438744297623636, 5: 17.385034179687501},
           'run_2688.root': {3: 19.996727696061136, 5: 19.909449061751367},
           'run_2689.root': {3: 22.393579511344438, 5: 22.313013105094438},
           'run_2690.root': {3: 24.770288495719438, 5: 24.676293945312505}}

fn_pmpnorms = {'run_2682.root': {5: 623.168},
               'run_2683_2684.root' : {5: 624.364},
               'run_2685.root' : {5: 631.567},
               'run_2686.root' : {5: 631.315},
               'run_2687.root' : {5: 638.598},
               'run_2688.root' : {5: 643.704},
               'run_2689.root' : {5: 645.513},
               'run_2690.root' : {5: 644.315}}

thresh = {straw: [] for straw in straws}
width = {straw: [] for straw in straws}
measured_thresh = {straw: [] for straw in straws}

draw = False

# calculate the run time for the lowest threshold runs (2683 and 2684)
# which will be used to scale the other runs
baselength = 0
f = ROOT.TFile("/data/HD2/rbonventre/prototypeData/lblPrototype/run_2683.root")
t = f.Get("T") 
t.GetEntry(0)
starttime = t.events.straws[0].timeGlobal
t.GetEntry(t.GetEntries()-1)
endtime = t.events.straws[0].timeGlobal
baselength += ((endtime-starttime)*2**28*15.625*10**-12)
f.Close()
f = ROOT.TFile("/data/HD2/rbonventre/prototypeData/lblPrototype/run_2684.root")
t = f.Get("T") 
t.GetEntry(0)
starttime = t.events.straws[0].timeGlobal
t.GetEntry(t.GetEntries()-1)
endtime = t.events.straws[0].timeGlobal
baselength += ((endtime-starttime)*2**28*15.625*10**-12)
f.Close()


lengths = []
norms = []
hs = {straw: [] for straw in straws}
for key, value in sorted(fn_peds.iteritems(), key=lambda (k,v): (v[5],k)):
  print "%s: %s" % (key, value)
  f = ROOT.TFile("/data/HD2/rbonventre/prototypeData/lblPrototype/%s" % key)
  t = f.Get("T") 
  t.GetEntry(0)
  starttime = t.events.straws[0].timeGlobal
  t.GetEntry(t.GetEntries()-1)
  endtime = t.events.straws[0].timeGlobal
  lengths.append((endtime-starttime)*2**28*15.625*10**-12)
  norms.append(fn_pmpnorms[key][5]/fn_pmpnorms["run_2683_2684.root"][5])

  for straw in straws:
    hs[straw].append(ROOT.TH1F("h_%s_%d" % (key,straw),"h_%s" % key,100,125,525))
    hs[straw][-1].Sumw2()
    t.Project("h_%s_%d" % (key,straw),"events.straws.peak-events.straws.pedestal","events.straws.channel == %d" % straw)
    hs[straw][-1].SetDirectory(0)
  f.Close()
  #print hs[straw][-1].Integral(),lengths[-1]
  #print hs[3][-1].Integral(),hs[5][-1].Integral(),lengths[-1]

lengths[0] = baselength

import sys
import array
import scipy.stats
from scipy.optimize import curve_fit
if draw:
  import matplotlib.pyplot as plt

def fcn(x,mu,sig):
  return scipy.stats.norm(loc=mu,scale=sig).cdf(x)

for straw in straws:
  ndt = hs[straw][0].Integral() / lengths[0]
  for i in range(len(hs[straw])):
    hs[straw][i].Scale(1.0/(ndt*lengths[i]))
    hs[straw][i].SetLineColor(i+1)
    if i==0:
      pass
  #    hs[i].Draw()
    else:
      hs[straw][i].Divide(hs[straw][0])
  #    hs[i].Draw("same")
      x = np.zeros(hs[straw][i].GetNbinsX())
      y = np.zeros(hs[straw][i].GetNbinsX())
      for j in range(hs[straw][i].GetNbinsX()):
        x[j] = hs[straw][i].GetXaxis().GetBinCenter(j+1)
        y[j] = hs[straw][i].GetBinContent(j+1)
  
      popt, pcov = curve_fit(fcn,x,y,p0=[200,10])
      yfit = scipy.stats.norm(loc=popt[0],scale=popt[1]).cdf(x)
  
      print popt
      print pcov
      print popt[1]/popt[0]
      thresh[straw].append(popt[0]/norms[i])
      width[straw].append(popt[1]/norms[i])
  
      if draw:
        plt.plot(x,y)
        plt.plot(x,yfit)
        plt.xlabel("Straw %d: %s" % (straw,sorted(fn_peds.iteritems(), key=lambda (k,v): (v[5],k))[i][0]))
        plt.show()

for key, value in sorted(fn_peds.iteritems(), key=lambda (k,v): (v[5],k)):
  for straw in straws:
    measured_thresh[straw].append(value[straw])

c1 = ROOT.TCanvas("c1","c1",800,600)
mg = ROOT.TMultiGraph("mg","mg")
l = ROOT.TLegend(0.15,0.65,0.35,0.85)
color_index = 3
for straw in straws:
  g = ROOT.TGraph(len(thresh[straw]),np.array(measured_thresh[straw][1:]),np.array(thresh[straw]))
  g.SetMarkerStyle(21)
  g.SetMarkerColor(color_index)
  g.SetLineColor(color_index)
  color_index += 1
  l.AddEntry(g,"Straw %d" % straw,"p")
  mg.Add(g)
 
fr = mg.Fit("pol1","SQ")
mg.SetTitle("Fe55 Threshold")
mg.Draw("ALP")
l.Draw()
mg.GetXaxis().SetTitle("Set threshold (mV)")
mg.GetYaxis().SetTitle("Output signal threshold (peak minus pedestal counts)")
params = fr.GetParams()
errors = fr.GetErrors()
print "Signal threshold (counts pmp): (",params[0],"+-",errors[0],") + (",params[1],"+-",errors[1],") * thresh [mV]"
p = ROOT.TPaveText(0.5,0.2,0.7,0.3,"NDC NB")
p.AddText("signal cutoff (pmp) = (%.1f #pm %.1f) + (%.1f #pm %.1f) * thresh [mV]" % (params[0],errors[0],params[1],errors[1]))
p.SetFillStyle(0)
p.SetBorderSize(0)
p.Draw()

c2 = ROOT.TCanvas("c2","c2",800,600)
mg2 = ROOT.TMultiGraph("mg2","mg2")
l2 = ROOT.TLegend(0.15,0.65,0.35,0.85)
color_index = 3
for straw in straws:
  g2 = ROOT.TGraph(len(thresh[straw]),np.array(measured_thresh[straw][1:]),np.array(width[straw]))
  g2.SetMarkerStyle(21)
  g2.SetMarkerColor(color_index)
  g2.SetLineColor(color_index)
  color_index += 1
  l2.AddEntry(g2,"Straw %d" % straw,"p")
  mg2.Add(g2)
 
#  plt.plot(measured_thresh[straw][1:],width[straw],label="Straw %d" % straw)
#plt.legend()
#plt.xlabel("measured threshold (mV)")
#plt.ylabel("fit width (pmp)")
#plt.show()
fr2 = mg2.Fit("pol1","SQ")
mg2.SetTitle("Fe55 Threshold width")
mg2.Draw("ALP")
l2.Draw()
mg2.GetXaxis().SetTitle("Set threshold (mV)")
mg2.GetYaxis().SetTitle("Signal cutoff width (ADC counts)")
params2 = fr2.GetParams()
errors2 = fr2.GetErrors()
print "Cutoff width (counts): (",params2[0],"+-",errors2[0],") + (",params2[1],"+-",errors2[1],") * thresh [mV]"
p2 = ROOT.TPaveText(0.5,0.2,0.7,0.3,"NDC NB")
p2.AddText("Cutoff width (counts) = (%.1f #pm %.1f) + (%.1f #pm %.1f) * thresh [mV]" % (params2[0],errors2[0],params2[1],errors2[1]))
p2.SetFillStyle(0)
p2.SetBorderSize(0)
p2.Draw()

fout = ROOT.TFile(sys.argv[1],"RECREATE")
t_straw = array.array('i',[0])
t_measured_thresh = array.array('d',[0])
t_thresh = array.array('d',[0])
t_width = array.array('d',[0])
tout = ROOT.TTree("fe55threshold","fe55threshold")
tout.Branch("straw",t_straw,"straw/I")
tout.Branch("measured_thresh",t_measured_thresh,"measured_thresh/D")
tout.Branch("thresh",t_thresh,"thresh/D")
tout.Branch("width",t_width,"width/D")
for straw in straws:
  for i in range(1,len(measured_thresh[straw])):
    t_straw[0] = straw
    t_measured_thresh[0] = measured_thresh[straw][i]
    t_thresh[0] = thresh[straw][i-1]
    t_width[0] = width[straw][i-1]
    tout.Fill()
tout.Write()

threshoffset = array.array('d',[params[0]])
threshslope = array.array('d',[params[1]])
widthoffset = array.array('d',[params2[0]])
widthslope = array.array('d',[params2[1]])
threshoffseterr = array.array('d',[errors[0]])
threshslopeerr = array.array('d',[errors[1]])
widthoffseterr = array.array('d',[errors2[0]])
widthslopeerr = array.array('d',[errors2[1]])

tout = ROOT.TTree("fe55thresholdfit","fe55thresholdfit")
tout.Branch("threshoffset",threshoffset,"threshoffset/D")
tout.Branch("threshslope",threshslope,"threshslope/D")
tout.Branch("widthoffset",widthoffset,"widthoffset/D")
tout.Branch("widthslope",widthslope,"widthslope/D")
tout.Branch("threshoffseterr",threshoffseterr,"threshoffseterr/D")
tout.Branch("threshslopeerr",threshslopeerr,"threshslopeerr/D")
tout.Branch("widthoffseterr",widthoffseterr,"widthoffseterr/D")
tout.Branch("widthslopeerr",widthslopeerr,"widthslopeerr/D")
tout.Fill()
tout.Write()

raw_input()
