import ROOT
import numpy as np
import sys
import array

if len(sys.argv) < 2:
  print "usage: bin/fe55adc_sim.py <input filename> <output filename>"

filename = sys.argv[1]

pmps = []
f = ROOT.TFile(filename)
t = f.Get("makeSD/sddiag")
for i in range(t.GetEntries()):
  t.GetEntry(i)
#  if t.panel == 1 and t.straw > 88:
  if 1:
    adcs = [t.adc[j] for j in range(15)]
    pmps.append(np.max(adcs)-np.mean(adcs[:2]))

h = ROOT.TH1F("h","h",100,0,np.mean(pmps)*1.5)
for i in range(len(pmps)):
  h.Fill(pmps[i])

c = ROOT.TCanvas("c1","c1",600,800)
h.Draw()
f1 = ROOT.TF1("f1","gaus(0)+gaus(3)+pol1(6)")
f1.SetParameter(0,h.GetMaximum())
f1.SetParameter(1,h.GetBinCenter(h.GetMaximumBin()))
f1.SetParameter(2,h.GetRMS())
f1.SetParameter(3,h.GetMaximum()/2.)
f1.SetParameter(4,h.GetBinCenter(h.GetMaximumBin()/2))
f1.SetParameter(5,h.GetRMS())
test = 0
while True:
  test2 = raw_input("Min val: (%f)" % test)
  if len(test2) > 0:
    test = float(test2)
  min_fit = test
  h.Fit("f1","","",min_fit,np.mean(pmps)*1.5)
  c.Update()
  test3 = raw_input("Redo?")
  if len(test3) == 0:
    break
  for i in range(8):
    val = raw_input("Parameter %d: (%f)" % (i,f1.GetParameter(i)))
    if len(val) > 0:
      f1.SetParameter(i,float(val))
print "DONE"

fout = ROOT.TFile(sys.argv[2],"RECREATE")
pmp = array.array('d',[0])
tout = ROOT.TTree("fe55adc","fe55adc")
tout.Branch("pmp",pmp,"pmp/D")
for i in range(len(pmps)):
  pmp[0] = pmps[i]
  tout.Fill()
tout.Write()

peak0scale = array.array('d',[f1.GetParameter(0)])
peak0mean = array.array('d',[f1.GetParameter(1)])
peak0width = array.array('d',[f1.GetParameter(2)])
peak1scale = array.array('d',[f1.GetParameter(3)])
peak1mean = array.array('d',[f1.GetParameter(4)])
peak1width = array.array('d',[f1.GetParameter(5)])
flatscale = array.array('d',[f1.GetParameter(6)])
flatslope = array.array('d',[f1.GetParameter(7)])
peak0scaleerr = array.array('d',[f1.GetParError(0)])
peak0meanerr = array.array('d',[f1.GetParError(1)])
peak0widtherr = array.array('d',[f1.GetParError(2)])
peak1scaleerr = array.array('d',[f1.GetParError(3)])
peak1meanerr = array.array('d',[f1.GetParError(4)])
peak1widtherr = array.array('d',[f1.GetParError(5)])
flatscaleerr = array.array('d',[f1.GetParError(6)])
flatslopeerr = array.array('d',[f1.GetParError(7)])

tout = ROOT.TTree("fe55adcfit","fe55adcfit")
tout.Branch("peak0scale",peak0scale,"peak0scale/D")
tout.Branch("peak0mean",peak0mean,"peak0mean/D")
tout.Branch("peak0width",peak0width,"peak0width/D")
tout.Branch("peak1scale",peak1scale,"peak1scale/D")
tout.Branch("peak1mean",peak1mean,"peak1mean/D")
tout.Branch("peak1width",peak1width,"peak1width/D")
tout.Branch("flatscale",flatscale,"flatscale/D")
tout.Branch("flatslope",flatslope,"flatslope/D")
tout.Branch("peak0scaleerr",peak0scaleerr,"peak0scaleerr/D")
tout.Branch("peak0meanerr",peak0meanerr,"peak0meanerr/D")
tout.Branch("peak0widtherr",peak0widtherr,"peak0widtherr/D")
tout.Branch("peak1scaleerr",peak1scaleerr,"peak1scaleerr/D")
tout.Branch("peak1meanerr",peak1meanerr,"peak1meanerr/D")
tout.Branch("peak1widtherr",peak1widtherr,"peak1widtherr/D")
tout.Branch("flatscaleerr",flatscaleerr,"flatscaleerr/D")
tout.Branch("flatslopeerr",flatslopeerr,"flatslopeerr/D")
tout.Fill()
tout.Write()
