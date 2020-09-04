import ROOT
import array
import numpy as np
import sys

# input file should have CosmicTrackDiag TTree
if len(sys.argv) < 2:
  print("Usage: recalibrate.py <CosmicTrackDiag TTree filename> <input fcl filename> <output fcl filename>  [use unbiased?]")
  sys.exit()

unbiased = False
if len(sys.argv) > 4:
  unbiased = True

f = ROOT.TFile(sys.argv[1])
t = f.Get("CosmicTrackDiag/hitT")

fold = open(sys.argv[2])
lold = fold.readlines()
old_fcl = {line.split(":")[0].rstrip() : line.split(":")[1].rstrip() for line in lold}
old_timeOffsetHV = eval(old_fcl["services.ProditionsService.strawResponse.timeOffsetStrawHV"])
old_timeOffsetCal = eval(old_fcl["services.ProditionsService.strawResponse.timeOffsetStrawCal"])

hresids = ROOT.TH1F("hresids","hresids",96,0,96)
hpvs = ROOT.TH1F("hpvs","hpvs",96,0,96)

habsresids = ROOT.TProfile("habsresids","habsresids",96,0,96)
if unbiased:
  t.Project("habsresids","ubtresid:straw","long != -999 && abs(long) < slen*1.1 && ublong != -999 && abs(ublong) < slen*1.1 && abs(deltat) < 20 && doca < 10 && hitused && ubdoca < 10 && abs(ubtresid) < 100")
else:
  t.Project("habsresids","tresid:straw","long != -999 && abs(long) < slen*1.1 && ublong != -999 && abs(ublong) < slen*1.1 && abs(deltat) < 20 && doca < 10 && hitused && ubdoca < 10 && abs(ubtresid) < 100")

perStrawHPV = []
timeOffsetCal = []
timeOffsetHV = []
hp = ROOT.TProfile("hp","hp",100,-500,500)
for i in range(96):
  if unbiased:
    t.Project("hp","deltat:ublong","long != -999 && abs(long) < slen*1.1 && ublong != -999 && abs(ublong) < slen*1.1 && abs(deltat) < 20 && straw == %d" % i)
  else:
    t.Project("hp","deltat:long","long != -999 && abs(long) < slen*1.1 && abs(deltat) < 20 && straw == %d" % i)
  
  try:
    hp.Approximate()
    fr = hp.Fit("pol1","SQN") 
    perStrawHPV.append(1.0/fr.GetParams()[1])
    timeOffsetCal.append(0)
    timeOffsetHV.append(fr.GetParams()[0])
  except:
    perStrawHPV.append(100)
    timeOffsetCal.append(0)
    timeOffsetHV.append(0)
    continue

  if unbiased:
    t.Project("hp","recolong:ublong","long != -999 && abs(long) < slen*1.1 && ublong != -999 && abs(ublong) < slen*1.1 && abs(deltat) < 20 && straw == %d" % i)
  else:
    t.Project("hp","recolong:long","long != -999 && abs(long) < slen*1.1 && abs(deltat) < 20 && straw == %d" % i)
  
  try:
    hp.Approximate()
    fr = hp.Fit("pol1","SQN") 
    hresids.SetBinContent(i+1,fr.GetParams()[0]/perStrawHPV[-1])
    hresids.SetBinError(i+1,fr.GetErrors()[0]/perStrawHPV[-1])
#    hresids.SetBinContent(i+1,fr.GetParams()[0])
#    hresids.SetBinError(i+1,fr.GetErrors()[0])
    hpvs.SetBinContent(i+1,fr.GetParams()[1]-1)
    hpvs.SetBinError(i+1,fr.GetErrors()[1])
  except:
    continue

if unbiased:
  t.Project("hp","deltat:ublong","long != -999 && abs(long) < slen*1.1 && ublong != -999 && abs(ublong) < slen*1.1 && abs(deltat) < 20")
else:
  t.Project("hp","deltat:long","long != -999 && abs(long) < slen*1.1 && abs(deltat) < 20")
  
hp.Approximate()
fr = hp.Fit("pol1","SQN") 
hpv = 1.0/fr.GetParams()[1]



#print("perStrawHPV:",perStrawHPV)
#print("timeOffsetStrawCal:",timeOffsetCal)
#print("timeOffsetStrawHV:",[timeOffsetHV[i] + oldtimehv[i] for i in range(96)])

fout = open(sys.argv[3],"w")
unchanged_params = ["eDep","halfPropVelocity","centralWirePos","tdCentralRes","tdResSlope"]
for param in unchanged_params:
  fout.write("services.ProditionsService.strawResponse.%s : %s\n" % (param,old_fcl["services.ProditionsService.strawResponse.%s" % param]))
fout.write("services.ProditionsService.strawResponse.timeOffsetStrawCal : [")
for i in range(len(timeOffsetCal)):
  fout.write(" %f" % (timeOffsetCal[i] + old_timeOffsetCal[i] + habsresids.GetBinContent(i+1)))
  if i != len(timeOffsetCal)-1:
    fout.write(",")
  else:
    fout.write("]\n")
fout.write("services.ProditionsService.strawResponse.timeOffsetStrawHV : [")
for i in range(len(timeOffsetHV)):
  fout.write(" %f" % (timeOffsetHV[i] + old_timeOffsetHV[i] + habsresids.GetBinContent(i+1)))
  if i != len(timeOffsetHV)-1:
    fout.write(",")
  else:
    fout.write("]\n")
fout.close()

hresids.GetXaxis().SetTitle("Straw number")
hresids.GetYaxis().SetTitle("Longitudinal offset (mm)")
hresids.GetYaxis().SetRangeUser(-1,1)
#hpvs.GetXaxis().SetTitle("Straw number")
#hpvs.GetYaxis().SetTitle("Propagational velocity correction (fractional)")
#hpvs.GetYaxis().SetRangeUser(-0.1,0.1)
c1 = ROOT.TCanvas("c1","c1",600,600)
hresids.Draw()
#c2 = ROOT.TCanvas("c2","c2",600,600)
#hpvs.Draw()
c3 = ROOT.TCanvas("c3","c3",600,600)
if unbiased:
  t.Draw("ubtresid","long != -999 && abs(long) < slen*1.1 && ublong != -999 && abs(ublong) < slen*1.1 && abs(deltat) < 20 && doca < 10 && hitused && ubdoca < 10 && abs(ubtresid) < 100")
else:
  t.Draw("tresid","long != -999 && abs(long) < slen*1.1 && ublong != -999 && abs(ublong) < slen*1.1 && abs(deltat) < 20 && doca < 10 && hitused && ubdoca < 10 && abs(ubtresid) < 100")
c4 = ROOT.TCanvas("c4","c4",600,600)
habsresids.Draw()

input()
