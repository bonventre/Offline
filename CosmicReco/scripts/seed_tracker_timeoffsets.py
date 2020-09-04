import ROOT
import numpy as np
import sys

# input file should have StrawHitDiag TTree
if len(sys.argv) < 2:
  print("Usage: seed_tracker_timeoffsets.py <StrawHitDiag TTree filename> [plot?]")
  sys.exit()

fn = sys.argv[1]
plot = False
if len(sys.argv) > 2:
  plot = True

f = ROOT.TFile(fn)
t = f.Get("SHD/shdiag")
h = ROOT.TH1F("h","h",96,0,96)
t.Project("h","straw","dead == 0")

# get closest functional neighbor to each straw
pairs = {}
strawi = 0
while strawi < 95:
  if h.GetBinContent(strawi+1) == 0:
    strawi += 1
    continue
  strawj = strawi + 1
  while strawj < 95:
    if h.GetBinContent(strawj+1) == 0:
      strawj += 1
      continue
    break
  if h.GetBinContent(strawj+1) == 0:
    break
  pairs[strawi] = strawj
  strawi += 1

# get delta t distribution for each channel to find center of each straw
print("Calibrating delta t")
hwidths = ROOT.TH1F("hwidths","hwidths",96,0,96)
ha = ROOT.TH1F("ha","ha",250,-0.024414625*1000,0.024414625*1000)
lowedges = []
highedges = []
for i in range(96):
  ha.Scale(0)
  t.Project("ha","tcal-thv","straw == %d && dead == 0 && esel == 1" % i)
  if ha.Integral() == 0:
    lowedges.append(0)
    highedges.append(0)
    continue

  maxval = ha.GetMaximum()
  minbin = 1
  while True:
    if ha.GetBinContent(minbin) > maxval/10.:
      break
    minbin += 1
  maxbin = ha.GetNbinsX()
  while True:
    if ha.GetBinContent(maxbin) > maxval/10.:
      break
    maxbin -= 1

  mindeltat = ha.GetBinLowEdge(minbin)
  maxdeltat = ha.GetBinLowEdge(maxbin+1)

  lowedges.append(mindeltat)
  highedges.append(maxdeltat)

  hwidths.SetBinContent(i+1,(maxdeltat-mindeltat))

  if plot:
    ha.Draw("hist")
    l0 = ROOT.TLine(mindeltat,0,mindeltat,ha.GetMaximum()*1.1)
    l1 = ROOT.TLine(maxdeltat,0,maxdeltat,ha.GetMaximum()*1.1)
    l0.Draw()
    l1.Draw()




# get time between straw hits for all events with coincidences between neighbors
# this distribution will average over drift times
print("Calibrating straw to straw time")
hs = [ROOT.TH1F("h_%d" % i,"h_%d" % i,200,-100,100) for i in range(96)]
last_event = -1
straws = []
times = []
failed = False
for i in range(t.GetEntries()):
    if i % 100000 == 0:
      print(i,t.GetEntries())
    t.GetEntry(i)
    if t.eventid != last_event and last_event != -1:
      if not failed:
        for indexi in range(len(straws)):
          if not straws[indexi] in pairs:
            continue
          try:
            indexj = straws.index(pairs[straws[indexi]])
          except:
            continue
          straw_to_straw_deltat = times[indexi] - times[indexj]
          hs[straws[indexi]].Fill(straw_to_straw_deltat)
      straws = []
      times = []
      failed = False
    last_event = t.eventid
    if t.esel == 0:
      failed = True
      continue
    straws.append(t.straw)
    times.append((t.tcal + t.thv)/2.)

offsets = [0 for i in range(96)]
for j in range(96):
    if not j in pairs:
        continue
    if hs[j].Integral() > 100:
        offsets[j] = hs[j].GetMean()
        if plot:
          print(j,pairs[j],hs[j].GetMean(),hs[j].GetRMS())
          hs[j].Draw()
          input()
    else:
        print(j,pairs[j],"NOT ENOUGH HITS")

total_offsets = [0 for i in range(96)]
for i in range(95):
  if not i in pairs:
    continue
  total_offsets[pairs[i]] = total_offsets[i] + offsets[i]

timeOffsetCal = []
timeOffsetHV = []
for i in range(96):
  center = (highedges[i]+lowedges[i])/2. # positive center means hv time offset is positive
  timeOffsetCal.append(total_offsets[i])
  timeOffsetHV.append(total_offsets[i] + center)


fout = open("seed_table.txt","w")
fout.write("TABLE TrkPreampRStraw\n")
fout.write("# index,delay_hv,delay_cal,threshold_hv,threshold_cal,gain\n")
for i in range(len(timeOffsetHV)):
  fout.write("%d,%f,%f,0.0,0.0,0.0\n" % (i,timeOffsetHV[i],timeOffsetCal[i]))

# FIXME now write other tables in order to make sure useDb doesn't fail due to run/subrun validity
fout.write("\n")
fout.write("TABLE TrkDelayPanel\n")
for i in range(216):
  fout.write("%d,0.0\n" % (i))
fout.write("\n")
fout.write("TABLE TrkPreampStraw\n")
for i in range(216*96):
  fout.write("%d,0.0,0.0,0.0,0.0,0.0\n" % (i))
fout.write("\n")
fout.write("TABLE TrkThresholdRStraw\n")
for i in range(96):
  fout.write("%d,12.0,12.0\n" % (i))
fout.close()
