import ROOT
import array
import numpy as np
import sys

def read_tables(fn):
  f = open(fn)
  lines = f.readlines()
  tables = {}

  current_tablename = None
  current_table = []
  for line in lines:
    line = line.strip()
    if len(line) == 0 or line[0] == "#":
      continue
    if line.startswith("TABLE"):
      if current_tablename:
        tables[current_tablename] = current_table
      current_tablename = line.split()[1]
      current_table = []
    else:
      if not current_tablename:
        print("ERROR: data without a table?")
        continue
      current_table.append(list(map(float,line.split(","))))
  if current_tablename:
    tables[current_tablename] = current_table
  return tables


# input file should have CosmicTrackDiag TTree
if len(sys.argv) < 2:
  print("Usage: recalibrate.py <CosmicTrackDiag TTree filename> <input table filename> <output table filename>  [use unbiased?]")
  sys.exit()

unbiased = False
if len(sys.argv) > 4:
  unbiased = True

f = ROOT.TFile(sys.argv[1])
t = f.Get("CosmicTrackDiag/hitT")

old_tables = read_tables(sys.argv[2])

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
fout.write("TABLE TrkPreampRStraw\n")
fout.write("# index,delay_hv,delay_cal,threshold_hv,threshold_cal,gain\n")
for i in range(len(timeOffsetHV)):
  index = old_tables["TrkPreampRStraw"][i][0]
  toHV = timeOffsetHV[i] + habsresids.GetBinContent(i+1) + old_tables["TrkPreampRStraw"][i][1]
  toCal = timeOffsetCal[i] + habsresids.GetBinContent(i+1) + old_tables["TrkPreampRStraw"][i][2]
  threshHV = old_tables["TrkPreampRStraw"][i][3]
  threshCal = old_tables["TrkPreampRStraw"][i][4]
  gain = old_tables["TrkPreampRStraw"][i][5]
  fout.write("%d,%f,%f,%f,%f,%f\n" % (index,toHV,toCal,threshHV,threshCal,gain))

# FIXME now copy other tables
for table in old_tables:
  if table == "TrkPreampRStraw":
    continue
  fout.write("\n")
  fout.write("TABLE %s\n" % table)
  for i in range(len(old_tables[table])):
    fout.write("%d,%s\n" % (int(old_tables[table][i][0]),",".join(map(str,old_tables[table][i][1:]))))
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
