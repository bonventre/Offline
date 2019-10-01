#include <iostream>
#include <TFile.h>
#include <TTree.h>

#include "prototype/Parser.hh"

bool has_non_neighbor_straw(Event *event, int channel)
{
  for (int i=0;i<event->straws.size();i++){
    if (abs(event->straws[i].channel - channel) > 1)
      return true;
  }
  return false;
}

bool has_straw_shower(Event *event)
{
  for (int i=0;i<event->straws.size();i++){
    for (int j=i+1;j<event->straws.size();j++){
      if (abs(event->straws[i].channel - event->straws[j].channel) > 1)
        return true;
    }
  }
  return false;
}

bool valid_pmt_data(Event *event, long &pmtTopTime, long &pmtBotTime, long &pmtTrigTime, long &lastPmtTopTime, long &lastPmtTrigTime)
{
  pmtTopTime = (event->pmtTop.timeGlobal * 268435456 + event->pmtTop.timeCal);
  pmtBotTime = (event->pmtTop.timeGlobal * 268435456 + event->pmtTop.timeHV);
  ///pmtBotTime = (event->pmtBot.timeGlobal * 268435456 + event->pmtBot.timeCal);
  pmtTrigTime = (event->pmtTrig.timeGlobal * 268435456 + event->pmtTrig.timeCal);
  //if (pmtTopTime == lastPmtTopTime || pmtTrigTime == lastPmtTrigTime || event->pmtTop.channel == -1){
  if (pmtTopTime == lastPmtTopTime || pmtTrigTime == lastPmtTrigTime){
    lastPmtTopTime = pmtTopTime;
    lastPmtTrigTime = pmtTrigTime;
    return false;
  }
  lastPmtTopTime = pmtTopTime;
  lastPmtTrigTime = pmtTrigTime;
  return true;
}


bool valid_pixel_data(Event *event, double &topavg_x, double &topavg_y, double &botavg_x, double &botavg_y, bool old_data)
{
  //FIXME bcid, tot
  //FIXME check if columns are adjacent!
  topavg_x = topavg_y = botavg_x = botavg_y = 0;
  if (event->topPixel.size() == 0 || event->botPixel.size() == 0)
    return false;

  // top and bot are no longer flipped!
  for (int i=0;i<event->botPixel.size();i++){
    botavg_x += event->botPixel[i].row;
    botavg_y += event->botPixel[i].column;
    if (i > 0){
      if (abs(event->botPixel[i].row - event->botPixel[i-1].row) > 1)
        return false;
    }
  }
  botavg_x /= event->botPixel.size();
  botavg_y /= event->botPixel.size();

  for (int i=0;i<event->topPixel.size();i++){
    topavg_x += event->topPixel[i].row;
    topavg_y += event->topPixel[i].column;
    if (i > 0){
      if (abs(event->topPixel[i].row - event->topPixel[i-1].row) > 1)
        return false;
    }
  }
  topavg_x /= event->topPixel.size();
  topavg_y /= event->topPixel.size();

  // FLIP OLD DATA
  if (old_data){
    double temp_x = topavg_x;
    double temp_y = topavg_y;
    topavg_x = botavg_x;
    topavg_y = botavg_y;
    botavg_x = temp_x;
    botavg_y = temp_y;
  }

  return true;
}

bool valid_tot_data(Event *event, int channel, double &tothv, double &totcal)
{
  if (event->straws.size() > 0){
    for (int i=0;i<event->straws.size();i++){
      if (channel == event->straws[i].channel){
        tothv = event->straws[i].totHV * 4 + (0xFF - (event->straws[i].timeHV & 0xFF))*0.015625;
        totcal = event->straws[i].totCal * 4 + (0xFF - (event->straws[i].timeCal & 0xFF))*0.015625;
        return true;
      }
    }
  }
  return false;
}

bool valid_straw_data(Event *event, int channel, long &strawtimecal, long &strawtimehv, double &maxval)
{
  if (event->straws.size() > 0){
    for (int i=0;i<event->straws.size();i++){
      if (channel == event->straws[i].channel){
        strawtimecal = (event->straws[i].timeGlobal * 268435456 + event->straws[i].timeCal);
        strawtimehv = (event->straws[i].timeGlobal * 268435456 + event->straws[i].timeHV);
        double pedestal = (event->straws[i].samples[0] + event->straws[i].samples[1])/2.;
        maxval = 0;
        for (int j=2;j<event->straws[i].samples.size();j++)
          if (event->straws[i].samples[j] > maxval)
            maxval = event->straws[i].samples[j];
        maxval -= pedestal;
//        maxval = event->straws[i].peak - event->straws[i].pedestal;
        return true;
      }
    }
  }
  return false;
}

void parse_data_files(std::vector<std::string> &filenames, int strawnum, std::vector<double> &alltop_x, std::vector<double> &allbot_x, std::vector<double> &alltop_y, std::vector<double> &allbot_y, std::vector<std::vector<bool> > &has_straw, std::vector<std::vector<bool> > &has_shower, std::vector<std::vector<double> > &alltimes, std::vector<std::vector<double> > &alldts, std::vector<std::vector<double> > &allmaxvals, std::vector<std::vector<double> > &alltots, bool old_data, double pmpminimum, int use_pmt_time){

  double total_time = 0; // total time length of these run before HV trips


  Event *event = new Event();

  for (int i=0;i<filenames.size();i++){
//    std::cout << filenames[i] << std::endl;
    TFile *f = new TFile(filenames[i].c_str());
    TTree *t = (TTree *) f->Get("T");
    t->GetBranch("events")->SetAddress(&event);

    // At some point the HV trips and we get no more straw hits, determine when that was
    int last_event_with_straw = 0;
    double lasttime = 0;
    for (int j=0;j<t->GetEntries();j++){
      t->GetEntry(j);

      if (event->straws.size() > 0){
        last_event_with_straw = j;
        lasttime = event->straws[0].timeGlobal*268435456*15.625e-12;
      }
    }
    total_time += lasttime;


    long lastPmtTopTimeCount = -1, lastPmtTrigTimeCount = -1;
    for (int j=0;j<last_event_with_straw;j++){
      t->GetEntry(j);

      // check if pmt time is valid
      long pmtTopTimeCount, pmtBotTimeCount, pmtTrigTimeCount;
      if (!valid_pmt_data(event,pmtTopTimeCount,pmtBotTimeCount,pmtTrigTimeCount,lastPmtTopTimeCount,lastPmtTrigTimeCount))
        continue;

      // check if shower
      if (has_straw_shower(event))
        continue;

      // check if pixel data is valid
      double topavg_x,topavg_y,botavg_x,botavg_y;
      if (!valid_pixel_data(event,topavg_x,topavg_y,botavg_x,botavg_y,old_data))
        continue;

      alltop_x.push_back(topavg_x);
      alltop_y.push_back(topavg_y);
      allbot_x.push_back(botavg_x);
      allbot_y.push_back(botavg_y);

      for (int k=0;k<8;k++){
        if (strawnum == k || strawnum < 0){
          // check if hit in any non-neighbor straw
          if (has_non_neighbor_straw(event,strawnum)){
            has_shower[strawnum].push_back(true);
          }else{
            has_shower[strawnum].push_back(false);
          }

          // check if hit in this straw
          long strawTimeCountCal, strawTimeCountHV;
          double totHV,totCal;
          double maxval;
          if (!valid_straw_data(event,strawnum,strawTimeCountCal,strawTimeCountHV,maxval) || !valid_tot_data(event,strawnum,totHV,totCal)){
            has_straw[strawnum].push_back(false);
            alltimes[strawnum].push_back(0);
            alldts[strawnum].push_back(0);
            allmaxvals[strawnum].push_back(0);
            alltots[strawnum].push_back(0);
            continue;
          }
          
          has_straw[strawnum].push_back(true);
          //long strawTimeCount = (strawTimeCountHV + strawTimeCountCal)/2.;
          long strawTimeCount = strawTimeCountCal;
//          double temp = (double) (strawTimeCountCal - pmtTopTimeCount) * 15.625*1e-3; // convert to ns
//          double temp2 = (double) (strawTimeCountHV - pmtTrigTimeCount) * 15.625*1e-3; // convert to ns
          double driftTimeTop = (double) (strawTimeCount - pmtTopTimeCount) * 15.625*1e-3; // convert to ns
          double driftTimeBot = (double) (strawTimeCount - pmtBotTimeCount) * 15.625*1e-3; // convert to ns
          double driftTimeTrig = (double) (strawTimeCount - pmtTrigTimeCount) * 15.625*1e-3; // convert to ns
          if (fabs(driftTimeTop) > 1000 || fabs(driftTimeBot) > 1000 || fabs(driftTimeTrig) > 1000)
            std::cout << filenames[i] << " " << k << " " << driftTimeTop << " " << driftTimeBot << " " << driftTimeTrig << " " << j << " " << last_event_with_straw << std::endl;
          //  std::cout << filenames[i] << " " << k << " " << temp << " " << lastPmtTopTimeCount - pmtTopTimeCount << " " << lastPmtTrigTimeCount - pmtTrigTimeCount << " " << pmtTopTimeCount-pmtTrigTimeCount << std::endl;

          //std::cout << filenames[i] << " " << k << " " << temp << " : " << lastPmtTopTimeCount << " " << pmtTopTimeCount << " : " << lastPmtTrigTimeCount << " " << pmtTrigTimeCount << " " << pmtTopTimeCount-pmtTrigTimeCount << std::endl;
          if (use_pmt_time == USE_TRIG_TIME){
            alltimes[strawnum].push_back(driftTimeTrig);
          }else if (use_pmt_time == USE_TOP_TIME){
            alltimes[strawnum].push_back(driftTimeTop);
          }else{
            alltimes[strawnum].push_back((driftTimeTop+driftTimeBot)/2.);
          }
          alldts[strawnum].push_back((strawTimeCountCal-strawTimeCountHV)*0.015625);
          alltots[strawnum].push_back((totHV+totCal)/2.);
          allmaxvals[strawnum].push_back(maxval);
          if (has_straw[strawnum][has_straw[strawnum].size()-1] && maxval < pmpminimum){
            has_shower[strawnum][has_shower[strawnum].size()-1] = true;
          }
        }
      } // loop over straws

      lastPmtTopTimeCount = pmtTopTimeCount;
      lastPmtTrigTimeCount = pmtTrigTimeCount;
      
    } // loop over events
  } // loop over files

  std::cout << " TOTAL RUN TIME: " << total_time/60./60. << " hours, " << alltop_x.size() << " events" << std::endl;
}
