#include <vector>
#include <string>

#include "event/Hit.h"

struct ParsedData{
  ParsedData(){
    has_straw = std::vector<std::vector<bool> >(8, std::vector<bool>());
    has_shower = std::vector<std::vector<bool> >(8, std::vector<bool>());

    alltimes = std::vector<std::vector<double> > (8, std::vector<double>());
    alldts = std::vector<std::vector<double> > (8, std::vector<double>());
    allmaxvals = std::vector<std::vector<double> > (8, std::vector<double>());
    alltots = std::vector<std::vector<double> > (8, std::vector<double>());
  }
  std::vector<double> alltop_x; // all hits top pixel avg perpendicular to straw [pixel number (50um)]
  std::vector<double> allbot_x; // all hits bot pixel avg perpendicular to straw [pixel number (50um)]
  std::vector<double> alltop_y; // all hits top pixel avg along straw [pixel number (250um)]
  std::vector<double> allbot_y; // all hits bot pixel avg along straw [pixel number (250um)]

  std::vector<std::vector<bool> > has_straw; // has a hit in straw i
  std::vector<std::vector<bool> > has_shower; // has a hit in a straw not neighbors with straw i

  std::vector<std::vector<double> > alltimes; // drift time measurement in straw i
  std::vector<std::vector<double> > alldts; // NOT USED
  std::vector<std::vector<double> > allmaxvals; // NOT USED
  std::vector<std::vector<double> > alltots; // NOT USED
};


#define USE_TRIG_TIME 0
#define USE_TOP_TIME 1
#define USE_TOPBOT_TIME 2

// checks if non-neighboring straw to this channel are hit
bool has_non_neighbor_straw(Event *event, int channel);

// checks if non-neighboring straws are hit
bool has_straw_shower(Event *event);

// checks if pmt times are valid
bool valid_pmt_data(Event *event, long &pmtTopTime, long &pmtBotTime, long &pmtTrigTime, long &lastPmtTopTime, long &lastPmtTrigTime);

// checks if both pixels hit
// and hits in each pixel are adjacent
bool valid_pixel_data(Event *event, double &topavg_x, double &topavg_y, double &botavg_x, double &botavg_y, bool old_data);

// checks if hit in this straw
bool valid_straw_data(Event *event, int channel, long &strawtimecal, long &strawtimehv, double &maxval);
bool valid_tot_data(Event *event, int channel, double &tothv, double &totcal);

void parse_data_files(std::vector<std::string> &filenames, int strawnum, std::vector<double> &alltop_x, std::vector<double> &allbot_x, std::vector<double> &alltop_y, std::vector<double> &allbot_y, std::vector<std::vector<bool> > &has_straw, std::vector<std::vector<bool> > &has_shower, std::vector<std::vector<double> > &alltimes, std::vector<std::vector<double> > &alldts, std::vector<std::vector<double> > &allmaxvals, std::vector<std::vector<double> > &alltots, bool old_data=false, double pmpminimum = -9e9, int use_pmt_time=USE_TOP_TIME);
