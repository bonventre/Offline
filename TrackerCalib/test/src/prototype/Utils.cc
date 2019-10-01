#include <iostream>
#include <TFile.h>
#include <TTree.h>

#include "prototype/Utils.hh"

double calculate_longitudinal_distance(double top_y, double bot_y, double longitudinal, double vertical, double separation)
{
    // longitudinal position in pixel counts at vertical position of wire
    double pixelrow = bot_y + (top_y-bot_y)*(vertical/separation);

    // longitudinal distance from track to wire in mm
    return pixelrow*0.250 - longitudinal;
}

double calculate_DOCA(double top, double bot, double horizontal, double vertical, double separation)
{
    // angle from vertical (hypotenus is track, opposite is horizontal distance, adjacent is vertical distance
    double angle = atan(((top-bot)*0.05)/separation);

    // horizontal position in pixel counts at vertical position of wire
    double pixelrow = bot + (top-bot)*(vertical/separation);

    // horizontal distance from track to wire in mm
    double horizontal_dist = pixelrow*0.05 - horizontal;

    // angle between vector of closest approach and horizontal
    // so cos angle x horizontal dist gives dist of closest approach
    double dist = fabs(horizontal_dist * cos(angle));

    return dist;
}

void calculate_relative_position(double top, double bot, double horizontal, double vertical, double separation, double &x, double &y)
{
    // angle from vertical (hypotenus is track, opposite is horizontal distance, adjacent is vertical distance
    double angle = atan(((top-bot)*0.05)/separation);
    //double angle = atan(PIXEL_VERTICAL_SEPARATION/((top-bot)*0.05+1e-5));

    // horizontal position in pixel counts at vertical position of wire
    double pixelrow = bot + (top-bot)*(vertical/separation);

    // horizontal distance from track to wire in mm
    //double horizontal_dist = 16.8 - pixelrow*0.05 - horizontal;
    double horizontal_dist = pixelrow*0.05 - horizontal;

    // angle between vector of closest approach and horizontal
    // so cos angle x horizontal dist gives dist of closest approach
    double dist = horizontal_dist * cos(angle);
    //double dist = fabs(horizontal_dist * sin(angle));
    
    x = dist * cos(angle);
    y = -1 * dist * sin(angle);
}

double calculate_max_drift(double top, double bot, double horizontal, double vertical, double separation, double circle_x, double circle_y)
{
  if (calculate_DOCA(top,bot,circle_x,circle_y,separation) > 2.5)
    return -1;
  
  double top_x = top*0.05 - circle_x;
  double bot_x = bot*0.05 - circle_x;
  double top_y = separation - circle_y;
  double bot_y = - circle_y;

  double dx = bot_x-top_x;
  double dy = bot_y-top_y;
  double dr = sqrt(dx*dx+dy*dy);
  double D = top_x*bot_y-bot_x*top_y;

  double x1,x2,y1,y2;
  if (dy > 0){
    x1 = (D*dy + dx*sqrt(2.5*2.5*dr*dr-D*D))/(dr*dr);
    x2 = (D*dy - dx*sqrt(2.5*2.5*dr*dr-D*D))/(dr*dr);
  }else{
    x1 = (D*dy - dx*sqrt(2.5*2.5*dr*dr-D*D))/(dr*dr);
    x2 = (D*dy + dx*sqrt(2.5*2.5*dr*dr-D*D))/(dr*dr);
  }
  y1 = (-D*dx + fabs(dy)*sqrt(2.5*2.5*dr*dr-D*D))/(dr*dr);
  y2 = (-D*dx - fabs(dy)*sqrt(2.5*2.5*dr*dr-D*D))/(dr*dr);

  double dist1 = sqrt(pow(x1-(horizontal-circle_x),2)+pow(y1-(vertical-circle_y),2));
  double dist2 = sqrt(pow(x2-(horizontal-circle_x),2)+pow(y2-(vertical-circle_y),2));

  if (dist1 > dist2)
    return dist1;
  
  return dist2;
}

  // x[0] is distance from back of pixel to wire horizontally perpendicular to the straw in mm
  // x[1] is distance from bot pixel to wire in mm
void draw_vplot(TH2F *vplot, std::vector<double> &top, std::vector<double> &bot, std::vector<double> &times, double timeoffset, double hdist, double vdist, double pixelseparation)
{
  for (int j=0;j<top.size();j++){
    double dist = calculate_DOCA(top[j],bot[j],hdist,vdist,pixelseparation);
    vplot->Fill(times[j]-timeoffset,dist);
  }
}
