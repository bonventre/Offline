#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>

#include "fit/StrawDrift.hh"

StrawDrift::StrawDrift(std::string filename, float wirevoltage, int _bins)
{
  this->bins = _bins;
  std::ifstream myfile(filename.c_str());

  std::vector<double> dataEField;
  std::vector<double> dataVInst;
  std::vector<double> dataDistances; //the set of distances that correspond to a given E-field

  //read the file
  std::string line;
  double a, b;
  while (getline(myfile, line))
  {
    //handle file formating
    std::istringstream iss(line);
    if ((iss >> a >> b)){
      a = a*100.0; //convert the E-field from KV/cm to V/mm
      b = b*0.01;//convert the speed from cm/us to mm/ns (10/1000 = 0.01)
      //fill the vector of structs
      dataEField.push_back(a);
      dataVInst.push_back(b);
    }
    else {
      //cout << "skipping line"<<"\n";
    }
  }

  //Use the E:insta-velc tables to build d:insta-veloc tables based on the voltage input 

  //define the wire and straw radius in mm (remove the hard coding?) 
  float wireradius = 12.5/1000.; //12.5 um in mm 
  float strawradius = 2.5; //2.5 mm in mm 

  // calculate the distances that correspond to the efields listed in the table (fix units!!) 
  for (size_t i=0; i < dataEField.size(); i++) { 
    dataDistances.push_back(wirevoltage/((dataEField[i])*log(strawradius/wireradius))); //in mm 
  } 

  // interpolate to get 50x more points for distances and instantSpeeds 
  for (size_t i=0;i < bins;i++){
    double this_distance = 0.001 * i;
    int j=1;
    while (j < dataDistances.size()-1){
      if (dataDistances[j] < this_distance){
        break;
      }
      j++;
    }
    double v_high = dataVInst[j-1];
    double v_low = dataVInst[j];
    double d_high = dataDistances[j-1];
    double d_low = dataDistances[j];
    this->distances.push_back(this_distance);
    this->instantSpeeds.push_back(v_low + (v_high-v_low)/(d_high-d_low) * (this_distance-d_low));
  }

  //numerically integrate distances and instantSpeeds
  double total_time = 0;
  this->averageSpeeds.push_back(this->instantSpeeds[0]);
  this->times.push_back(0);
  for (size_t i=1;i<this->distances.size();i++){
    double slice_distance = this->distances[i]-this->distances[i-1];
    double total_distance = this->distances[i];
    double this_speed = this->instantSpeeds[i];
    double slice_time = slice_distance/this_speed;
    total_time += slice_time;
    this->averageSpeeds.push_back(total_distance/total_time);
    this->times.push_back(total_time);
//    if (i%50 == 0)
//      std::cout << i << " " << this->distances[i] << " " << this->instantSpeeds[i] << " " << this->averageSpeeds[i] << " " << this->times[i] << std::endl;
  }
}

double StrawDrift::D2T(double distance)
{
  size_t bin = distance/0.001;
  if (bin >= bins)
    bin = bins;
  return this->times[bin];
  return 0;
}

double StrawDrift::T2D(double time)
{
  return 0;
}
