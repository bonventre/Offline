#include <vector>
#include <string>

class StrawDrift {
  public:
    StrawDrift(std::string filename, float wirevoltage, int _bins);
    ~StrawDrift(){};

    double D2T(double dist);
    double T2D(double time);

  private:

    int bins;
    std::vector<double> distances;
    std::vector<double> instantSpeeds;
    std::vector<double> averageSpeeds;
    std::vector<double> times;
};
