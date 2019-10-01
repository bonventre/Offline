#include <vector>
#include <string>

#include <TH2F.h>

void calculate_relative_position(double top, double bot, double horizontal, double vertical, double separation, double &x, double &y);
double calculate_max_drift(double top, double bot, double horizontal, double vertical, double separation, double circle_x, double circle_y);

double calculate_longitudinal_distance(double top_y, double bot_y, double longitudinal, double vertical, double separation = 20.215);
double calculate_DOCA(double top, double bot, double horizontal, double vertical, double separation = 20.215);
void draw_vplot(TH2F *vplot, std::vector<double> &top, std::vector<double> &bot, std::vector<double> &times, double timeoffset, double hdist, double vdist, double pixelseparation);

