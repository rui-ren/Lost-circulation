#include "wellbore_pressure_drop.h"
#include <math.h>
#include <iostream>
using namespace std;

double v_av = Q/A;		     // average velocity

int r_s = ((1 + 2 * m) / 3 / m)*(12 * v_av) / (D_po - D_pi);
double t_w = t_y + pow(r_s, m);
double x = t_y / t_w;
double Ca = 1 - x / (1 + m) - m* pow(x, 2) / (1 + m);
double d = (1 + 2 * m) / (Ca * r_s);
double t_w1 = t_y + K * pow(d,m);

