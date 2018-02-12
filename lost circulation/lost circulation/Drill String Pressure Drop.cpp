#include <math.h>
#include <iostream>
using namespace std;

// Calculate the wall shear stress in the pipe

/*

double Pipe_Pressure_Drop()
{
double v_av = Q / A_pip;									  // average velocity in the pipe
double r_s = ((1 + 3 * m) / 4 / m)*(8 * v_av) / (D_pi);     // shear rate for the pipe
double t_w1 = t_y + K * pow(r_s, m);		                // initial t_w1

double e = 1;
do
{
double x = t_y / t_w1;
double Ca = 1 - x / (1 + m) - m* pow(x, 2) / (1 + m);
double t_w2 = t_y + K * pow((r_s / Ca), m);
e = (t_w2 - t_w1) / t_w1;
t_w1 = t_w2;
} while (e < 0.01);


return 0;
}

*/
