#include "wellbore_pressure_drop.h"
#include <math.h>
#include <iostream>
using namespace std;

double Pressure_Drop( )
{

	double v_av = Q / A;							// average velocity
	double r_s = ((1 + 2 * m) / 3 / m)*(12 * v_av) / (D_po - D_pi);
	double t_w1 = t_y + K * pow(r_s, m);		// initial t_w1
	double e = 1;

	do
	{
		double x = t_y / t_w1;
		double Ca = 1 - x / (1 + m) - m* pow(x, 2) / (1 + m);
		double t_w2 = t_y + K * pow((r_s / Ca), m);
		e = (t_w2 - t_w1) / t_w1;
		t_w1 = t_w2;
	} while (e < 0.01);

	cout << "wall shear stress is: " << t_w1 << endl;
	double P_drop = 4 * t_w1 / (D_po - D_pi);
	cout << "the pressure drop for annulus: " << P_drop << endl;

}





