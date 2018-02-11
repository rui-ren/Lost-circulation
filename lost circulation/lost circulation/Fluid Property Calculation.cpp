#include "wellbore_pressure_drop.h"
#include <math.h>
#include <iostream>
using namespace std;

double Annulus_Pressure_Drop( )
{
	double v_av = Q / A_anu;							// average velocity
	double r_s = ((1 + 2 * m) / 3 / m)*(12 * v_av) / (D_hl - D_po);     // shear rate
	double t_w1 = t_y + K * pow(r_s, m);				// initial t_w1
	double e = 1;

	do
	{
		double x = t_y / t_w1;
		double Ca = 1 - x / (1 + m) - m* pow(x, 2) / (1 + m);
		double t_w2 = t_y + K * pow((r_s / Ca), m);
		e = (t_w2 - t_w1) / t_w1;
		t_w1 = t_w2;
	} while (e < 0.01);

	// Calculate Generalized flow-behavior index. (Ahmed and Miska 2009)

	cout << "wall shear stress is: " << t_w1 << endl;
	double x = t_y / t_w1;

	double P = 3 * m / (1 + 2 * m)*(1 - 1 / (1 + m)* x - (m / (1 + m))*pow(x, 2));
	double N = 1 / (3 / P - 2);					// Generalize fluid flow index
	cout << "the generalized fluid flow index(0.15 < N < 0.4): " << N << endl;

	// Use the correlation to calculate the critical Reynolds number for laminar and turbulent flow

	double k = D_pi / D_pi;			   // Geometry index of the pipe and wellbore
	// Reynolds number for laminar flow
	double Re1 = 2100 * (pow(N, 0.33)*(1 + 1.402*k - 0.977*pow(k, 2)) - 0.019*stdoff*pow(N, -0.868)*k); 
	// Reynolds number for turbulent flow
	double Re2 = 2900 * (pow(N, (-0.039*pow(Re1, 0.307))));
	// Reynolds number for Yield Power Law flow
	double Re_YPL = 12 * e_den*v_av / t_w1;

	double f_l = 24 / Re_YPL;
	// calculate the  flow friction factor
	if (Re1 > Re_YPL)      // for laminar flow
	{
		if (stdoff != 0)
		{
			double R = (1 - 0.072*stdoff / N*pow(k, 0.8454) - 1.5*pow(stdoff, 2)*sqrt(N)*pow(k, 0.1852) + 0.96*pow(stdoff, 2)*sqrt(N)*pow(k, 0.2527));
			f_l = f_l * R;
			double f = f_l;
			double P_drop = 4 * f * t_w1 / (D_po - D_pi);
			cout << "the pressure drop for annulus: " << P_drop << endl;
		}
	}
	else
	{
		if (Re_YPL > Re2)	// calculate the turbulent flow friction factor
		{
			// ??? how to calculate the friction factor.
			double f_t;
			double f = f_t;

			if (Re_YPL > Re1 & Re_YPL < Re2)         // calculate the transitional fluid flow friction factor
			{
				double f_tran = f_l + (Re_YPL - Re1)*(f_t - f_l) / (Re2 - Re1);
				double f = f_tran;
			}
			double P_drop = 4 * f * t_w1 / (D_po - D_pi);
			cout << "the pressure drop for annulus: " << P_drop << endl;
		}
	}
}