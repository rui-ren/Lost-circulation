#pragma once
#include "wellbore_pressure_drop.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
using namespace std;

int main() {

	///////////////////////////////calculate the pressure drop////////////////////////////

	// calculate the wall shear stress
	double v_av = Q / A_anu;							// average velocity
	double r_s = ((1 + 2 * m) / 3 / m)*(12 * v_av) / (D_hl - D_po);     // shear rate
	double t_w1 = t_y + K * pow(r_s, m);				// initial t_w1
	double e = 1;
	// calculate the wall shear stress
	do
	{
		double x = t_y / t_w1;
		double Ca = 1 - x / (1 + m) - m* pow(x, 2) / (1 + m);
		double t_w2 = t_y + K * pow((r_s / Ca), m);
		e = t_w2 - t_w1;
		t_w1 = t_w2;
	} while (e > 0.00000001);
	cout << "wall shear stress is: " << t_w1 << endl;
	system("pause");

	// Calculate Generalized flow-behavior index. (Ahmed and Miska 2009)

	double x = t_y / t_w1;
	double P = 3 * m / (1 + 2 * m)*(1 - x / (1 + m) - (m / (1 + m))*pow(x, 2));
	double N = m * P / (1 + 2 * m * (1 - P));
	cout << "the generalized fluid flow index(0.15 < N < 0.4): " << N << endl;
	system("pause");

	// Use the correlation to calculate the critical Reynolds number for laminar and turbulent flow

	// Geometry index of the pipe and wellbore
	double k = D_pi / D_po;
	// Critical Reynolds number for laminar flow
	double Re1 = 2100 * (pow(N, 0.331)*(1 + 1.402*k - 0.977*pow(k, 2)) - 0.019*stdoff*pow(N, -0.868)*k);
	// Critical Reynolds number for turbulent flow
	double Re2 = 2900 * (pow(N, (-0.039*pow(Re1, 0.307))));
	// Reynolds number for Yield Power Law flow
	double Re_YPL = 12 * e_den * v_av / t_w1;
	// Laminar flow friction factor
	double f_Lam = 24 / Re_YPL;

	cout << f_Lam << endl;
	cout << Re2 << endl;
	cout << Re_YPL << endl;
	system("pause");
	double f_1 = f_Lam;
	// calculate the  flow friction factor
	if (Re1 > Re_YPL)        // for laminar flow
	{
		if (stdoff != 0)
		{
			double R = (1 - 0.072*stdoff / N*pow(k, 0.8454) - 1.5*pow(stdoff, 2)*sqrt(N)*pow(k, 0.1852) + 0.96*pow(stdoff, 2)*sqrt(N)*pow(k, 0.2527));
			f_1 = f_Lam * R;
		}
	}
	else
	{
		if (Re_YPL > Re2)	// calculate the turbulent flow friction factor
		{
			// ??? how to calculate the friction factor.
			double f_t = 0;
			double f = f_t;

			if (Re_YPL > Re1 && Re_YPL < Re2)         // calculate the transitional fluid flow friction factor
			{
				double f_tran = f_Lam + (Re_YPL - Re1)*(f_t - f_Lam) / (Re2 - Re1);
				double f_1 = f_tran;
			}
		}
	}
	double P_drop = 4 * f_1 * t_w1 / (D_po - D_pi);
	cout << "the pressure drop for annulus: " << P_drop << endl;
	system("pause");
	
	// notice the pressure drop should include the pressure drop due to gravity
	for (int i = 1; i <= (MD/50); i++)
	{
		int md ;
		md  = i * 50;
		if (md < H1)
			cout <<"Annulus pressure drop for H1: "<< (md * P_drop + e_den * g * md) / 1000000 << endl;
		if (md < H2 && md > H1)
		{
			cout << "Annulus pressure drop for H2:" << (md * P_drop*cos(theta*PI / 180) + e_den * g * md) / 1000000 << endl;
		}
		if (md > H2) // Horizontal well 
			cout << "Annulus pressure drop for H3: " << (md * P_drop + e_den * g * H2 + H2 * P_drop*cos(theta*PI / 180)) / 1000000 << endl;
	}
	system("pause");

	//////////////////////// Calculate the pressure drop in pipe //////////////////////////



}
