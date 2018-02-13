#pragma once
#include "wellbore_pressure_drop.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
using namespace std;

// define the wellbore model and associate with the formation model!

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

	// Calculate Generalized flow-behavior index. (Ahmed and Miska 2009)

	double x = t_y / t_w1;
	double P = 3 * m / (1 + 2 * m)*(1 - x / (1 + m) - (m / (1 + m))*pow(x, 2));
	double N = m * P / (1 + 2 * m * (1 - P));
	cout << "the generalized fluid flow index(0.15 < N < 0.4): " << N << endl;

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
	
	// notice the pressure drop should include the pressure drop due to gravity
	//for (int i = 1; i <= (MD/10); i++)
	//{
	//	int md ;
	//	md  = i * 10;
	//	if (md < H3)
	//		cout <<"Annulus pressure drop for H1: "<< (md * P_drop ) / 1000000 << endl;
	//	if (md < (H3+H2) && md > (H3))
	//	{
	//		cout << "Annulus pressure drop for H2:" << (md * P_drop + e_den * g * (md-H3) *cos(theta*PI / 180)) / 1000000 << endl;
	//	}
	//	if (md > (H3+H2)&& md <(H1+H2+H3)) // Horizontal well 
	//		cout << "Annulus pressure drop for H3: " << (md * P_drop + e_den * g * H2 * cos(theta*PI / 180)+ e_den*g*(md-H2-H3)) / 1000000 << endl;
	//}
	//system("pause");

	//////////////////////// Calculate the pressure drop in pipe //////////////////////////

	//Step 1: Calculate the wall shear stress
		double v_pipe = Q / A_pip;									  // average velocity in the pipe
		double r_s_p = ((1 + 3 * m) / 4 / m)*(8 * v_pipe) / (D_pi);     // shear rate for the pipe
		double t_w1_p = t_y + K * pow(r_s_p, m);		                  // initial t_w1
		double e_p = 1;
		do
		{
			double x = t_y / t_w1_p;
			double Ca = 1 - x / (2*m+1);
			double t_w2 = t_y + K * pow((r_s_p / Ca), m);
			e = (t_w2 - t_w1) / t_w1_p;
			t_w1_p = t_w2;
		} while (e_p < 0.01);
	// Step 2: Calculate the Reynolds number
		double y = (log(m)+3.93)/50;
		double z = (1.75 - log(m)) / 7;
		double Cc = 1 - (1 / (2 * m + 1))*(t_y/t_w1);
		// Critical Reynolds Number
		double Re_cr = pow((4 * (3 * m + 1) / m / y), (1 / (1 - z)));
		double Re_eq = ((6 * m + 2) / m)*(e_den*pow(v_pipe, (2 - m))*(D_pi/2))/((t_y*pow((D_pi/2/v_pipe),m))+K*((3*m+1)/m/Cc));

		double f = 16 / Re_eq*(3 * m + 1) / 4 / m;

		if (Re_cr > Re_eq)  // Laminar flow
		{
			double f1 = y * pow(Cc*Re_eq, (-z));
			f = f1;
		}
		double pressure_drop_pipe = 2 * f * e_den*pow(v_pipe, 2) / D_pi;

		// Step 3: Calculate the pressure drop cross the bit
		float Db = 0.5;
		double At = 3 * PI * pow(Db, 2) / 4;
		double P_bit = e_den_1 * pow(Q_1,2)/(pow(At,2)*10858)/145;
		cout << "bit pressure drop: "<< P_bit << endl;

		double P_f_H3 = 0;
		double P_f_H2 = 0;
		double P_f_H1 = 0;

		// Step 4: Calculate the pressure drop along the pipe and get the standpipe pressure.
		for (int i = 1; i <= (MD / 10); i++)
		{
			double md;
			double section;
			section = i * 10;
			md = MD - section;    // the depth to wellhead.
			if (md > (H2 + H1) && md < (H1 + H2 + H3))
			{
				P_f_H3 = (section * pressure_drop_pipe) / 1000000;
				cout << "Pipe pressure drop for H3: " << P_f_H3 + P_bit << endl;
			}
			if (md > H1 && md < (H2 + H1))
			{
				P_f_H2 = (section * pressure_drop_pipe )/1000000;  // friction pressure drop
				cout << "Pipe pressure drop for H2:" << P_f_H2 + e_den * g * (section - H3)*cos(theta*PI / 180) / 1000000 + P_bit << endl;
			}
			if (md < H1 && md > 0) // Horizontal well 
			{
				// H3 section pressure drop acoording to the gravity and friction, plus H2 section
				P_f_H1 = (section * pressure_drop_pipe) / 1000000;
				cout << "Pipe pressure drop for H1: " << P_f_H1 + e_den * g * H2*cos(theta*PI / 180) / 1000000 + e_den*g* (section - H3 - H2) / 1000000 + P_bit << endl;
			}
		}
		system("pause");
}
