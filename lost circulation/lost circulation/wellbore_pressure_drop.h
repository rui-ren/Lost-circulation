#pragma once
#ifndef WELLBORE_PRESSURE_DROP_H
#define WELLBORE_PRESSURE_DROP_H
#include <string>
#include <list>
#include <math.h>

using namespace std;
double PI = 3.1415926;

	//Define wellbore pressure drop model

	class wellbore
	{
	protected:

	public:
		wellbore();      // constructor
		~wellbore();		 // destructor
		double wellbore::Annulus_pressure_drop(double theta, double g, double D_po, double stdoff, double m, double D_pi, double D_hl, double H1, double H2, double H3, double MD, double e_den, double t_y, double K, double w, double h, double A_anu, double A_pip, double Q);
		double wellbore::Drillpipe_pressure_drop(double D_po, double m, double D_pi, double H1, double H2, double H3, double MD, double e_den, double t_y, double K, double w, double h, double A_anu, double A_pip, double Q);

	};
 // bit pressure drop

	//const double Q_1 = 300;				// the injection flow rate
	// Transform the Field Unit to SI
	//double Q = Q_1 * (3.785 * pow(10, -3) / 60);
	// tansform inch to m
	// transform ft to m
	//double H1 = H1_1 * 0.3048;
	//double H2 = H2_1 * 0.3048;
	//double H3 = H3_1 * 0.3048;
	//double MD = H1 + H2 + H3;
	// transform lbm/gal to kg/cm^3
	//double e_den = e_den_1 * 119.83;
	//double t_y = t_y_1 * 4.44822/ 100 / 0.092903;
	//double K = K_1 * 4.44822 / 100 / 0.092903;
	//double w = PI*(D_hl + D_po) / 2;      
	//double h = (D_hl - D_po) / 2;	      
	//double A_anu = w * h;
	//double A_pip = PI*pow(D_pi, 2) / 4;
#endif // !WELLBORE_PRESSURE_DROP_H
