#pragma once
#ifndef WELLBORE_PRESSURE_DROP_H
#define WELLBORE_PRESSURE_DROP_H
#include <string>
#include <list>
#include <math.h>

using namespace std;
double PI = 3.1415926;
/*
	Define the wellbore pressure drop model for lost circulation.
*/

	// Fluid Property

	 double K_1 = 0.57;        // consistensy of fluid
	 double m = 0.59;        // fluid flow index
	 double t_y_1 = 20.3;      // yield stress of the fluid
	 double e_den_1 = 9;       // fluid density
	 double g = 9.81;

	// Geometry of the drill string

	double D_po_1 = 4.5;		 // the outer diameter of drill pipe
	double D_pi_1 = 3.826;		 // the inner diameter of drill pipe
	double D_hl_1 = 8.5 ;		 // the diameter of wellbore
	double H1_1 = 2000; 
	float stdoff = 0;
	int H2_1 = 5000; 
	int H3_1 = 5000;						// the length of the well section
	int theta = 50;						// wellbore inclination for the H2 section

	// Boundary and Initial condition

	const double Q_1 = 300;				// the injection flow rate

	// Transform the Field Unit to SI

	double Q = Q_1 * (3.785 * pow(10, -3) / 60);

	// tansform inch to m
	double D_po = D_po_1 * 25.4 / 1000;
	double D_pi = D_pi_1 * 25.4 / 1000;
	double D_hl = D_hl_1 * 25.4 / 1000;

	// transform ft to m
	double H1 = H1_1 * 0.3048;
	double H2 = H2_1 * 0.3048;
	double H3 = H3_1 * 0.3048;
	double MD = H1 + H2 + H3;

	// transform lbm/gal to kg/cm^3
	double e_den = e_den_1 * 119.83;

	// pressure transform
	double t_y = t_y_1 * 4.44822/ 100 / 0.092903;

	// flow behavior index
	double K = K_1 * 4.44822 / 100 / 0.092903;
	double w = PI*(D_hl + D_po) / 2;      // the width of the annulus
	double h = (D_hl - D_po) / 2;	      
	double A_anu = w * h;
	double A_pip = PI*pow(D_pi, 2) / 4;

#endif // !WELLBORE_PRESSURE_DROP_H