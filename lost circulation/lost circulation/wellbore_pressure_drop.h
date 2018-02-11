#pragma once
#include <string>
#include <list>
#include <math.h>

using namespace std;
double PI = 3.1415926;
/*
	Define the wellbore pressure drop model for lost circulation.
*/

// Fluid Property

	const double K = 0.57;        // consistensy of fluid
	const double m = 0.59;        // fluid flow index
	const double t_y = 20.3;      // yield stress of the fluid
	const double e_den = 9;       // fluid density

// Geometry of the drill string

	const double D_po = 4.5;		 // the outer diameter of drill pipe
	const double D_pi = 3.826;		 // the inner diameter of drill pipe
	const double D_hl = 8.5 ;		 // the diameter of wellbore
	const double H1 = 2000; 
	float stdoff = 0;
	int H2 = 5000; 
	int H3 = 5000;						// the length of the well section
	int theta = 50;						// wellbore inclination for the H2 section
	double w = PI*(D_hl - D_po)/2;      // the width of the annulus
	double h =  (D_hl - D_po) / 2;	    // the 
	double A_anu = w * h;
	double A_pip = PI*pow(D_pi, 2) / 4;

// Boundary and Initial condition

	const double Q = 300;				// the injection flow rate