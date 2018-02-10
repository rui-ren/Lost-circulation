#pragma once
#include <string>
#include <list>
using namespace std;
/*
	Define the wellbore pressure drop model for lost circulation.
*/

// Fluid Property

class Fluid_Property
{
public:
	const double K = 0.57;       // consistensy of fluid
	const double m = 0.59;       // fluid flow index
	const double t_y = 20.3;     // yield stress of the fluid
};

// Geometry of the drill string

class Well_Geometry
{
public:
	const double D_po = 4.5;		 // the outer diameter of drill pipe
	const double D_pi = 3.826;		 // the inner diameter of drill pipe
	const double D_hl = 8.5 ;		 // the diameter of wellbore
	const double H1 = 2000, H2 = 5000, H3 = 5000;   // the length of the well section
	double theta = 50;		 // wellbore inclination for the H2 section
};

// Boundary and Initial condition

class Boundary
{
protected:
	const double Q = 9.00;		// the injection flow rate
};
