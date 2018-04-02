#include <math.h>
#include <iostream>

#pragma once
#include"wellbore1.h"
#include<stdlib.h>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<fstream>
#include<string>
#include<vector>
#include<stdio.h>
#include<math.h>

using namespace std;

int rowA = 0;
int colA = 0;

int main(){
	////////////////////////////////////// read the data file of the well trajectoy/////////////////////////////////

	string lineA;
	float x;
	string filename;
	float arrayA[10000][3] = { { 0 } };   // define the array of the data.
	ifstream fileIN;

	// Intro
	cout << "input the file of well log data: measure depth, inclination angle and azimuth" << endl;
	fileIN.open("1.txt");

	// Error Check
	if (fileIN.fail())
	{
		cout << "this file cannot be open or access" << endl;
	}

	// reading the well log data
	cout << "\n" << endl;
	while (fileIN.good()){
		while (getline(fileIN, lineA)){
			istringstream streamA(lineA);
			colA = 0;
			while (streamA >> x){
				arrayA[rowA][colA] = x;
				colA++;
			}
			rowA++;
			}
		}

	// display data
	cout << "# of row ------>" << rowA << endl;
	cout << "# of colums----->" << colA << endl;
	cout << " " << endl;
	for (int i = 0; i < rowA; i++){
		for (int j = 0; j < colA; j++)
		{
			cout << left << setw(6) << arrayA[i][j] << " ";
		}
		cout << endl;
	}
	fileIN.close();
	
	//std::string line_;
	//ifstream file_("1.txt");
	//int depth;
	//int inc;
	//int azi;
	//if (file_.is_open())
	//{
	//	while (file_>> depth >> inc >> azi);
	//	{
	//		std::cout << " the varaible in file " << depth << " " << inc << " " << azi << endl;
	//	}
	//	file_.close();
	//}
	//else
	//	std::cout << "file is not open" << "\n";

	/////////////////////////////// input the parameter of the wellbore and formation data/////////////////////////////
	// caution!! this is field unit.

	double D_po = 5.875;         // pipe outer diameter
	D_po = D_po * 25.4 / 1000;
	double D_pi = 5.0 - 0.3 * 2;           // pipe inner diameter
	double h_s1_length = 6200;   // first section of the hole
	double h_s2_length = 8365;   // second section of the hole
	double h_s3_length = 1287;   // third section of the hole

	double D_h_s1 = 9.56;    // first section inner diameter of the hole
	double D_h_s2 = 9.76;   // second section inner diameter of the hole
	double D_h_s3 = 9.5;    // third section inner diameter of the hole

	double e_den = 12.6;   // ppg
	double t_y = 10.06;    // yield shear stress
	double K = 0.3285;      // flow consistency
	double m = 0.7323;      // flow index
	double g = 9.8;

	// cutting transport parameter
	double ROP = 20;                     // Check the unit ft/hour 
	double ds = 2/1000000;                      // Particle diameter of the      cm^2
	double h = 1/100;                       // The thickness of the cutting   1cm
	double e_cutting = 2.8;             // The density of the cutting
	double cutting_concentration = 1;   // cuting concentration calculation model .   Feifei zhang or API 13 RD
	// Change the unit to SI unit

	// tool joint information
	double D_to = 6.2;                   // assume 6 inch
	double L_tooljoint = 24;            // assume 24 cm

	 D_to = D_to * 25.4 / 1000;
	 L_tooljoint = L_tooljoint / 100;    // change to meter

	// ????????????????????? flow rate need to calculate!

	double Q = 12;     // gpm
    	Q = Q * (3.785 * pow(10, -3) / 60);

	 //tansform inch to m
	 //transform ft to m

	  D_h_s1 = D_h_s1 * 25.4 / 1000;
	  D_h_s2 = D_h_s2 * 25.4 / 1000;
	  D_h_s3 = D_h_s3 * 25.4 / 1000;

	  h_s1_length = h_s1_length * 0.3048;   // first section of the hole
	  h_s2_length = h_s2_length * 0.3048;   // second section of the hole
	  h_s3_length = h_s3_length * 0.3048;   // third section of the hole

	 double s_length = h_s1_length + h_s2_length + h_s3_length;
	 //transform lbm/gal to kg/cm^3

	 e_den = e_den * 119.83;        // change ppg to kg/cm^3 
	 t_y = t_y * 4.44822/ 100 / 0.092903;
	 K = K * 4.44822 / 100 / 0.092903;

	double A_pip = PI*pow(D_pi, 2) / 4;     // drill string inner pipe diameter!!!

	double ToolJoint = 1;            // take tooljoint effect into consideration..
	double Tortuosity = 0;          // take the tortuosity into consideration..
	double Roughness_wellbore = 0;  // take roughness of wellbore into consideration..

	//////////////////////////////////////// frictional pressure drop /////////////////////////////////////////////////

	// define an vector for storing the frictional pressure data
	vector<int> frictional_pressure_drop;

	// inner annulus pressure drop for the well

	for (int i = 0; i < s_length; i++)
	{
		if (i < h_s1_length)
		{
				double stdoff = 0;    // vertical well, we assume the drill pipe in concentric. 
				double D_h = D_h_s1;  // the inner diameter of the casing
				 wellboreAnnulus_pressure_drop(D_po, D_to, L_tooljoint, stdoff, m, D_h, e_den, t_y, K, Q, ToolJoint, Roughness_wellbore, Tortuosity);   // pressure drop function.
			}
		if (i> h_s1_length && i < h_s2_length)
		{

				double stdoff = 0.70;    // input the stdoff of the drill string
				double D_h = D_h_s2;     // the inner diameter of the casing
				wellboreAnnulus_pressure_drop(D_po, D_to, L_tooljoint, stdoff, m, D_h, e_den, t_y, K, Q, ToolJoint, Roughness_wellbore, Tortuosity);   // pressure drop function.
		}
		else
		{
			double stdoff = 0.70;    // input the stdoff of the drill string
			double D_h = D_h_s3;
			wellboreAnnulus_pressure_drop(D_po, D_to, L_tooljoint, stdoff, m, D_h, e_den, t_y, K, Q, ToolJoint, Roughness_wellbore, Tortuosity);   // pressure drop function.
		}
	}

	//////////////////////////////////////// gravity pressure drop //////////////////////////////////////////////////

	// caution calculate the ROP and wellbore cutting concentration!! /////  

	// calculate the cutting concentration in the wellbore
	double C_cutting = 0;
	double e_mix = 0;

	if (cutting_concentration == 1)
	{
		// calculate the cutting concentration according to (Feifei Zhang. 2014)
		double A_w = PI * pow(D_h_s3, 2) / 4;
		double A_annulus = PI* (pow(D_h_s3, 2) - pow(D_po, 2)) / 4;
		// 1. superficial cutting concentration in annuli
		double v_sc = ROP * A_w / A_annulus;
		// 2. superficial liquid velocity in annuli
		double v_sl = Q / A_annulus;
		// 3. cutting feed concentration
		double CF = v_sc / (v_sc + v_sl);
		// 4. calculate the Reynold's number of the annuli (according to Ozbayoglu, 2002)
		double d_hyd = D_h_s3 - D_po;
		double d_eq = sqrt(pow(D_h_s3, 2) - pow(D_po, 2));
		double d_dim = d_hyd / d_eq;
		double a = 1.012 - 0.25* pow(d_dim, 1.563);
		double b = 0.51 - 0.242* pow(d_dim, 1.366);
		double K_prim = K * (a + b / m);
		double Re_general = pow(8, (1 - m))*e_den*pow(v_sl, (2 - m))*pow(d_hyd, m) / K_prim;
		// 5. calculate the particle slip velocity in a vertical well
		double v_slip = 0;
		if (Re_general > 1 && Re_general < 800)
		{
			double w = PI*(D_h_s3 + D_po) / 2;
			double h = (D_h_s3 - D_po) / 2;
			double A_anu = w * h;
			double v_av = Q / A_anu;
			double r_s = ((1 + 2 * m) / 3 / m)*(12 * v_av) / (D_h_s3 - D_po);     // shear rate      there is a problem!!!!  should be the generalized flow index in the formula.
			double t_w1 = t_y + K * pow(r_s, m);                                  // initial t_w1
			double e = 1000;
			do
			{
				double x = t_y / t_w1;
				double Ca = 1 - x / (1 + m) - m* pow(x, 2) / (1 + m);
				double t_w2 = t_y + K * pow((r_s / Ca), m);
				e = t_w2 - t_w1;
				t_w1 = t_w2;
			} while (e > 0.00000001);

			double x = t_y / t_w1;          // the eventually data from calculation.
			double Ca = (1 - x / (1 + m) - m* pow(x, 2) / (1 + m)) * 2 / w;
			double N = Ca*m / (1 + 2 * m*(1 - Ca));
			r_s = ((1 + 2 * N) / 3 / N)*(12 * v_av) / (D_h_s3 - D_po);
			double v_slip = 0.2 * pow(9.8 *(e_cutting - e_den) / e_den, 0.72)* pow(ds, 1.18) / pow(t_w1 / r_s / e_den, 0.45);
		}
		else
		{
			double v_slip = 1.74 * pow(9.8 *(e_cutting - e_den) / e_den, 0.5)* pow(ds, 0.5);
		}
		// 6. calculate the cutting concentration of the wellbore
		double v_mix = v_sl + v_sc;
		double v_mean = (v_mix - v_slip) / 2;
		double C_cuting = v_mean + pow((pow(v_mean, 2) + v_mix*CF / v_slip), 2);
	} 

	// calculate the mix flow density (According to API RD 13)
	//????????????????? transfer unit ???????????????????????????

	if (cutting_concentration == 2)
	{
		// 1. according to the correlation of Sooh
		double e_eff = h *(e_cutting - e_den) / e_den;
		double v_slip = 2.19 * pow(e_eff, 0.5);

		// 2. particle Reynolds number
		double t_s = 7.9 * pow(h*(8.345*e_cutting - e_cutting), 0.5);

		// 3. the apparent viscosity is calculated by this shear rate
		double w = PI*(D_h_s3 + D_po) / 2;
		double h = (D_h_s3 - D_po) / 2;
		double A_anu = w * h;
		double v_av = Q / A_anu;
		double r_s = ((1 + 2 * m) / 3 / m)*(12 * v_av) / (D_h_s3 - D_po);     // shear rate      there is a problem!!!!
		double t_w1 = t_y + K * pow(r_s, m);                                  // initial t_w1
		double e = 1000;
		do
		{
			double x = t_y / t_w1;
			double Ca = 1 - x / (1 + m) - m* pow(x, 2) / (1 + m);
			double t_w2 = t_y + K * pow((r_s / Ca), m);
			e = t_w2 - t_w1;
			t_w1 = t_w2;
		} while (e > 0.00000001);
		double u_app = 479 * t_w1 / r_s;

		// 4. calculate the Reynolds number 
		double N_rep = 928 * e_den * v_slip;
		if (N_rep < 100)
		{
			v_slip = 0.0203*t_s * pow(ds* r_s / pow(e_den, 1 / t_w1), 0.5);
		}
		double A_annulus = PI* (pow(D_h_s3, 2) - pow(D_po, 2)) / 4;
		double v_ann = Q / A_annulus;
		double v_up = v_ann - v_slip;
		double Rt = v_up / v_ann;        // cutting transport ratio
		C_cutting = pow(D_h_s3, 2)* ROP / 448 / Q / Rt;
	}

	   e_mix = e_den * (1 - C_cutting) + e_cutting * C_cutting;	
  // define an array for storing the frictional pressure data.
		vector<double> gravity_pressure_drop;
		 //consider the well trajectory!
		for (int i = 0; i < count; i++)
		{
			double theta = arrayA[i][2];
			double MD = arrayA[i][1] * cos(arrayA[i][2]);
			double pressure_drop_gravity = e_mix * g * MD;    // MD is the step length of the calculation
			gravity_pressure_drop.push_back(pressure_drop_gravity);
		}

	/////////////////////////////////// pressure profile along the wellbore ////////////////////////////////////
	vector<double> Annulus_Pressure_Profile;
	for (int i = 0; i < count; i++)
	{
		Annulus_Pressure_Profile[i] = frictional_pressure_drop[i] + gravity_pressure_drop[i];

	///////////////////////////////////// output the data file/////////////////////////////////////////////////////

		ofstream outfile;
		outfile.open("result the calculation");
		outfile << "Annulus Pressure drop" << Annulus_Pressure_Profile[i] << endl;
		outfile.close();
		return 0;
	}

}
