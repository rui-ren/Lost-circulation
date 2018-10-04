#pragma once
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <string>

using namespace std;
double PI = 3.1415926;

void calculate(vector< double> & time, vector< double> & r_fd, vector <double> & r_frac, vector<double> & Q_loss, vector<double> & Volume, vector <double> & pressure_drop, double e_den, const double t_max, const double dt, double m, double rw, double w_i, double K, double delta_P, double t_y)
{
	double iterations(t_max / dt);
	//cout << "iteration number*************************" << iterations << endl;
	vector <double> pressure_drop1 = { {0} };
	
	vector<double> leak_off = { { 0 } };
	double w_i_initial = w_i;

	for (int i(1); i < iterations; i++)
	{
		time.push_back(time[i - 1] + dt);
		if (delta_P > 0)
		{
			// Part NO.1                                                   Calculate the mud penetration lenght along the wellbore
			double alpha = (2 * m + 1) / (m + 1) * (2 * rw / w_i)*(t_y / delta_P);
			double beta = (m / (2 * m + 1)) * pow((w_i / rw), ((1 + m) / m)) * pow((delta_P / K), (1 / m));
			double c = beta * pow((1 - alpha*(r_fd[i - 1] - 1)), (1 / m)) / pow(2, ((m + 1) / m)) / r_fd[i - 1] / pow(((pow(r_fd[i - 1], (1 - m)) - 1) / (1 - m)), (1 / m));

			r_fd.push_back((c * dt + r_fd[i - 1]));         // the invasion length (dimensionless)
			//	cout << pow(2, ((m + 1) / m)) <<"    ~~~~~~~~~~~~~~~~   " << endl;
			//	cout << r_fd[i - 1]<< "!!!!!!!!!!!!!!!!!!!!" << endl;
			//	cout << pow(((pow(r_fd[i - 1], (1 - m)) - 1) / (1 - m)), (1 / m)) << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
			//cout << "alpha    " << alpha << endl;
			//	cout << "beta     " << beta << endl;
			//	cout << " c       " << c << endl;    // really small number!!!
			//	cout << "r_fd     " << r_fd[i-1] << endl;
			/*cout << "invasion radius " << r_fd[i] << endl;*/
			///*	cout << r_fd[i] * rw <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< endl;*/
			r_frac.push_back(r_fd[i] * rw);

			//cout << "invasion radius " << r_frac[i] << endl;
			//system("pause");

			// calculate the mud loss flow rate !!!
			Q_loss.push_back(2 * PI * w_i * r_frac[i] * (r_frac[i] - r_frac[i - 1]) / dt);          //    velocity    r_w_i  is the radius of wellbore!

			// calculate the mud loss total volume !!!

			Volume.push_back(Volume[i - 1] + PI * w_i * (pow(r_frac[i], 2) - pow(r_frac[i - 1], 2)));
			//  w_i is changable !!!!

			// Part 1. 
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ The pressure drop along the fracture ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			/*double v_av_frac = Q_loss[i - 1] / 2 / PI / r_w / w_i;	*/			    // average velocity   the Radius is changing !!! be ware of that!!!

			double v_av_frac = Q_loss[i] / PI / 2 / r_frac[i] / w_i;
			double r_s_fra = ((1 + 2 * m) / 4 / m)*(8 * v_av_frac) / w_i;        // shear rate  at fracture  assume that m = N
			double t_w1 = t_y + K * pow(r_s_fra, m);								    // initial t_w1
			/*std::cout << "wall shear stress is: ~~~~~~~~~~" << t_w1 << endl;*/
			double e = 1;

			// calculate the wall shear stress
			do
			{
				double x = t_y / t_w1;
				double Ca = (1 - x / (1 + m) - m* pow(x, 2) / (1 + m)) * 2;
				double t_w2 = t_y + K * pow((r_s_fra / Ca), m);
				e = fabs(t_w2 - t_w1);
				t_w1 = t_w2;
			} while (e > 0.000000001);
			// Calculate Generalized flow-behavior index. 

			/*std::cout << "wall shear stress is: " << t_w1 << endl;*/

			// .............................................................. Calculate the velocity profile of fracture .................................................

			double yp = w_i / 2 * t_y / t_w1;            // plug zone area.

			/*cout << yp << "!!!!!!!!!!" << endl;*/
			double vp = w_i / 2 / t_w1 / pow(K, (1 / m))* m / (1 + m)* pow((t_w1 - t_y), ((1 + m) / m));
			/*cout << vp << "plug zone velocity of the fracture" << endl << '\n';*/


			for (double y = yp; y < (w_i / 2); y += ((w_i - yp) / 20))
			{
				double vv = w_i / 2 / t_w1 / pow(K, (1 / m))* m / (1 + m)* (pow((t_w1 - t_y), ((1 + m) / m)) - pow((2 * t_w1 / w_i *y - t_y), ((1 + m) / m)));
				/*cout <<y * 1000 <<  " " <<  vv  << endl;*/
			}
			
			double x = t_y / t_w1;          // the eventually data from calculation.
			double Ca = (1 - x / (1 + m) - m* pow(x, 2) / (1 + m));
			double N = Ca*m / (1 + 2 * m*(1 - Ca));                                              // Generalized N , Calculate !!
			/*	std::cout << "the generalized fluid flow index(0.15 < N < 0.4): " << N << endl;*/

			//// Reynold's number of the calculation !!!!!!!!!!!!!!!!!!!!!!!!!!!!     check the Reynold's number, to see if we can get Turbulent Flow.

			r_s_fra = ((1 + 2 * N) / 2 / N)*(8 * v_av_frac) / w_i;
			double u_app = t_w1 / r_s_fra;
			// calculate the Reynolds number            according to this penny shape, and calculate the Reynold's number
			double D_hy = 4 * PI * r_frac[i] * w_i / (w_i + 2 * PI * r_frac[i]);
			double Re = e_den * v_av_frac * D_hy / u_app;

			// check the Reynold's number !!!!!
			/*std::cout << " !!!!! the Reynold's number for the fracture is !!!!!" << Re << endl;*/

			double f1 = 0;
			// calculate the frictional factor of the fracture
			if (Re < 2000)
			{
				f1 = 16 / Re;
			}
			else
			{
				// how to calculate the friction factor with Colebrook correlation.
				double rough = 0.0003;   // roughness  degree !!!
				double D_eq = w_i;       // equivalent  diameter 
				double f_1 = 0.001;
				double abc = 100;
				do
				{
					double right = 1.14 - 2 * log(0.000001 + 9.34 / (Re*pow(f_1, 0.5)));              // A trial and error solution was used to solve colebrook equation with finning friction factor
					double f2 = pow(1 / right, 2);
					abc = fabs(f2 - f_1);
					f_1 = f2;
				} while (abc > 1e-6);
			}
			// the pressure drop along the fracture can be calculated by
			double pressure_drop_fracture = 2 * t_w1 / w_i * (r_frac[i] - r_frac[i - 1]);

			if (i == 1)
			{
				pressure_drop1.push_back(241.85);   //pascal
			}
			else
			{
				pressure_drop1.push_back(pressure_drop_fracture);
			}
		
			delta_P = delta_P - pressure_drop1[i];                       // the new pressure value in the fracture

			pressure_drop.push_back( delta_P );
			

			/*cout << pressure_drop[i] << " the pressure distribution in the fracture" << endl;*/

			//// Part 2.  (according to Carter. 1969)
			////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ The fracture width change along the fracture ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// take the leak off into account!!!! 

			//resolvtion integration
			double C_leak = 3 * pow(10, -6);        // unit ... in / sqrt(s)
			C_leak = C_leak * 25.4 / 1000;      // unit ... m / sqrt(s)

			double z = 0;                       // store the leak_off volume.
			double r = 0;
			double o = 0;

			if (i < 11)
			{
				for (int j = (i - 1); j > 0; j--)
				{
					double k = PI * (pow(r_frac[j + 1], 2) - pow(r_frac[j], 2)) * 2 * C_leak / sqrt((i - j) * dt) * dt;     // this is lost rate !!!
					z = z + k;
				}
			}
			else
			{
				for (int j = (i - 1); j > (i - 11); j--)
				{
					double k = PI * (pow(r_frac[j + 1], 2) - pow(r_frac[j], 2)) * 2 * C_leak / sqrt((i - j) * dt) * dt;   // this is the loss volume!!

					? ? ? ? ? ? lost volume ... possibly it this is worng  !!!!!
					z = z + k;
				}
			}

			/*cout << z << "zzzzzzzzzzzzzzzzzzzzzzzzzzzz" << endl;*/
			leak_off.push_back(leak_off[i - 1] + z);        // the leak off volume

			Volume[i] = Volume[i] + leak_off[i - 1];
			Q_loss[i] = Q_loss[i] - r;

			//cout << Q_loss[i] << " the total loss volume of the mud" << endl;

			//// Part 3.  (according to Lavrov. 2012)
			////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ The pressure leak-off along the fracture ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			double K_n = 3 * 100000 * 6894.75729 / 25.4 * 1000;        // stiffness of the fracture.    3e5 psi / in
			w_i = w_i_initial + delta_P / K_n;
		}
		else
		{
			r_fd.push_back ( r_fd[i-1]);
			r_frac.push_back( r_frac[i-1]);
			Volume.push_back( Volume[i-1]);
			Q_loss.push_back(Q_loss[i-1]);
			pressure_drop.push_back(delta_P);
		}
	}
}


void save(const vector <double> & time, const vector <double> & r_frac, const vector <double> & Q_loss, const vector <double> & Volume, vector <double>& pressure_drop,const string & filename)
{
	ofstream file_out;
	file_out.open(filename);
	for(int i(0); i < time.size(); i++)
	{
		/*file_out << (time[i]/60) << "  " << r_frac[i] <<"  " << Q_loss[i]<< "  " <<  Volume [i]<< endl;*/
		file_out << (time[i]/60) << "  " << pressure_drop[i] << endl; 
	}
	file_out.close();
}


double fracture_loss_model(double t_max,double dt,double t_y, double K, double e_den, double rw, double w_i, double m, double delta_P)
{	
	/////////////////////////////////////////////////*** calculate the pressure distribution in the fracture***///////////////////////////////////////////

	vector<double> time{ { 0 } };        // initial time = 0
	vector<double> r_fd{ { 50} };      // initial r_fd = rw
	vector<double> r_frac{ { 50 * rw } };
	double a = PI * pow((49 * rw), 2) * w_i;
	vector<double> Q_loss{ { 0 } };      // initial Q_loss       this can be known from the wellhead!
	vector<double> Volume{ { a } };
	vector <double> pressure_drop = { { delta_P } };
	calculate(time, r_fd, r_frac, Q_loss, Volume, pressure_drop, e_den,t_max, dt, m, rw, w_i, K, delta_P, t_y);
	save(time, r_frac, Q_loss, Volume,pressure_drop, "pressure.xls");

	//// test code
	//cout << Volume[0] << endl;
	//cout << Volume[1] << endl;
	//cout << Volume[2] << endl;
	//cout << r_frac[2] << "~~~~~~~~~"<< endl;
	//cout << r_frac[1]<<"!!!!!!!!!!!" << endl;
	//cout << PI * (pow(r_frac[t_max], 2) - pow(r_frac[t_max - 1], 2))* w_i << "eeeeeeeeee"<< endl;
	//cout << Q_loss[t_max] << endl;
	cout << Volume[double (t_max)] << endl;

	system("pause");
	return 0;
}
