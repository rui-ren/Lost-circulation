#include"formation.h"
#include<stdlib.h>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<fstream>
#include<string>
#include<vector>
#include<stdio.h>
#include<math.h>
#include<cstdlib>



int main(){

//  flow rate need to calculate!
	   double e_den = 12.6;    //ppg

	   double Q_loss = 1;   // unit gpm 
	   Q_loss = Q_loss * (3.785 * pow(10, -3) / 60);       // convert gpm to m^3 /s.

	   double t_y_test = 8.4;     // test t_y
	   t_y_test = t_y_test * 4.44822 / 100 / 0.092903;

	   cout << " yield stress of the mud " << t_y_test << endl;

	   double K_test = 0.08;     //lbf/ 100 ft^2 

	   K_test = K_test * 4.44822 / 100 / 0.092903;

	   //cout << " the flow index of the mud " << K_test << endl;

	   double m_test = 0.75;
	   double w_i = ((double)2 /(double) 1000);   //   the initial width of fracture   m                2***mm
	   //cout << w_i << "~~~~~~~~~~~~~~~~~~~~" << endl;
	   double e_den_test = 1000;
	   double  delta_P= 145;    // psi
	   delta_P = delta_P / 145 * 1000000;     // pascal
	  /* cout << delta_P << "~~~~~~~~~~~~~~ pressure" << endl;*/

	   double   rw = 4.6 *25.4 / 1000; // input the wellbore diameter  m. 
	   double   t_max(1000);          // Time to end simulation
	   double   dt(0.1);             // Time step

	   fracture_loss_model(t_max,dt,t_y_test, K_test, e_den, rw, w_i, m_test, delta_P);
}
