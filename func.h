#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include "nr.h"
#include <iomanip>
#include <ctime>

double pi = 3.1415;
double epsilon = pow(10,-8);

using namespace std;

vector<double> PressureHRG, BaryonDensityHRG, EntropyDensityHRG;
double muB,m,B,d,Q,S,T,f0,f1,n0, a, b, endpt,endpt_n,sumint,sumint_dens,integral, dens_integral;
int muB_grid_min_MeV, muB_grid_max_MeV, T_grid_min_MeV, T_grid_max_MeV;

// Gaussian Quadrature method for taking integrals -- from Numerical Recipes in C++
DP NR::qgaus(DP func(const DP), const DP a, const DP b)
{
	static const DP x[]={0.1488743389816312,0.4333953941292472,
		0.6794095682990244,0.8650633666889845,0.9739065285171717};
	static const DP w[]={0.2955242247147529,0.2692667193099963,
		0.2190863625159821,0.1494513491505806,0.0666713443086881};
	int j;
	DP xr,xm,dx,s;

	xm=0.5*(b+a);
	xr=0.5*(b-a);
	s=0;
	for (j=0;j<5;j++) {
		dx=xr*x[j];
		s += w[j]*(func(xm+dx)+func(xm-dx));
	}
	return s *= xr;
}

// Particle list variables class
class Plist {
public:

    static double element3, element4, element5, element6, element7, element8, element9, element10, element11, element12;	
    static std::string element1, element2;
	static vector <double> mass;
 	static vector <int> degeneracy;
	static vector <int> baryonnumber;
	static vector <int> strangeness;
	static vector <int> charge;


};


// Function that reads particle list data
void get_particlelist(string particle_list){
	ifstream data(particle_list);

	if(!data.is_open()){
    cerr << "ERROR: could not read particle list" << endl;
    exit(1);
	}	

	while (!data.eof()){

	data >> Plist::element1 >> Plist::element2 >> Plist::element3 >> Plist::element4 >> Plist::element5 >> Plist::element6 >> Plist::element7 >> Plist::element8 >> Plist::element9 >> Plist::element10 >> Plist::element11 >> Plist::element12;
	//Mass is converted from GeV to MeV
		double massmev = 1000*Plist::element3;
		Plist::mass.push_back(massmev);

		Plist::degeneracy.push_back(Plist::element5);
		Plist::baryonnumber.push_back(Plist::element6);
		Plist::strangeness.push_back(Plist::element7);
		Plist::charge.push_back(Plist::element11);

	}

	data.close();
}

//Pressure in Ideal HRG
DP press(const DP x)
{

    return (d*pow(-1,B+1))/(2*pow(pi,2))*pow(x,2)*log(1+pow(-1,B+1)*exp(-sqrt(pow(x,2)+pow(m,2))/T+muB*B/T));
}

//Two entropy density integrals in Ideal HRG

DP entropydensity1(const DP x){

    return pow(x,2)*log(1+ pow((-1),B+1)*exp(-sqrt(pow(x,2)+pow(m,2))/T+muB*B/T));
}

DP entropydensity2(const DP x){


	return pow(x,2)*(sqrt(pow(x,2) + pow(m,2))-B*muB)*exp((B*muB - sqrt(pow(x,2) + pow(m,2)))/T)/(1 + pow(-1,B+1)*exp((B*muB - sqrt(pow(x,2)+pow(m,2)))/T));
}

//Baryon Density in Ideal HRG

DP baryondensity(const DP x){

	return pow(x,2)/(exp(sqrt(pow(x,2)+pow(m,2))/T - muB*B/T)+ pow(-1,B+1));
}

