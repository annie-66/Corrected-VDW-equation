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

vector <int> id;
vector <string> name;
vector <double> mass;
vector <int> width;
vector <int> degeneracy;
vector <int> baryonnumber;
vector <int> strangeness;
vector <int> charmness;
vector <int> bottomness;
vector <int> isospin;
vector <int> charge;
vector <int> decay;
vector<double> PressureHRG, BaryonDensityHRG, EntropyDensityHRG;
double muB,m,B,d,Q,S,T,f0,f1,endpt,sumint,integral;
double element3, element4, element5, element6, element7, element8, element9, element10, element11, element12;	
std::string element1, element2;
int muB_grid_size_MeV, T_grid_size_MeV;

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


void get_particlelist(data){
	while (!data.eof()){

	data >> element1 >> element2 >> element3 >> element4 >> element5 >> element6 >> element7 >> element8 >> element9 >> element10 >> element11 >> element12;
	//Mass is converted from GeV to MeV
		double massmev = 1000*element3;
		mass.push_back(massmev);

		degeneracy.push_back(element5);
		baryonnumber.push_back(element6);
		strangeness.push_back(element7);
		charge.push_back(element11);

}
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