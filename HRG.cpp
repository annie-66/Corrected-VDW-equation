//
//  HRG.cpp
// 
//
//  Created by Debora Mroczek on 8/1/17.
//  Copyright © 2017 Debora Mroczek. All rights reserved.
//


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


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



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



DP press(const DP x)
{

    return (d*pow(-1,B+1))/(2*pow(pi,2))*pow(x,2)*log(1+pow(-1,B+1)*exp(-sqrt(pow(x,2)+pow(m,2))/T+muB*B/T));
}

DP entropydensity1(const DP x){

    return pow(x,2)*log(1+ pow((-1),B+1)*exp(-sqrt(pow(x,2)+pow(m,2))/T+muB*B/T));
}

DP entropydensity2(const DP x){


	return pow(x,2)*(sqrt(pow(x,2) + pow(m,2))-B*muB)*exp((B*muB - sqrt(pow(x,2) + pow(m,2)))/T)/(1 + pow(-1,B+1)*exp((B*muB - sqrt(pow(x,2)+pow(m,2)))/T));
}

DP baryondensity(const DP x){

	return pow(x,2)/(exp(sqrt(pow(x,2)+pow(m,2))/T - muB*B/T)+ pow(-1,B+1));
}




int main(){
int start_s= clock();
cout << "* Importing PDG data" << endl;
cout << "..." << endl;

string filename = "PDG2016Plus.dat";
ifstream data(filename);

if(!data.is_open()){
    cerr << "could not open file" << endl;
    exit(1);
}

while (!data.eof()){

data >> element1 >> element2 >> element3 >> element4 >> element5 >> element6 >> element7 >> element8 >> element9 >> element10 >> element11 >> element12;
	
	double massmev = 1000*element3;
	mass.push_back(massmev);
	degeneracy.push_back(element5);
	baryonnumber.push_back(element6);
	strangeness.push_back(element7);
	charge.push_back(element11);

}

data.close();

string myfile = "press.dat";
ofstream funcfile(myfile);
    
if(!funcfile.is_open()){
    cerr << "could not open file" << endl;
    exit(1);
}




// // CHECK
// 	m = mass.at(10);
// 	d = degeneracy.at(10);
// 	B = baryonnumber.at(10);
// 	S = strangeness.at(10);
// 	Q = charge.at(10);
// 	muB = 100;
// 	T=150;

// 	cout << "Mass: " << m << endl;
// 	cout << "B: " << B << endl;
// 	cout << "d: " << d << endl;

	// for(int h= 0; h<2000;h++){

	// 	double yy = (entropydensity2(h))/pow(T,3);
	// 	funcfile << h << " " << yy << endl;
	// 	}





for(int l=0; l<601; l++){
	muB = l;
	cout << "muB = " << muB << "\r" << endl; 	/// DENSITY LOOP
for(int j=1; j<821; j++){

	double sum = 0;
	T = j;

	// cout << "Calculating HRG pressure for T = " << T << endl;			/// TEMPERATURE LOOP
	

for(int k=0; k< mass.size() ;k++){										///PARTICLE LOOP

	m = mass.at(k);
	d = degeneracy.at(k);
	B = baryonnumber.at(k);
	S = strangeness.at(k);
	Q = charge.at(k);

	for(int h= 10; h<100;h++){

		endpt = h;
		f0 = press(h);

		if(f0<epsilon){
			endpt = h;
			break;
		}	
	}

	sumint = 0.0;											///INTEGRAL LOOP
	double X2 = 0.0;
		for(int p = 0; p<10;p++){

			double X1 = X2;	
			X2 = X2 + endpt/10.0;									

        	double ss;
        	ss=NR::qgaus(press,X1,X2);
 			sumint = sumint+ss;
     	
	
		integral = (sumint/pow(T,3));
	}
	
	sum = sum + integral;
	}
		  
    funcfile << muB << " " << T << " " << sum << endl;        
 }

}

funcfile.close();
int stop_s=clock();

cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << endl;


return 0;	
}