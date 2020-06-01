//
//  HRG.cpp
// 
//
//  Created by Debora Mroczek on 8/1/17.
// 	Updated on 06/01/2020.
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
#include <iomanip>
#include <ctime>
#include <iostream> 
#include "func.h"
#include "nr.h"



using namespace std;


double Plist::element3, Plist::element4, Plist::element5, Plist::element6, Plist::element7, Plist::element8, Plist::element9, Plist::element10, Plist::element11, Plist::element12;	
std::string Plist::element1, Plist::element2;
vector <double> Plist::mass;
vector <int> Plist::baryonnumber;
vector <int> Plist::degeneracy;
vector <int> Plist::strangeness;
vector <int> Plist::charge;



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(){

muB_grid_min_MeV = 0;
muB_grid_max_MeV = 600;

T_grid_min_MeV = 1;
T_grid_max_MeV = 820;

//We set up an optional clock to keep track of how long the program takes to run. Uncomment if output is desired.
// int start_s= clock();


// Input particle list
string particle_list = "PDG2016Plus.dat";

//Create substring with particle list name
// size_t pos = particle_list.find(".dat");
// string foldername = particle_list.substr(0,pos);



get_particlelist(particle_list);


//Create output file name 

string myfile = "PRESS_HRG_IDEAL_MUB_" + to_string(muB_grid_min_MeV) + "_" + to_string(muB_grid_max_MeV) + "_T_" + to_string(T_grid_min_MeV) + "_" + to_string(T_grid_max_MeV) + ".dat";

ofstream funcfile(myfile);
    
if(!funcfile.is_open()){
    cerr << "could not open file" << endl;
    exit(1);
}


//The integral loops over 4 variables: Baryon density, temperature, particle, and subdivisions of the integral. 

for(int l= muB_grid_min_MeV; l<muB_grid_max_MeV+1; l++){ /// DENSITY LOOP
	muB = l;
	
	for(int j=T_grid_min_MeV; j<T_grid_max_MeV+1; j++){ /// TEMPERATURE LOOP

		double sum = 0;
		T = j;		
	

		for(int k=0; k<  Plist::mass.size() ;k++){ ///PARTICLE LOOP

			m = Plist::mass.at(k);
			d =  Plist::degeneracy.at(k);
			B =  Plist::baryonnumber.at(k);
			S =  Plist::strangeness.at(k);
			Q =  Plist::charge.at(k);

			//Integral contribution check
			for(int h= 1; h<100;h++){

				endpt = h;
				f0 = press(h);

				// Checks where the integral stops contributing. This is the end point of the function. 
				if(f0<epsilon){
				endpt = h;
				break;
				}	
			}

			/// INTEGRAL LOOP

			//Must reset the sum after each particle.
			sumint = 0.0;	
			//Each integral is broken into 10 smaller integrals. 
			double X2 = 0.0;
			for(int p = 0; p<10;p++){

				double X1 = X2;	
				X2 = X2 + endpt/10.0;									

        		double ss;
        		ss=NR::qgaus(press,X1,X2);
 				sumint = sumint+ss;
     	
	
				integral = (sumint/pow(T,3)); //Integral normalized by the correct power of temperature.
			}
	
			sum = sum + integral; // updates the total sum after each particle.
		}
		  
    //Writes the output to the file 
    funcfile << muB << " " << T << " " << sum << endl;        
 	}

}

funcfile.close();

// Uncomment for runtime output
// int stop_s=clock();
// cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << endl;


return 0;	
}