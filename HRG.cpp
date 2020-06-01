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
#include "func.h"



using namespace std;



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(){

muB_grid_size_MeV = 600;
T_grid_size_MeV = 820;

//We set up an optional clock to keep track of how long the program takes to run. Uncomment if output is desired.
// int start_s= clock();


// Input particle list
string particle_list = "PDG2016Plus.dat";
ifstream data(particle_list);

if(!data.is_open()){
    cerr << "ERROR: could not read particle list" << endl;
    exit(1);
}

get_particlelist(data);

data.close();


string myfile = "press.dat";
ofstream funcfile(myfile);
    
if(!funcfile.is_open()){
    cerr << "could not open file" << endl;
    exit(1);
}


//The integral loops over 4 variables: Baryon density, temperature, particle, and subdivisions of the integral. 

for(int l=0; l<muB_grid_size_MeV+1; l++){ /// DENSITY LOOP
	muB = l;
	
	for(int j=1; j<T_grid_size_MeV+1; j++){ /// TEMPERATURE LOOP

		double sum = 0;
		T = j;		
	

		for(int k=0; k< mass.size() ;k++){ ///PARTICLE LOOP

			m = mass.at(k);
			d = degeneracy.at(k);
			B = baryonnumber.at(k);
			S = strangeness.at(k);
			Q = charge.at(k);

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