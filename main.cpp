
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                                                                    ///
/// Evolutionary model ///  Model 1: insertions + deletions Model 2: insertions + deletions + translocations                           ///
///                         Model 3: insertions + deletions + selection Model 4: insertions + deletions + translocations + selection   ///
/// ------------------------------------------------------------------------------------------------------------                       ///
///                                                                                                                                    ///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <string>
#include <cstring>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <random>
#include <sstream>
#include <ctime>
#include "simulation.h"
#include "abcsmc.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

using namespace std;

// to build the program:
//   g++ main.cpp abcsmc.cpp simulation.cpp -o SVModel -L/usr/local/lib/ -lgsl -I/usr/local/include
//
// to run the model:
//   ./SVModel Seed=1 Model=2


//////////////////////////////////////////////////////////
///                                                    ///
///   MAIN                                             ///
///                                                    ///
//////////////////////////////////////////////////////////
int main (int argc, char * const argv[]) {

  int set_seed = 0;
  double Final_ep = 1;
  int Npart = 1;
  int Job = 1;
  double alph = 1;
  int nparam = 1;
  int position = -1;
  
  ////////////////////////////////////////////////////////////
  ///Generating parameters and random seed from bash script///
  ////////////////////////////////////////////////////////////
	
  for (int CommandLineCounter=1; CommandLineCounter < argc; CommandLineCounter++) {
				
    for (int StringPosition=0; StringPosition < strlen(argv[CommandLineCounter]); StringPosition++) {
						
      int EqualPosition = -1;
			
      if (argv[CommandLineCounter][StringPosition]=='=') {

	EqualPosition = StringPosition;				
				
      }
			
      if (EqualPosition > 1) {
	
	int VariableNameLength = EqualPosition;
	
	char VariableNameChar[VariableNameLength];
				
	for (int i=0; i<EqualPosition; i++) {
					
	  VariableNameChar[i] = argv[CommandLineCounter][i];

	}
				
	int VariableStringLength = strlen(argv[CommandLineCounter])-EqualPosition;

	char VariableChar[VariableStringLength];
				
	for (int i=0; i<VariableStringLength; i++) {
					
	  VariableChar[i] = argv[CommandLineCounter][i+1+EqualPosition];
					
	}
								
	double t1,t2; 
	char Bunk[2];
	Bunk[0] = VariableNameChar[0];
	Bunk[1] = VariableNameChar[1];
	
	if(Bunk[0]=='S'){
	  
	  t1 = strtod(VariableChar,NULL);
	  set_seed = int(t1);
	  
	  //std::cout << "Seed " << set_seed << std::endl;					
					
	}else if(Bunk[0]=='P'){
	  
	  Npart = atoi(VariableChar);

	}else if(Bunk[0]=='F'){
	  
	  Final_ep = atof(VariableChar);
	 				
	}else if(Bunk[0]=='N'){
	  
	  nparam = atoi(VariableChar);
	 				
	}else if(Bunk[0]=='J'){

	  Job = atoi(VariableChar);

	}else if(Bunk[0]=='A'){

	  alph = atof(VariableChar);

	}else if(Bunk[0]=='T'){

	  position = atoi(VariableChar);

	}
      }
    }
  }

  std::cout << "Read input arguments" << endl;
  std::cout << "\tNparticles : " << Npart << std::endl;
  std::cout << "\tset_seed  : " << set_seed << std::endl;
  std::cout << "\tFinal_ep      : " << Final_ep << std::endl;
  std::cout << "\talpha       : " << alph << std::endl;
  std::cout << "\tposition       : " << position << std::endl;
	
   ////////////////////////////////////////////////////////////
   ///   Loop over the particles                            ///
   ////////////////////////////////////////////////////////////
  
  abcsmc mabc;
  
  mabc.run_ams(Final_ep, alph, Job, Npart, nparam, set_seed, position, position);
  
  return(0); 
}














