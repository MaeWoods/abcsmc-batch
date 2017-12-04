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
#include <dirent.h>

using namespace std;

int main (int argc, char * const argv[]) {

  int Job = 1;
  
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
	
	if(Bunk[0]=='J'){

	  Job = atoi(VariableChar);

	}
      }
    }
  }

	
   ////////////////////////////////////////////////////////////
   ///   Loop over the particles                            ///
   ////////////////////////////////////////////////////////////

if(Job==1){

int position;

char Dat[100];
    	int Dn;

	  //this file contains, position and epsilon
      Dn=sprintf(Dat,"position%d.dat",Job);

ifstream infile (Dat);
  int counter = 0; 
  if (infile.is_open()){

    std::string line;
    while(!getline(infile,line).eof()){
      if(line.empty()) continue;

		std::vector<std::string> split;
      std::string buf; 
      stringstream ss(line);
      while (ss >> buf) split.push_back(buf);

		if(counter==0){
			 //NChr
			 position = atoi( split[0].c_str() );
	     }
   
      counter++;
    }
      
  }
  else{
    std::cerr << "Error: open of constant position file unsuccessful: " <<  Dat << std::endl;
  }

      char DRem[100];
      int Dmm;
      Dmm=sprintf(DRem,"rm -r %d/params_curr",position);
	const int dir1= system(DRem);
  
      char DMak[100];
      int Dma;
      Dma=sprintf(DMak,"mkdir -p %d/params_curr",position);
        const int dir2= system(DMak);

      char Dcopy[100];
      int Dcp;
      Dcp=sprintf(Dcopy,"cp %d/params_prev/*.* %d/params_curr",position,position);
        const int dir3= system(Dcopy);

  }

  return(0); 
}
